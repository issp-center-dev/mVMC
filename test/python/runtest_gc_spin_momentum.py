"""Deterministic combined spin/momentum replay, based on the Phase 3b B1 oracle.

No statistical thresholds: rational Pfaffians use all dumped binary64 QP
coefficients, independently checked against Fock rotations and frozen inputs.
"""
import argparse
import cmath
import json
import math
import os
from types import SimpleNamespace

import numpy as np
import grandcanonical_ap_oracle as ap
import grandcanonical_spin_oracle as gc
import runtest_gc_mp as mp
import runtest_gc_spin as spin
from gc_antiparallel_failure import require_expected_failure
from runtest_gc_antiparallel import parse_sr_dump

def addq(first, second):
    return gc.qadd(first, second)

def combined_wave(F, dump, mappings, log_charge):
    """Exact-Q+iQ oracle using the dumped combined binary64 QP weights."""
    gauss = dump["gauss"]
    base = gc.antiparallel_matrix(F)
    translated = []
    for r in range(dump["count"] // gauss):
        matrix = ap.translated_matrix(base, mappings, r, spin.NSITE)
        translated.append([[matrix[i][j + spin.NSITE]
                            for j in range(spin.NSITE)]
                           for i in range(spin.NSITE)])
    exact = dict((mask, gc.ZERO) for mask in gc.balanced_basis(spin.NSITE))
    for node in dump["nodes"]:
        r = node["q"] // gauss
        component = gc.rational_pfaffians(
            gc.rotated_block_matrix(translated[r], node))
        weight = gc.rational(node["weight"])
        for mask in exact:
            exact[mask] = addq(exact[mask], gc.qmul(weight, component[mask]))
    raw = dict((mask, 0.0j) for mask in gc.even_masks(spin.NSITE))
    charged = dict(raw)
    for mask, value in exact.items():
        if value != gc.ZERO:
            raw[mask] = gc.qfloat(value)
            charged[mask] = raw[mask] * math.exp(log_charge(mask))
    return exact, raw, charged

def validate_dump(dump, weights):
    gauss = dump["gauss"]
    if dump["count"] != gauss * len(weights):
        raise AssertionError("combined QP count differs")
    reference = []
    worst = 0.0
    for r, trans_weight in enumerate(weights):
        if trans_weight == 0:
            raise AssertionError("zero translation weight")
        group = []
        for node in dump["nodes"][r * gauss:(r + 1) * gauss]:
            item = dict(node)
            item["weight"] = complex(node["weight"]) / trans_weight
            group.append(item)
        worst = max(worst, gc.validate_coefficients(group, gauss))
        if r == 0:
            reference = group
        else:
            for actual, expected in zip(group, reference):
                for key in ("cs", "cc", "ss", "weight"):
                    worst = max(worst, abs(actual[key] - expected[key]))
    if worst > gc.COEFFICIENT_TOLERANCE:
        raise AssertionError("combined coefficient mismatch {}".format(worst))
    return worst

def mean(masks, values):
    return sum(values[mask] for mask in masks) / float(len(masks))

def sr_row(masks, values, energy_values, energy):
    mean_o = mean(masks, values)
    mean_oh = sum(energy_values[mask] * values[mask]
                  for mask in masks) / float(len(masks))
    gradient = 2.0 * (mean_oh.real - energy.real * mean_o.real)
    return (mean_o.real, mean_o.imag, mean_oh.real, mean_oh.imag, gradient)

def expected_outputs(fixture, F, dump, mappings, masks, calmode):
    exact, raw, wave = combined_wave(F, dump, mappings, fixture.log_charge)
    if any(wave[mask] == 0 for mask in masks):
        raise AssertionError("production sampled an exact combined node")
    unique = sorted(set(masks))
    tables = spin.local_tables(wave, fixture.terms(), unique, calmode)
    energy = mean(masks, tables["energy"])
    energy2 = sum(abs(tables["energy"][mask]) ** 2 for mask in masks) / \
        float(len(masks))
    numbers = dict((mask, spin.popcount(mask)) for mask in unique)
    number = mean(masks, numbers)
    number2 = sum(numbers[mask] ** 2 for mask in masks) / float(len(masks))
    expected = {
        "energy": energy,
        "energy2": energy2,
        "relative_variance": ((energy2 - energy * energy) /
                              (energy * energy)).real,
        "number": number,
        "number2": number2,
        "variance_number": number2 - number * number,
        "sz": 0.0,
        "sz2": 0.0,
        "s2": mean(masks, tables["s2"]),
    }
    for group in ("onebody", "twobody", "nbody", "anomalous"):
        if group in tables:
            expected[group] = dict((key, mean(masks, values))
                                   for key, values in tables[group].items())
    if calmode == 0:
        expected["sr"] = expected_sr(
            fixture, F, dump, mappings, masks, wave, tables["energy"], energy)
    return exact, raw, wave, expected

def translated_f(F, mappings, r):
    matrix = ap.translated_matrix(gc.antiparallel_matrix(F), mappings, r,
                                  spin.NSITE)
    return [row[spin.NSITE:] for row in matrix[:spin.NSITE]]


def fock_evaluator(dump, mappings, weights):
    """Independent creator-rotation operators, cached for finite differences."""
    groups = []
    for r, weight in enumerate(weights):
        nodes = [dict(node, weight=node['weight'] / weight)
                 for node in dump['nodes'][r * dump['gauss']:(r + 1) * dump['gauss']]]
        basis, operator = gc.rotation_operator(spin.NSITE, nodes)
        groups.append((weight, operator))

    def evaluate(F):
        total = np.zeros(len(basis), dtype=complex)
        for r, (weight, operator) in enumerate(groups):
            pair = gc.pair_wave(translated_f(F, mappings, r))
            total += weight * (operator @ np.array([pair[mask] for mask in basis]))
        return dict(zip(basis, total))
    return evaluate


def expected_sr(fixture, F, dump, mappings, masks, wave, energy_values, energy):
    unique = sorted(set(masks))
    result = {}
    projection = mp.projection_derivatives(fixture.mp, unique)
    for parameter in range(spin.NPROJ):
        for imaginary in (0, 1):
            values = {mask: complex(projection[parameter][mask])
                      if not imaginary else 0.0j for mask in unique}
            result[2 * parameter + imaginary] = sr_row(masks, values, energy_values, energy)
    evaluate = fixture.combined_evaluate
    base = evaluate(F)
    for parameter in range(spin.NSITE ** 2):
        for imaginary in (0, 1):
            estimates = []
            for epsilon in (1.0e-5, 5.0e-6):
                plus, minus = [list(row) for row in F], [list(row) for row in F]
                i, j = divmod(parameter, spin.NSITE)
                step = epsilon * (1j if imaginary else 1)
                plus[i][j] += step
                minus[i][j] -= step
                positive, negative = evaluate(plus), evaluate(minus)
                estimates.append({m: (positive[m] - negative[m]) / (2 * epsilon)
                                  for m in unique})
            if any(abs(estimates[0][m] - estimates[1][m]) >
                   1e-6 * (1 + abs(estimates[1][m])) for m in unique):
                raise AssertionError('unstable combined finite difference')
            values = {m: (4 * estimates[1][m] - estimates[0][m]) / (3 * base[m])
                      for m in unique}
            result[2 * (spin.NPROJ + parameter) + imaginary] = sr_row(
                masks, values, energy_values, energy)
    return result


def check_frozen(fixture, k_index, F, raw):
    path = os.path.join(os.path.dirname(__file__), 'data', 'gc_spin_momentum_fixture.json')
    with open(path) as stream:
        frozen = json.load(stream)
    encoded = [[[z.real.hex(), z.imag.hex()] for z in row] for row in F]
    if encoded != frozen['fixture_F']:
        raise AssertionError('combined fixture differs from pre-implementation freeze')
    case = next(c for c in frozen['cases'] if c['boundary'] == fixture.boundary and
                c['G'] == fixture.gauss and c['k_index'] == k_index)
    for text_mask, pair in case['raw_amplitude_binary64_display'].items():
        reference = complex(*(float.fromhex(word) for word in pair))
        spin.compare('frozen combined mask ' + text_mask, raw[int(text_mask)],
                     reference, gc.FOCK_TOLERANCE * (1 + abs(reference)))


def run_case(rootdir, boundary, gauss, k_index, calmode, store, procs=1,
             threads=1, namespace="replay", sample_count=spin.REPLAY_SAMPLES):
    fixture = spin.SpinFixture(boundary, gauss)
    label = 'gc_smp_{}_g{}_k{}_m{}_o{}_np{}_t{}'.format(
        boundary, gauss, k_index, calmode, store, procs, threads) + '_' + namespace
    workdir = spin.prepare_work(rootdir, label)
    spin.write_fixture(workdir, fixture, calmode=calmode, store=store,
                       nsplit=procs, samples=sample_count,
                       seed=spin.REPLAY_SEED, iterations=2, freeze=True,
                       nmptrans=-4 if boundary == 'apbc' else 4, momentum=True)
    mp.write_translations(workdir, fixture.mp.lattice, k_index)
    environment = {'MVMC_GC_STATE_DUMP': 'state.dat', 'MVMC_GC_STATE_AMPLITUDE': '1',
                   'OMPI_MCA_mpi_yield_when_idle': '1'}
    if calmode == 0:
        environment['MVMC_GC_SR_DUMP'] = 'sr.dat'
    if threads > 1:
        environment.update(OMP_NUM_THREADS=str(threads), OMP_DYNAMIC='FALSE',
                           MVMC_GC_OMP_AUDIT='1')
    output = spin.run_with_dump(rootdir, workdir, procs, environment)
    dump = gc.read_qp_dump(os.path.join(workdir, spin.DUMP_NAME), allow_complex_weights=True)
    weights, mappings = mp.projection_data(fixture, workdir)
    validate_dump(dump, weights)
    F = spin.read_pair_matrix(workdir)
    spin.check_readback(F, fixture.parameters())
    exact, raw, wave = combined_wave(F, dump, mappings, fixture.log_charge)
    check_frozen(fixture, k_index, F, raw)
    fixture.combined_evaluate = fock_evaluator(dump, mappings, weights)
    independent = fixture.combined_evaluate(F)
    for mask in exact:
        spin.compare('combined Fock mask {}'.format(mask), raw[mask], independent[mask],
                     gc.FOCK_TOLERANCE * (1 + abs(raw[mask])))
    if k_index == 1 and (exact[0] != gc.ZERO or exact[255] != gc.ZERO):
        raise AssertionError('boundary reference is not exact zero')
    per_rank = []
    for rank in range(procs):
        records = spin.state_dump_records(os.path.join(workdir, 'state.dat.rank{}'.format(rank)))
        for record in records:
            if record[0] != 'SAMPLE':
                continue
            pieces = record[4]
            if len(pieces) != 5 or pieces[4].split()[0] != 'LOGIP':
                raise AssertionError('missing sample amplitude')
            words = pieces[0].split()
            order = [int(word) for word in words[3:3 + record[2]]]
            inversions = sum(order[i] > order[j] for i in range(len(order))
                             for j in range(i + 1, len(order)))
            sign = -1 if inversions % 2 else 1
            logip = complex(*(float.fromhex(word) for word in pieces[4].split()[1:]))
            actual = cmath.exp(logip)
            expected = sign * raw[record[3]]
            spin.compare('sampler amplitude rank{} sample{}'.format(rank, record[1]),
                         actual, expected, 5e-12 * (1 + abs(expected)))
        per_rank.append(mp.replay_iterations(records, samples=sample_count))
    for iterations in per_rank[1:]:
        if iterations != per_rank[0]:
            raise AssertionError('QP split ranks disagree on configuration/order/burn state')
    if calmode == 0:
        mp.assert_frozen_parameters(workdir, fixture.parameters())
    all_masks = []
    for iteration, records in enumerate(per_rank[0]):
        masks = [r[3] for r in records]
        all_masks.extend(masks)
        if any(m not in exact or exact[m] == gc.ZERO for m in masks):
            raise AssertionError('sample outside combined support')
        unused, unused_raw, unused_wave, expected = expected_outputs(
            fixture, F, dump, mappings, masks, calmode)
        spin.compare_replay(label, workdir, expected, calmode, iteration, len(fixture.parameters()))
        if iteration == 1:
            mp.check_iteration_one_guards(workdir, expected, calmode, len(fixture.parameters()))
    # This dense fixture need not enter the singular-component cofactor path.
    # Audit actual parallel Green work in measurement mode; optimization mode
    # still compares every SR row and sample against the OMP=1 reference.
    # Dedicated GC_SP/GC_MP singular fixtures audit parallel cofactor work.
    if threads > 1 and calmode == 1:
        mp.assert_omp_audit(output, threads, 'green')
    print('{}: all {} QP, {} samples/rank, physical/Green/SR PASS'.format(
        label, dump['count'], len(all_masks)))
    sr = parse_sr_dump(os.path.join(workdir, 'sr.dat')) if calmode == 0 else None
    return all_masks, sr


def check_store(first, second):
    if first[0] != second[0]:
        raise AssertionError('NStoreO changes combined sample chain')
    for a, b in zip(first[1], second[1]):
        if set(a['p']) != set(b['p']) or len(a['p']) != 2 * (spin.NPROJ + spin.NSITE ** 2):
            raise AssertionError('combined SR row inventory')
        for key in a['p']:
            for x, y in zip(a['p'][key], b['p'][key]):
                spin.compare('combined NStoreO row {}'.format(key), x, y, 5e-10 * (1 + abs(x)))


def negative_cases(rootdir, procs):
    # Reader order, malformed dimensions and remaining gates in combined mode.
    spin.INPUT_CASES['combined_workspace'] = (
        {'gauss_token': '200000000'}, 'NSPGaussLeg workspace')
    spin.INPUT_CASES['combined_r_over'] = (
        {'momentum': True, 'nmptrans': -2147483648},
        'NMPTrans must be a nonzero integer')
    for order in ('early', 'late'):
        for name in ('general', 's_nonzero', 'spin_jastrow', 'opttrans',
                     'g_zero', 'g_fraction', 'g_workspace_sum_over',
                     'combined_workspace', 'combined_r_over'):
            options, message = spin.INPUT_CASES[name]
            combined_name = 'smp_' + name
            combined_options = dict(options, momentum=True, nmptrans=4)
            if name == 'combined_r_over':
                combined_options = options
            spin.INPUT_CASES[combined_name] = (combined_options, message)
            try:
                spin.input_case(rootdir, SimpleNamespace(case=combined_name, order=order, np=procs))
            finally:
                del spin.INPUT_CASES[combined_name]
    for boundary in ('pbc', 'apbc'):
        for n in (0, 8):
            fixture = spin.SpinFixture(boundary, 8)
            workdir = spin.prepare_work(rootdir, 'gc_smp_init_{}_n{}_np{}'.format(boundary, n, procs))
            spin.write_fixture(workdir, fixture, calmode=1, store=1, nsplit=procs,
                               samples=8, seed=11, iterations=1, freeze=True,
                               nmptrans=-4 if boundary == 'apbc' else 4, momentum=True)
            mp.write_translations(workdir, fixture.mp.lattice, 1)
            path = os.path.join(workdir, 'modpara.def')
            with open(path) as stream:
                content = stream.read()
            spin.write(path, content.replace('NGCInitNelec   4', 'NGCInitNelec   {}'.format(n)))
            code, output = spin.run_vmc(rootdir, workdir, procs, [])
            require_expected_failure(code, output, ['zero projected amplitude',
                                                   'after 1 attempt(s)', 'Ncur={}'.format(n)])
    print('combined negative/initial boundary np{} PASS'.format(procs))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('mode', choices=('replay', 'mpi', 'negative', 'omp'))
    parser.add_argument('--boundary', choices=('pbc', 'apbc'), default='pbc')
    parser.add_argument('--k-index', type=int, choices=(0, 1), default=1)
    parser.add_argument('--np', type=int, default=2)
    args = parser.parse_args()
    rootdir = os.getcwd()
    if args.mode == 'negative':
        negative_cases(rootdir, 1)
        if args.np > 1:
            negative_cases(rootdir, args.np)
    elif args.mode == 'replay':
        for gauss in (2, 8):
            run_case(rootdir, args.boundary, gauss, args.k_index, 1, 0)
            first = run_case(rootdir, args.boundary, gauss, args.k_index, 0, 0)
            second = run_case(rootdir, args.boundary, gauss, args.k_index, 0, 1)
            check_store(first, second)
    else:
        procs = args.np if args.mode == 'mpi' else 1
        threads = 4 if args.mode == 'omp' else 1
        gauss = 2 if procs > 3 else 8
        # Empty-rank coverage is deterministic; a shorter chain avoids MPI
        # oversubscription dominating the suite while retaining all output rows.
        sample_count = 64 if procs > 3 else spin.REPLAY_SAMPLES
        for calmode in (1, 0):
            namespace = args.mode + str(procs)
            reference = run_case(rootdir, args.boundary, gauss, args.k_index,
                                 calmode, 1, namespace=namespace + '_reference',
                                 sample_count=sample_count)
            actual = run_case(rootdir, args.boundary, gauss, args.k_index,
                              calmode, 1, procs, threads, namespace=namespace,
                              sample_count=sample_count)
            if reference[0] != actual[0]:
                raise AssertionError('serial/MPI/OMP sample chains differ')
            if calmode == 0:
                check_store(reference, actual)


if __name__ == '__main__':
    main()
