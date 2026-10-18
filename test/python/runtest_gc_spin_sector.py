"""Regression for disconnected S=0 sectors, including endpoint observables.

Reference amplitudes come from explicit Fock-space spin rotation, never the
production Pfaffian. Only the proven triplet parity (and translation boundary)
zeros are imposed; near-triplet amplitudes are retained. Statistical checks use
block means, while output replay on the actual samples is deterministic.
"""
import argparse
from collections import Counter
import json
import math
from pathlib import Path
import os
import subprocess
import tempfile
import time

import numpy as np
import grandcanonical_spin_oracle as spin
from grandcanonical_exact_oracle import one_body
from runtest_gc import mpi_command, state_dump_records
from runtest_gc_spin import local_of


def table(name, count, rows, complex_type=False):
    return '\n'.join(['====================', '{} {}'.format(name, count),
                      'ComplexType 1' if complex_type else '====================',
                      '====================', '===================='] + rows) + '\n'


def fixture(kind):
    if kind in ('minimal', 'near', 'singlet'):
        eps = 1e-4 if kind == 'near' else 0
        F = np.array([[eps, 1], [-1, eps]], dtype=complex)
        if kind == 'singlet':
            F = np.array([[0.25, 1], [1, 0.5]], dtype=complex)
    else:
        size = 6 if kind == 'six' else 4
        rng = np.random.default_rng(7)
        raw = rng.normal(size=(size, size)) + 1j*rng.normal(size=(size, size))
        F = (raw-raw.T)/4
        if kind == 'general': F += np.eye(size)*0.2
    return F


def write_input(work, F, gauss, initial, seed, samples, np_, combined=False,
                apbc=False):
    work.mkdir(parents=True, exist_ok=False)
    L = len(F)
    options = dict(NVMCCalMode=1, NLanczosMode=0, NDataIdxStart=1,
                   NDataQtySmp=1, Nsite=L, Ne=L//2, NSPGaussLeg=gauss,
                   NSPStot=0, NMPTrans=(-L if apbc else L) if combined else 1,
                   NSROptItrStep=1, NSROptItrSmp=1, DSROptRedCut=1e-12,
                   DSROptStaDel=0.01, DSROptStepDt=0.002,
                   NVMCWarmUp=256, NVMCInterval=1, NVMCSample=samples,
                   NExUpdatePath=0, RndSeed=seed, NSplitSize=np_, NStore=1,
                   NSRCG=0, NGrandCanonical=1, NGCInitNelec=initial)
    options['2Sz'] = 0
    (work/'modpara.def').write_text(table('Model_Parameters', 0, [
        'CDataFileHead zvo', 'CParaFileHead zqp', '===================='] + [
        '{} {}'.format(k, v) for k, v in options.items()]))
    (work/'locspn.def').write_text(table('NlocalSpin', 0,
        ['{} 0'.format(i) for i in range(L)]))
    signs = np.ones((L, L), dtype=int)
    if apbc:
        signs[np.tril_indices(L, -1)] = -1
    (work/'orbitalidx.def').write_text(table('NOrbitalIdx', L*L,
        ['{} {} {} {}'.format(i, j, i*L+j, signs[i, j])
         for i in range(L) for j in range(L)] +
        ['{} 1'.format(i) for i in range(L*L)], True))
    params = F/signs
    initial_values = [0.0]*6
    for z in params.flat: initial_values.extend([z.real, z.imag, 0])
    (work/'initial.def').write_text(' '.join(format(x, '.17g') for x in initial_values)+'\n')
    # k=0 combined projection exercises signed APBC translations without
    # removing the vacuum/full support of the minimal regression.
    translations = list(range(L)) if combined else [0]
    rows = ['{} {:.17g} 0'.format(r, 1/len(translations)) for r in translations]
    for r in translations:
        rows += ['{} {} {} {}'.format(r, i, (i+r)%L,
                 -1 if apbc and i+r >= L else 1) for i in range(L)]
    (work/'qptransidx.def').write_text(table('NQPTrans', len(translations), rows))
    hopping = []
    for i in range(L):
        j = (i+1)%L
        if L == 2 and i == 1: continue
        t = -1 if apbc and j == 0 else 1
        for s in range(2):
            hopping += [(i+s*L, j+s*L, t), (j+s*L, i+s*L, t)]
    for i in range(2*L): hopping.append((i, i, 0.25))
    (work/'trans.def').write_text(table('NTransfer', len(hopping),
        ['{} {} {} {} {:.17g} 0'.format(i%L, i//L, j%L, j//L, t)
         for i, j, t in hopping]))
    (work/'coulomb.def').write_text(table('NCoulombIntra', L,
        ['{} 1.5'.format(i) for i in range(L)]))
    keys = [(i, j) for i in range(2*L) for j in range(2*L)]
    (work/'greenone.def').write_text(table('NCisAjs', len(keys),
        ['{} {} {} {}'.format(i%L, i//L, j%L, j//L) for i, j in keys]))
    (work/'anomalousg.def').write_text(table('NAnomalousG', 2,
        ['1 0 0 1 1', '0 1 1 0 0']))
    names = [('ModPara', 'modpara'), ('LocSpin', 'locspn'),
             ('OrbitalAntiParallel', 'orbitalidx'), ('TransSym', 'qptransidx'),
             ('Trans', 'trans'), ('CoulombIntra', 'coulomb'),
             ('OneBodyG', 'greenone'), ('AnomalousG', 'anomalousg')]
    (work/'namelist.def').write_text('\n'.join(k+' '+v+'.def' for k, v in names)+'\n')
    return hopping, translations


def reference(F, gauss, translations, apbc):
    L = len(F)
    wave = {mask: 0j for mask in spin.even_masks(L)}
    for r in translations:
        translated = [[F[(i+r)%L, (j+r)%L] *
            (-1 if apbc and i+r >= L else 1) *
            (-1 if apbc and j+r >= L else 1)
            for j in range(L)] for i in range(L)]
        component = spin.finite_fock(translated, spin.independent_nodes(gauss))
        for mask in component: wave[mask] += component[mask]/len(translations)
    triplet = np.array_equal(F, -F.T)
    if triplet and gauss > 1:
        for mask in wave:
            if mask.bit_count() % 4 == 2: wave[mask] = 0j
    norm = sum(abs(z)**2 for z in wave.values())
    return wave, {m: abs(z)**2/norm for m, z in wave.items()}


def compare(label, actual, expected, tolerance=2e-8):
    if not np.isfinite(actual) or abs(actual-expected) > tolerance*(1+abs(expected)):
        raise AssertionError('{}: {} != {}'.format(label, actual, expected))


def statistical(label, values, expected):
    values = np.asarray(values, dtype=float)
    blocks = values.reshape(-1, 256).mean(axis=1)
    se = blocks.std(ddof=1)/math.sqrt(len(blocks))
    if abs(values.mean()-expected) > 5*se+0.005:
        raise AssertionError('{} mean={} exact={} se={}'.format(
            label, values.mean(), expected, se))
    variance = values.var()
    tau = len(values)*se*se/variance if variance > 0 else 0
    return dict(mean=float(values.mean()), exact=float(expected), se=float(se),
                variance=float(variance), tau_block=float(tau))


def run(binary, root, kind, np_, samples=8192):
    F = fixture(kind)
    L = len(F)
    combined = kind == 'combined'
    apbc = combined
    initials = sorted(set((0, 4, 2*L)))
    output = []
    for gauss in ((2, 8) if kind in ('minimal', 'triplet') else (8,)):
        wave, probability = reference(F, gauss, list(range(L)) if combined else [0], apbc)
        for seed in (17, 29):
            for initial in initials:
                work = root/('{}-g{}-N{}-s{}-np{}'.format(kind, gauss, initial, seed, np_))
                hopping, _ = write_input(work, F, gauss, initial, seed, samples, np_, combined, apbc)
                env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1',
                           MVMC_GC_STATE_DUMP='state.dat')
                command = [str(binary), '-e', 'namelist.def', 'initial.def']
                if np_ > 1: command = mpi_command(np_, command)
                start = time.monotonic()
                proc = subprocess.run(command, cwd=work, env=env, text=True,
                    stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=120)
                elapsed = time.monotonic()-start
                (work/'run.log').write_text(proc.stdout)
                if proc.returncode: raise AssertionError(proc.stdout[-4000:])
                records = state_dump_records(work/'state.dat.rank0')
                masks = [r[3] for r in records if r[0] == 'SAMPLE']
                if len(masks) != samples: raise AssertionError('sample count')
                for rank in range(1, np_):
                    if (work/('state.dat.rank{}'.format(rank))).read_bytes() != (work/'state.dat.rank0').read_bytes():
                        raise AssertionError('QP ranks followed different chains')
                counts = Counter(masks)
                for mask in counts:
                    if probability.get(mask, 0) == 0: raise AssertionError('sampled exact zero')
                numbers = np.array([m.bit_count() for m in masks])
                energies = {}
                for m in probability:
                    if probability[m] == 0: continue
                    energies[m] = 1.5*sum((m>>i&1)*(m>>(i+L)&1) for i in range(L)) - sum(
                        t*local_of(wave, m, one_body(i, j)) for i, j, t in hopping)
                exact_n = sum(m.bit_count()*p for m, p in probability.items())
                exact_e = sum(energies.get(m, 0)*p for m, p in probability.items())
                stats = statistical('N', numbers, exact_n)
                statistical('N squared', numbers**2,
                    sum(m.bit_count()**2*p for m, p in probability.items()))
                for n in range(0, 2*L+1, 2):
                    statistical('P(N={})'.format(n), numbers == n,
                        sum(p for m, p in probability.items() if m.bit_count() == n))
                statistical('energy', [energies[m].real for m in masks], exact_e.real)
                seen_gc = np.loadtxt(work/'zvo_gc.dat').reshape(-1)
                compare('replay N', seen_gc[0], numbers.mean())
                compare('replay N variance', seen_gc[2], numbers.var())
                seen_e = np.loadtxt(work/'output/zvo_out_001.dat').reshape(-1)
                mean_e = sum(energies[m]*count for m, count in counts.items())/samples
                compare('replay energy', seen_e[0]+1j*seen_e[1], mean_e)
                greens = np.loadtxt(work/'output/zvo_cisajs_001.dat')
                for row in greens:
                    i, s, j, t = map(int, row[:4])
                    expect = sum(count*local_of(wave, m, one_body(i+s*L, j+t*L))
                                 for m, count in counts.items())/samples
                    compare('replay Green', row[4]+1j*row[5], expect)
                for row in np.loadtxt(work/'output/zvo_anomalousg_001.dat'):
                    kind_, i, s, j, t = map(int, row[:5])
                    letter = 'c' if kind_ else 'a'
                    ops = ((letter, i+s*L), (letter, j+t*L))
                    expect = sum(count*local_of(wave, m, ops)
                                 for m, count in counts.items())/samples
                    compare('replay anomalous Green', row[5]+1j*row[6], expect)
                stats.update(case=kind, gauss=gauss, initial=initial, seed=seed,
                             np=np_, seconds=elapsed, histogram=dict(Counter(numbers.tolist())))
                output.append(stats)
                print(json.dumps(stats), flush=True)
    (root/'results.json').write_text(json.dumps(output, indent=2)+'\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('binary', type=Path)
    parser.add_argument('--case', choices=['minimal', 'triplet', 'near', 'singlet',
                                         'general', 'combined', 'six'], required=True)
    parser.add_argument('--np', type=int, default=1)
    parser.add_argument('--work', type=Path, required=True)
    args = parser.parse_args()
    args.work.mkdir(parents=True, exist_ok=True)
    root = Path(tempfile.mkdtemp(prefix='run-', dir=args.work.resolve()))
    run(args.binary.resolve(), root, args.case, args.np)
