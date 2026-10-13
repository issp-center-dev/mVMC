"""Independent tensor contraction and coefficient/gauge/projector audits.

The lists below represent directed single-center channels. Their outer products
are contracted with the BF parameter matrix and orbital matrix; no C loop
selection or reverse-derivative formula is reused.
"""
import copy
import random

import numpy as np

from bf_canonical_model import CanonicalModel, _complex_vector, _read_blocks


def active_channels(model, site, spin):
    channels = []
    nr = model.nsite*model.nrange
    for slot, neighbor in enumerate(model.posbf[site]):
        distance = model.rangeidx[site][neighbor]
        for mu in range(4):
            if distance == 0:
                if mu != 0:
                    continue
                index = 0
            else:
                if mu == 0:
                    continue
                index = 3*(distance-1)+mu
            value = model.count[spin*4*nr+mu*nr+site*model.nrange+slot]
            if value:
                channels.append((neighbor, index, value))
    return channels


def directed_coefficients(model, i, j, transform):
    mp, sg = model.transform[transform], model.sign[transform]
    endpoint = sg[i]*sg[j]
    orbital = np.asarray(model.slater)
    ia = np.zeros(model.nslater, complex)
    ib = np.zeros_like(ia)
    pa = np.zeros(model.nprojbf, complex)
    pb = np.zeros_like(pa)
    qa = np.zeros_like(pa)
    qb = np.zeros_like(pa)
    parts = ((ia, pa, qa, 0, 1, False), (ib, pb, qb, 1, 0, True))
    for oc, pc, qc, si, sj, reverse in parts:
        for k, mu, left in active_channels(model, i, si):
            for l, nu, right in active_channels(model, j, sj):
                if mu == 0 and nu == 0:
                    continue
                a, b = (mp[l], mp[k]) if reverse else (mp[k], mp[l])
                index, sign = model.orbital_idx[a][b], model.orbital_sign[a][b]
                scale = -endpoint*left*right*model.seam_phase[mp[i]][mp[k]]*model.seam_phase[mp[j]][mp[l]]*sign
                parameter = model.bfsubidx[mu][nu]
                oc[index] += scale*model.projbf[parameter]
                pc[parameter] += scale*orbital[index]
                qc[parameter] += 1j*scale*orbital[index]
        nr = model.nsite*model.nrange
        active_eta = any(model.count[nr+s*model.nrange+r]
                         for s in (i, j) for r in range(model.nrange))
        a, b = (mp[j], mp[i]) if reverse else (mp[i], mp[j])
        index, sign = model.orbital_idx[a][b], model.orbital_sign[a][b]
        oc[index] += endpoint*sign*(model.projbf[0].real if active_eta else 1)
        if active_eta:
            pc[0] += endpoint*sign*orbital[index]
    return ia, ib, pa, qa, pb, qb


def scaled_error(left, right):
    a, b = np.asarray(left), np.asarray(right)
    return float(np.max(np.abs(a-b)/(1+np.maximum(np.abs(a), np.abs(b)))))


def gauge_error(model, seed=20260929):
    rng = random.Random(seed)
    base = np.asarray(model.build())
    maximum = 0.0
    for _ in range(20):
        gauge = [rng.choice((-1, 1)) for _ in range(model.nsite)]
        transformed = copy.deepcopy(model)
        transformed.nslater = model.nsite**2
        transformed.orbital_idx = np.arange(transformed.nslater).reshape(model.nsite, model.nsite).tolist()
        transformed.orbital_sign = [[1]*model.nsite for _ in range(model.nsite)]
        transformed.slater = [gauge[i]*gauge[j]*model.slater[model.orbital_idx[i][j]]*model.orbital_sign[i][j]
                              for i in range(model.nsite) for j in range(model.nsite)]
        transformed.seam_phase = [[gauge[i]*gauge[k]*model.seam_phase[i][k]
                                   for k in range(model.nsite)] for i in range(model.nsite)]
        transformed.sign = [[sign[i]*gauge[i]*gauge[mapping[i]] for i in range(model.nsite)]
                            for mapping, sign in zip(model.transform, model.sign)]
        g = np.asarray(gauge+gauge)
        expected = base*g[:, None]*g[None, :]
        maximum = max(maximum, scaled_error(transformed.build(), expected))
    return maximum


def direct_configuration_error(model):
    """Apply U to the configuration and rebuild with identity independently.

    This identity holds for arbitrary orbitals; it assumes no symmetry of f.
    Electron labels are preserved, so no extra sorting permutation is present.
    """
    maximum = 0.0
    for transform in range(model.nmp):
        mapping, signs = model.transform[transform], model.sign[transform]
        translated = copy.deepcopy(model)
        translated.ele_idx = [mapping[i] for i in model.ele_idx]
        translated.count = translated._make_count(translated.ele_idx)
        translated.nmp = 1
        translated.nqpfull = model.nsp
        translated.transform = [list(range(model.nsite))]
        translated.sign = [[1]*model.nsite]
        direct = np.asarray(translated.build())
        base = np.asarray(model.build())[transform*model.nsp:(transform+1)*model.nsp]
        indices = np.asarray(mapping+[i+model.nsite for i in mapping])
        g = np.asarray(signs+signs)
        expected = direct[:, indices][:, :, indices]*g[:, None]*g[None, :]
        maximum = max(maximum, scaled_error(base, expected))
    return maximum


def check_coefficients(path):
    maximum = {key: 0.0 for key in ("coefficient", "value", "gauge", "direct_configuration")}
    count = 0
    names = ("orbital_a", "orbital_b", "proj_real_a", "proj_imag_a", "proj_real_b", "proj_imag_b")
    for block in _read_blocks(path):
        model = CanonicalModel(block)
        for mp in range(model.nmp):
            for i in range(model.nsite):
                for j in range(model.nsite):
                    parts = directed_coefficients(model, i, j, mp)
                    for name, values in zip(names, parts):
                        c = _complex_vector(block, "coef_{}_{}_{}_{}".format(name, mp, i, j), len(values))
                        maximum["coefficient"] = max(maximum["coefficient"], scaled_error(values, c))
                        count += len(values)
                    values = (np.dot(parts[0], model.slater), np.dot(parts[1], model.slater))
                    maximum["value"] = max(maximum["value"], scaled_error(values, model.directed(i, j, mp)))
        maximum["gauge"] = max(maximum["gauge"], gauge_error(model))
        maximum["direct_configuration"] = max(maximum["direct_configuration"], direct_configuration_error(model))
    assert count > 0
    assert maximum["coefficient"] <= 1e-12, maximum
    assert maximum["value"] <= 1e-13, maximum
    assert maximum["gauge"] <= 1e-12 and maximum["direct_configuration"] <= 1e-12, maximum
    return dict(maximum, coefficient_count=count)
