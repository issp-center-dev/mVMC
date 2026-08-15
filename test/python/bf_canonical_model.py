from __future__ import print_function

import math


VALUE_TOL = 2.0e-12
IP_TOL = 2.0e-10
DERIVATIVE_TOL = 2.0e-6
FD_STEP = 1.0e-6


def _read_blocks(path):
    blocks = []
    current = {}
    with open(path) as source:
        for line in source:
            fields = line.split()
            if not fields:
                if current:
                    blocks.append(current)
                    current = {}
                continue
            current[fields[0]] = fields[1:]
    if current:
        blocks.append(current)
    return blocks


def _integer(block, key):
    values = block.get(key, [])
    if len(values) != 1:
        raise ValueError("missing integer field {}".format(key))
    return int(values[0])


def _integers(block, key, count):
    values = block.get(key, [])
    if len(values) != count:
        raise ValueError("invalid integer vector {}".format(key))
    return [int(value) for value in values]


def _complex_value(block, key):
    values = block.get(key, [])
    if len(values) != 2:
        raise ValueError("invalid complex field {}".format(key))
    value = complex(float(values[0]), float(values[1]))
    if not math.isfinite(value.real) or not math.isfinite(value.imag):
        raise ValueError("non-finite complex field {}".format(key))
    return value


def _complex_vector(block, key, count):
    values = block.get(key, [])
    if len(values) != 2 * count:
        raise ValueError("invalid complex vector {}".format(key))
    result = []
    for idx in range(count):
        value = complex(float(values[2 * idx]), float(values[2 * idx + 1]))
        if not math.isfinite(value.real) or not math.isfinite(value.imag):
            raise ValueError("non-finite complex vector {}".format(key))
        result.append(value)
    return result


def _sub_index(mu, distance):
    packed = 4 * distance + mu
    result = packed - 3 - distance
    if packed % 4 == 0:
        result = -1
    if packed == 0:
        result = 0
    return result


class CanonicalModel(object):
    def __init__(self, block):
        self.nsite = _integer(block, "nsite")
        self.nsite2 = _integer(block, "nsite2")
        self.nsize = _integer(block, "nsize")
        self.ne = _integer(block, "ne")
        self.nrange = _integer(block, "nrange")
        self.nrangeidx = _integer(block, "nrangeidx")
        self.nqpfull = _integer(block, "nqpfull")
        self.nsp = _integer(block, "nspgaussleg")
        self.nmp = _integer(block, "nmptrans")
        self.nprojbf = _integer(block, "nprojbf")
        self.nslater = _integer(block, "nslater")
        if self.nsite2 != 2 * self.nsite or self.nqpfull != self.nmp * self.nsp:
            raise ValueError("inconsistent canonical-model dimensions")

        self.ele_idx = _integers(block, "ele_idx", self.nsize)
        self.count = _integers(
            block, "ele_projbf_cnt", 16 * self.nsite * self.nrange)
        self.projbf = [_complex_value(block, "projbf_param_{}".format(idx))
                       for idx in range(self.nprojbf)]
        self.slater = [_complex_value(block, "slater_param_{}".format(idx))
                       for idx in range(self.nslater)]
        self.transform = [
            _integers(block, "qp_transform_{}".format(idx), self.nsite)
            for idx in range(self.nmp)]
        self.sign = [
            _integers(block, "qp_sign_{}".format(idx), self.nsite)
            for idx in range(self.nmp)]
        self.orbital_idx = [
            _integers(block, "orbital_idx_{}".format(site), self.nsite)
            for site in range(self.nsite)]
        self.orbital_sign = [
            _integers(block, "orbital_sgn_{}".format(site), self.nsite)
            for site in range(self.nsite)]
        self.posbf = [
            _integers(block, "posbf_{}".format(site), self.nrange)
            for site in range(self.nsite)]
        self.rangeidx = [
            _integers(block, "rangeidx_{}".format(site), self.nsite)
            for site in range(self.nsite)]
        self.bfsubidx = [
            _integers(block, "bfsubidx_{}".format(idx), self.nrangeidx)
            for idx in range(self.nrangeidx)]
        self.spin = []
        for idx in range(self.nsp):
            values = block.get("spgl_{}".format(idx), [])
            if len(values) != 6:
                raise ValueError("invalid spin-projection row {}".format(idx))
            self.spin.append(tuple(
                complex(float(values[2 * column]), float(values[2 * column + 1]))
                for column in range(3)))
        self.weight = [_complex_value(block, "qp_weight_{}".format(idx))
                       for idx in range(self.nqpfull)]
        self.c_value = [[
            _complex_vector(block, "slater_elm_{}_{}".format(qp, row),
                            self.nsite2)
            for row in range(self.nsite2)]
            for qp in range(self.nqpfull)]
        self.c_proj_derivative = []
        for idx in range(self.nprojbf):
            values = block.get("projbf_derivative_{}".format(idx), [])
            if len(values) != 4:
                raise ValueError("invalid ProjBF derivative row")
            self.c_proj_derivative.append((
                complex(float(values[0]), float(values[1])),
                complex(float(values[2]), float(values[3]))))
        self.c_orbital_derivative = []
        for idx in range(self.nslater):
            values = block.get("orbital_derivative_{}".format(idx), [])
            if len(values) != 4:
                raise ValueError("invalid orbital derivative row")
            self.c_orbital_derivative.append((
                complex(float(values[0]), float(values[1])),
                complex(float(values[2]), float(values[3]))))
        self.c_ip = _complex_value(block, "ip_base")

    def orbital(self, raw_i, raw_j, transform_index,
                virtual_neighbor_signs=False):
        mapping = self.transform[transform_index]
        transformed_i = mapping[raw_i]
        transformed_j = mapping[raw_j]
        orbital_index = self.orbital_idx[transformed_i][transformed_j]
        sign = self.orbital_sign[transformed_i][transformed_j]
        if virtual_neighbor_signs:
            sign *= (self.sign[transform_index][raw_i] *
                     self.sign[transform_index][raw_j])
        return self.slater[orbital_index] * sign

    def directed(self, ri, rj, transform_index,
                 virtual_neighbor_signs=False):
        nsite_range = self.nsite * self.nrange
        slater_ij = 0.0j
        slater_ji = 0.0j
        ij_count = 0
        ji_count = 0
        for xn in range(4):
            for xm in range(4):
                if xn == 0 and xm == 0:
                    continue
                for xk in range(self.nrange):
                    rk = self.posbf[ri][xk]
                    idx_ik = ri * self.nrange + xk
                    idx_jk = rj * self.nrange + xk
                    nidx = _sub_index(xn, self.rangeidx[ri][rk])
                    if nidx < 0:
                        continue
                    ij_count += self.count[nsite_range + idx_ik]
                    ij_count += self.count[nsite_range + idx_jk]
                    ji_count += self.count[nsite_range + idx_ik]
                    ji_count += self.count[nsite_range + idx_jk]
                    count_ij_i = self.count[xn * nsite_range + idx_ik]
                    count_ji_i = self.count[
                        4 * nsite_range + xn * nsite_range + idx_ik]
                    if count_ij_i == 0 and count_ji_i == 0:
                        continue
                    for xl in range(self.nrange):
                        rl = self.posbf[rj][xl]
                        idx_jl = rj * self.nrange + xl
                        count_ij_j = self.count[
                            4 * nsite_range + xm * nsite_range + idx_jl]
                        count_ji_j = self.count[xm * nsite_range + idx_jl]
                        midx = _sub_index(xm, self.rangeidx[rj][rl])
                        if midx < 0:
                            continue
                        if ((count_ij_i == 0 or count_ij_j == 0) and
                                (count_ji_i == 0 or count_ji_j == 0)):
                            continue
                        parameter = self.bfsubidx[nidx][midx]
                        if count_ij_i != 0 and count_ij_j != 0:
                            slater_ij -= (self.projbf[parameter] *
                                           count_ij_i * count_ij_j *
                                           self.orbital(
                                               rk, rl, transform_index,
                                               virtual_neighbor_signs))
                        if count_ji_i != 0 and count_ji_j != 0:
                            slater_ji -= (self.projbf[parameter] *
                                           count_ji_i * count_ji_j *
                                           self.orbital(
                                               rl, rk, transform_index,
                                               virtual_neighbor_signs))
        eta_ij = 1.0 if ij_count == 0 else self.projbf[0].real
        eta_ji = 1.0 if ji_count == 0 else self.projbf[0].real
        slater_ij += eta_ij * self.orbital(
            ri, rj, transform_index, virtual_neighbor_signs)
        slater_ji += eta_ji * self.orbital(
            rj, ri, transform_index, virtual_neighbor_signs)
        if virtual_neighbor_signs:
            return slater_ij, slater_ji
        endpoint_sign = (self.sign[transform_index][ri] *
                         self.sign[transform_index][rj])
        return endpoint_sign * slater_ij, endpoint_sign * slater_ji

    def build(self, virtual_neighbor_signs=False):
        result = []
        for qp in range(self.nqpfull):
            transform_index = qp // self.nsp
            spin_index = qp % self.nsp
            cs, cc, ss = self.spin[spin_index]
            matrix = [[0.0j for _ in range(self.nsite2)]
                      for _ in range(self.nsite2)]
            for ri in range(self.nsite):
                for rj in range(self.nsite):
                    directed_ij, directed_ji = self.directed(
                        ri, rj, transform_index,
                        virtual_neighbor_signs)
                    matrix[ri][rj] = -(directed_ij - directed_ji) * cs
                    matrix[ri][rj + self.nsite] = directed_ij * cc + directed_ji * ss
                    matrix[ri + self.nsite][rj] = -directed_ij * ss - directed_ji * cc
                    matrix[ri + self.nsite][rj + self.nsite] = (
                        directed_ij - directed_ji) * cs
            result.append(matrix)
        return result

    def projected_ip(self, virtual_neighbor_signs=False):
        matrices = self.build(virtual_neighbor_signs)
        total = 0.0j
        for qp, slater_matrix in enumerate(matrices):
            occupied = [[0.0j for _ in range(self.nsize)]
                        for _ in range(self.nsize)]
            for left in range(self.nsize):
                left_orbital = self.ele_idx[left] + (left // self.ne) * self.nsite
                for right in range(self.nsize):
                    right_orbital = (self.ele_idx[right] +
                                     (right // self.ne) * self.nsite)
                    occupied[left][right] = -slater_matrix[left_orbital][right_orbital]
            total += self.weight[qp] * pfaffian(occupied)
        return total

    def _make_count(self, ele_idx):
        n0 = [0 for _ in range(self.nsite)]
        n1 = [0 for _ in range(self.nsite)]
        for particle, site in enumerate(ele_idx):
            (n0 if particle < self.ne else n1)[site] += 1
        nsite_range = self.nsite * self.nrange
        result = [0 for _ in range(16 * nsite_range)]
        for ri in range(self.nsite):
            xid = n0[ri] * n1[ri]
            xih = (1 - n0[ri]) * (1 - n1[ri])
            xidh = n0[ri] * (1 - n1[ri])
            xihd = n1[ri] * (1 - n0[ri])
            for slot, rk in enumerate(self.posbf[ri]):
                idx = ri * self.nrange + slot
                xkd = n0[rk] * n1[rk]
                xkh = (1 - n0[rk]) * (1 - n1[rk])
                xkdh = n0[rk] * (1 - n1[rk])
                xkhd = n1[rk] * (1 - n0[rk])
                diagonal = 1 if ri == rk else 0
                for group in range(4):
                    result[group * 4 * nsite_range + idx] = diagonal
                result[nsite_range + idx] = xid * xkh
                result[2 * nsite_range + idx] = xidh * xkhd
                result[3 * nsite_range + idx] = (
                    xid * xkhd + xidh * xkh)
                group1 = 4 * nsite_range
                result[group1 + nsite_range + idx] = xid * xkh
                result[group1 + 2 * nsite_range + idx] = xihd * xkdh
                result[group1 + 3 * nsite_range + idx] = (
                    xid * xkdh + xihd * xkh)
                group2 = 8 * nsite_range
                result[group2 + nsite_range + idx] = xkd * xih
                result[group2 + 2 * nsite_range + idx] = xkdh * xihd
                result[group2 + 3 * nsite_range + idx] = (
                    xkd * xihd + xkdh * xih)
                group3 = 12 * nsite_range
                result[group3 + nsite_range + idx] = xkd * xih
                result[group3 + 2 * nsite_range + idx] = xkhd * xidh
                result[group3 + 3 * nsite_range + idx] = (
                    xkd * xidh + xkhd * xih)
        return result

    def _translated_state(self, transform_index):
        mapping = self.transform[transform_index]
        signs = self.sign[transform_index]
        translated = []
        state_sign = 1
        for site in self.ele_idx:
            translated.append(mapping[site])
            state_sign *= signs[site]
        return translated, self._make_count(translated), state_sign

    def _composition(self, left, right):
        mapping = []
        signs = []
        for site in range(self.nsite):
            intermediate = self.transform[right][site]
            mapping.append(self.transform[left][intermediate])
            signs.append(self.sign[right][site] *
                         self.sign[left][intermediate])
        for candidate in range(self.nmp):
            if self.transform[candidate] != mapping:
                continue
            ratios = [signs[site] * self.sign[candidate][site]
                      for site in range(self.nsite)]
            if all(value == ratios[0] for value in ratios):
                return candidate, ratios[0]
        raise ValueError("QP transform rows are not closed under composition")

    def _projection_character(self, translated_by):
        ratios = []
        for left in range(self.nmp):
            composed, global_sign = self._composition(left, translated_by)
            many_body_sign = global_sign ** self.nsize
            for spin_index in range(self.nsp):
                left_weight = self.weight[left * self.nsp + spin_index]
                composed_weight = self.weight[
                    composed * self.nsp + spin_index]
                if abs(left_weight) <= 1.0e-14 and abs(composed_weight) <= 1.0e-14:
                    continue
                if abs(composed_weight) <= 1.0e-14:
                    raise ValueError("projection weights do not form a character")
                ratios.append(many_body_sign * left_weight / composed_weight)
        if not ratios:
            raise ValueError("projection character has no nonzero weight")
        if max(abs(value - ratios[0]) for value in ratios) > 2.0e-10:
            raise ValueError("projection weights are inconsistent with QP group")
        return ratios[0]

    def _with_state(self, ele_idx, count, evaluator):
        saved_idx = self.ele_idx
        saved_count = self.count
        self.ele_idx = ele_idx
        self.count = count
        try:
            return evaluator()
        finally:
            self.ele_idx = saved_idx
            self.count = saved_count

    def _raw_finite_difference(self, parameters, index, imaginary,
                               virtual_neighbor_signs=False):
        original = parameters[index]
        displacement = (1.0j if imaginary else 1.0) * FD_STEP
        parameters[index] = original + displacement
        plus = self.projected_ip(virtual_neighbor_signs)
        parameters[index] = original - displacement
        minus = self.projected_ip(virtual_neighbor_signs)
        parameters[index] = original
        return (plus - minus) / (2.0 * FD_STEP)

    def _most_active_derivative(self, parameters, imaginary, start=0):
        values = []
        for index in range(start, len(parameters)):
            values.append((abs(self._raw_finite_difference(
                parameters, index, imaginary)), index))
        magnitude, index = max(values)
        if magnitude <= 1.0e-12:
            raise ValueError("physical covariance derivative is vacuous")
        return index, self._raw_finite_difference(
            parameters, index, imaginary)

    def physical_covariance_errors(self):
        candidates = [index for index in range(1, self.nmp)
                      if self.transform[index] != list(range(self.nsite))]
        if not candidates:
            raise ValueError("physical covariance needs a nonidentity transform")
        transform_index = candidates[0]
        translated_idx, translated_count, state_sign = self._translated_state(
            transform_index)
        character = self._projection_character(transform_index)
        phase = state_sign * character
        base_value = self.projected_ip()
        translated_value = self._with_state(
            translated_idx, translated_count, self.projected_ip)
        scale = max(abs(base_value), abs(translated_value), 1.0e-12)
        errors = {"value": abs(translated_value - phase * base_value) / scale}
        derivative_cases = (
            ("projbf-real", self.projbf, False, 0),
            ("projbf-imag", self.projbf, True, 1),
            ("orbital-real", self.slater, False, 0),
            ("orbital-imag", self.slater, True, 0),
        )
        for name, parameters, imaginary, start in derivative_cases:
            index, base = self._most_active_derivative(
                parameters, imaginary, start)
            translated = self._with_state(
                translated_idx, translated_count,
                lambda p=parameters, i=index, im=imaginary:
                self._raw_finite_difference(p, i, im))
            scale = max(abs(base), abs(translated), 1.0e-12)
            errors[name] = abs(translated - phase * base) / scale

        virtual_base = self.projected_ip(True)
        virtual_translated = self._with_state(
            translated_idx, translated_count,
            lambda: self.projected_ip(True))
        virtual_scale = max(
            abs(virtual_base), abs(virtual_translated), 1.0e-12)
        virtual_error = abs(
            virtual_translated - phase * virtual_base) / virtual_scale
        return errors, virtual_error


def pfaffian(matrix):
    size = len(matrix)
    if size % 2 != 0:
        raise ValueError("Pfaffian matrix size must be even")
    if size == 0:
        return 1.0 + 0.0j
    value = 0.0 + 0.0j
    for partner in range(1, size):
        reduced = []
        keep = [idx for idx in range(size) if idx not in (0, partner)]
        for row in keep:
            reduced.append([matrix[row][column] for column in keep])
        sign = 1.0 if partner % 2 == 1 else -1.0
        value += sign * matrix[0][partner] * pfaffian(reduced)
    return value


def _max_matrix_difference(left, right):
    maximum = 0.0
    for left_qp, right_qp in zip(left, right):
        for left_row, right_row in zip(left_qp, right_qp):
            for left_value, right_value in zip(left_row, right_row):
                maximum = max(maximum, abs(left_value - right_value))
    return maximum


def _finite_difference(model, parameters, index, imaginary):
    original = parameters[index]
    displacement = (1.0j if imaginary else 1.0) * FD_STEP
    parameters[index] = original + displacement
    plus = model.projected_ip()
    parameters[index] = original - displacement
    minus = model.projected_ip()
    parameters[index] = original
    return ((plus - minus) / (2.0 * FD_STEP)) / model.c_ip


def check_dump(path):
    blocks = _read_blocks(path)
    if not blocks:
        raise ValueError("empty canonical BackFlow dump")
    maxima = {"value": 0.0, "ip": 0.0, "projbf": 0.0, "orbital": 0.0,
              "covariance": 0.0, "virtual_mutation": 0.0,
              "physical_checks": 0}
    for block in blocks:
        model = CanonicalModel(block)
        rebuilt = model.build()
        maxima["value"] = max(
            maxima["value"], _max_matrix_difference(rebuilt, model.c_value))
        model_ip = model.projected_ip()
        maxima["ip"] = max(maxima["ip"], abs(model_ip - model.c_ip))
        for idx in range(model.nprojbf):
            real_fd = _finite_difference(model, model.projbf, idx, False)
            maxima["projbf"] = max(
                maxima["projbf"], abs(real_fd - model.c_proj_derivative[idx][0]))
            if idx != 0:
                imag_fd = _finite_difference(model, model.projbf, idx, True)
                maxima["projbf"] = max(
                    maxima["projbf"], abs(imag_fd - model.c_proj_derivative[idx][1]))
        for idx in range(model.nslater):
            real_fd = _finite_difference(model, model.slater, idx, False)
            imag_fd = _finite_difference(model, model.slater, idx, True)
            maxima["orbital"] = max(
                maxima["orbital"],
                abs(real_fd - model.c_orbital_derivative[idx][0]),
                abs(imag_fd - model.c_orbital_derivative[idx][1]))
        if model.nmp == model.nsite:
            if model._make_count(model.ele_idx) != model.count:
                raise ValueError("Python/C BackFlow count reconstruction mismatch")
            covariance, virtual_mutation = model.physical_covariance_errors()
            maxima["covariance"] = max(
                maxima["covariance"], max(covariance.values()))
            maxima["virtual_mutation"] = max(
                maxima["virtual_mutation"], virtual_mutation)
            maxima["physical_checks"] += 1
    if maxima["value"] > VALUE_TOL:
        raise ValueError("Python/C canonical value mismatch: {:.3e}".format(
            maxima["value"]))
    if maxima["ip"] > IP_TOL:
        raise ValueError("Python/C projected amplitude mismatch: {:.3e}".format(
            maxima["ip"]))
    if maxima["projbf"] > DERIVATIVE_TOL:
        raise ValueError("Python/C ProjBF derivative mismatch: {:.3e}".format(
            maxima["projbf"]))
    if maxima["orbital"] > DERIVATIVE_TOL:
        raise ValueError("Python/C orbital derivative mismatch: {:.3e}".format(
            maxima["orbital"]))
    if maxima["physical_checks"] and maxima["covariance"] > DERIVATIVE_TOL:
        raise ValueError("physical BackFlow covariance mismatch: {:.3e}".format(
            maxima["covariance"]))
    if (maxima["physical_checks"] and
            maxima["virtual_mutation"] <= 1.0e-6):
        raise ValueError(
            "virtual-neighbor sign mutation was not detected: {:.3e}".format(
                maxima["virtual_mutation"]))
    return maxima
