"""Generate and audit a Z2 seam for rectangular expert-input fixtures.

Numbering is explicit: x-fast means site=x+Lx*y, y-fast means site=y+Ly*x.
Only a unique shortest displacement is supported; antipodal bonds are rejected.
Transfer verification takes the unseamed (real) hopping amplitude explicitly.
"""
import argparse
import math
from pathlib import Path
import re


INTEGER = re.compile(r"[+-]?[0-9]+\Z")


def definition_rows(path, widths):
    lines = Path(path).read_text().splitlines()
    start = next((i for i, line in enumerate(lines)
                  if line.split() and INTEGER.fullmatch(line.split()[0])), None)
    if start is None:
        raise ValueError("no data rows: {}".format(path))
    rows = []
    for number, line in enumerate(lines[start:], start+1):
        fields = line.split()
        if len(fields) not in widths:
            raise ValueError("{}:{}: wrong column count".format(path, number))
        rows.append(fields)
    return lines[:start], rows


class RectangularSeam:
    def __init__(self, lx, ly, ap_axes, numbering):
        if lx <= 0 or ly <= 0 or numbering not in ("x-fast", "y-fast"):
            raise ValueError("invalid lattice dimensions/numbering")
        self.lx, self.ly, self.numbering = lx, ly, numbering
        self.ap = set(ap_axes)
        if not self.ap <= {"x", "y"}:
            raise ValueError("AP axes must be x and/or y")

    def coordinates(self, site):
        if not 0 <= site < self.lx*self.ly:
            raise ValueError("site outside lattice: {}".format(site))
        return ((site % self.lx, site // self.lx) if self.numbering == "x-fast"
                else (site // self.ly, site % self.ly))

    def site(self, x, y):
        x, y = x % self.lx, y % self.ly
        return x+self.lx*y if self.numbering == "x-fast" else y+self.ly*x

    def phase(self, left, right):
        result = 1
        for axis, length, origin, target in zip(
                ("x", "y"), (self.lx, self.ly),
                self.coordinates(left), self.coordinates(right)):
            delta = target-origin
            if length > 1 and 2*abs(delta) == length:
                raise ValueError("ambiguous antipodal bond {}->{} on {}".format(left, right, axis))
            displacement = delta
            if 2*displacement > length:
                displacement -= length
            elif 2*displacement < -length:
                displacement += length
            winding = (origin+displacement-target)//length
            if axis in self.ap and winding % 2:
                result = -result
        return result

    def translation(self, dx, dy):
        mapping, signs = [], []
        for site in range(self.lx*self.ly):
            x, y = self.coordinates(site)
            winding = ((x+dx)//self.lx if "x" in self.ap else 0)
            winding += ((y+dy)//self.ly if "y" in self.ap else 0)
            mapping.append(self.site(x+dx, y+dy))
            signs.append(-1 if winding % 2 else 1)
        return mapping, signs


def range_table(path):
    header, raw = definition_rows(path, (3, 4))
    rows, table = [], {}
    for fields in raw:
        if not all(INTEGER.fullmatch(x) for x in fields):
            raise ValueError("BFRange requires integer columns")
        row = tuple(map(int, fields))
        i, k, shell = row[:3]
        phase = row[3] if len(row) == 4 else 1
        if (i, k) in table or min(i, k, shell) < 0 or phase not in (-1, 1):
            raise ValueError("duplicate/invalid BFRange entry")
        if i == k and phase != 1:
            raise ValueError("self phase must be +1")
        rows.append((i, k, shell, phase))
        table[i, k] = (shell, phase)
    if len({len(fields) for fields in raw}) != 1:
        raise ValueError("mixed BFRange columns")
    for i, k in table:
        if table.get((k, i)) != table[i, k]:
            raise ValueError("asymmetric BFRange")
    return header, rows, table


def validate_transform(table, mapping, signs):
    size = len(mapping)
    if sorted(mapping) != list(range(size)) or len(signs) != size or any(s not in (-1, 1) for s in signs):
        raise ValueError("invalid signed permutation")
    for (i, k), (shell, phase) in table.items():
        if not (0 <= i < size and 0 <= k < size):
            raise ValueError("BFRange site outside transform")
        if table.get((mapping[i], mapping[k])) != (shell, signs[i]*signs[k]*phase):
            raise ValueError("seam covariance violation at {}->{}".format(i, k))


def validate_transfer(path, table, bare_hopping):
    if not math.isfinite(bare_hopping) or bare_hopping == 0:
        raise ValueError("bare hopping must be finite and nonzero")
    _, rows = definition_rows(path, (6,))
    checked = 0
    for row in rows:
        i, si, k, sk = map(int, row[:4])
        value = complex(float(row[4]), float(row[5]))
        if not (math.isfinite(value.real) and math.isfinite(value.imag)):
            raise ValueError("nonfinite Transfer value")
        if si != sk or i == k or (i, k) not in table:
            continue
        expected = bare_hopping*table[i, k][1]
        if abs(value-expected) > 1e-12*max(1, abs(expected)):
            raise ValueError("Transfer/seam mismatch at {}->{}".format(i, k))
        checked += 1
    if not checked:
        raise ValueError("no overlapping nonlocal Transfer bonds")
    return checked


def validate_orbital_subgroup(path, mapping, signs):
    size = len(mapping)
    # The index table is followed by two-column optimization flags.
    _, rows = definition_rows(path, (2, 3, 4))
    table = {}
    for fields in rows:
        if len(fields) == 2:
            continue
        row = list(map(int, fields))
        i, j, idx = row[:3]
        if (i, j) in table or not (0 <= i < size and 0 <= j < size):
            raise ValueError("invalid orbital table")
        table[i, j] = (idx, row[3] if len(row) == 4 else 1)
    if len(table) != size*size:
        raise ValueError("incomplete non-FSZ orbital table")
    for (i, j), (idx, sign) in table.items():
        if table[mapping[i], mapping[j]] != (idx, signs[i]*signs[j]*sign):
            raise ValueError("orbital is not covariant under the declared subgroup at {},{}".format(i, j))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=("generate", "check", "subgroup"))
    parser.add_argument("--lx", type=int, required=True)
    parser.add_argument("--ly", type=int, required=True)
    parser.add_argument("--ap", choices=("none", "x", "y", "xy"), required=True)
    parser.add_argument("--numbering", choices=("x-fast", "y-fast"), required=True)
    parser.add_argument("--range", type=Path, required=True)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--transfer", type=Path)
    parser.add_argument("--bare-hopping", type=float)
    parser.add_argument("--translation", nargs=2, type=int, action="append", default=[])
    parser.add_argument("--orbital", type=Path)
    args = parser.parse_args()
    grid = RectangularSeam(args.lx, args.ly, "" if args.ap == "none" else args.ap, args.numbering)
    header, rows, table = range_table(args.range)
    if args.command == "generate":
        if args.output is None or args.output.resolve() == args.range.resolve():
            parser.error("generate requires a separate --output")
        rows = [(i, k, shell, grid.phase(i, k)) for i, k, shell, _ in rows]
        table = {(i, k): (shell, phase) for i, k, shell, phase in rows}
    else:
        for i, k, _, phase in rows:
            if phase != grid.phase(i, k):
                raise ValueError("phase differs from specified seam at {}->{}".format(i, k))
    for dx, dy in args.translation:
        mapping, signs = grid.translation(dx, dy)
        validate_transform(table, mapping, signs)
        if args.command == "subgroup":
            if args.orbital is None:
                parser.error("subgroup requires --orbital")
            validate_orbital_subgroup(args.orbital, mapping, signs)
    if args.command == "subgroup" and not args.translation:
        parser.error("subgroup requires at least one --translation")
    if args.transfer:
        if args.bare_hopping is None:
            parser.error("Transfer audit requires explicit --bare-hopping")
        print("Transfer bonds checked:", validate_transfer(args.transfer, table, args.bare_hopping))
    if args.command == "generate":
        args.output.write_text("\n".join(header+["{} {} {} {}".format(*row) for row in rows])+"\n")
    print("seam entries:", len(rows), "negative:", sum(p < 0 for _, _, _, p in rows))


if __name__ == "__main__":
    main()
