#!/usr/bin/env python3
"""Write a grand-canonical S=0 spin-projection input (expert mode).

Four-site Hubbard ring with chemical potential, nearest-neighbour Coulomb
interaction and an up-down pair field, anti-parallel pairing
(OrbitalAntiParallel), NGrandCanonical=1, 2Sz=0, NSPStot=0 and one
translation.  NSPGaussLeg=G>=2 applies the finite-G S=0 projection;
NSPGaussLeg=1 keeps the unprojected state.  The two-body Green-function
rows give S^2 = <Sz^2> + (1/2) sum_ij <S+_i S-_j + S-_i S+_j>.

Example:
    python3 generate.py --gauss 8 --calmode 1 --samples 256 --seed 4321 \
        --boundary pbc --outdir run_g8
"""
import argparse
import os

NSITE = 4
NORBITAL = 2 * NSITE
GUTZWILLER = -0.31
JASTROW = (0.23, -0.12)  # nearest and next-nearest neighbour classes
HOPPING = 1.11 - 0.13j
MU = 0.35
COULOMB_INTRA = 0.35
COULOMB_INTER = 0.27
DELTA = 0.33 - 0.21j     # pair field on c+_{3 up} c+_{0 dn} and its h.c.
BAR = "=============================================\n"


def fmt(value):
    return "{:.18e}".format(value)


def fields(orbital):
    """Fused orbital index site + spin*NSITE -> (site, spin)."""
    return orbital % NSITE, orbital // NSITE


def def_text(keyword, count, rows, complex_type=None):
    second = ("ComplexType {}\n".format(complex_type)
              if complex_type is not None else BAR)
    return BAR + "{} {}\n".format(keyword, count) + second + BAR + BAR + \
        "".join(rows)


def pair_amplitudes():
    """Non-singlet anti-parallel pair amplitudes F_ij (S=1 and S=2 parts)."""
    F = [[complex((7 + 3 * i - 2 * j) / 16.0, (1 + i + 2 * j) / 32.0)
          for j in range(NSITE)] for i in range(NSITE)]
    F[0][1] += 0.7
    F[2][0] -= 0.3j
    return F


def factors(pairs):
    out = []
    for out_orbital, in_orbital in pairs:
        out.extend(fields(out_orbital))
        out.extend(fields(in_orbital))
    return " ".join(str(value) for value in out)


def two_body_rows():
    """S+_i S-_j, S-_i S+_j, and three further products (one Sz-changing)."""
    pairs = [((i, i + NSITE), (j + NSITE, j)) for i in range(NSITE)
             for j in range(NSITE)]
    pairs += [((i + NSITE, i), (j, j + NSITE)) for i in range(NSITE)
              for j in range(NSITE)]
    pairs += [((0, 0), (5, 5)), ((0, 1), (1, 0)), ((0, 4), (1, 5)),
              ((0, 4), (5, 5))]
    return ["{}\n".format(factors(pair)) for pair in pairs]


def modpara(args):
    data_qty = 1
    return """--------------------
Model_Parameters   0
--------------------
VMC_Cal_Parameters
--------------------
CDataFileHead  zvo
CParaFileHead  zqp
--------------------
NVMCCalMode    {mode}
NLanczosMode   0
--------------------
NDataIdxStart  1
NDataQtySmp    {data_qty}
--------------------
NGrandCanonical 1
NGCInitNelec   4
Nsite          {nsite}
Ne             2
Ncond          4
2Sz            0
NSPGaussLeg    {gauss}
NSPStot        0
NMPTrans       {nmptrans}
NSROptItrStep  1
NSROptItrSmp   1
DSROptRedCut   0.000000000001
DSROptStaDel   0.01
DSROptStepDt   {step}
NVMCWarmUp     {warmup}
NVMCInterval   4
NVMCSample     {samples}
NExUpdatePath  0
RndSeed        {seed}
NSplitSize     1
NStore         1
NSRCG          0
""".format(mode=args.calmode, data_qty=data_qty, nsite=NSITE,
           gauss=args.gauss, nmptrans=-1 if args.boundary == "apbc" else 1,
           step=fmt(0.002), warmup=args.warmup, samples=args.samples,
           seed=args.seed)


def write(outdir, name, text):
    with open(os.path.join(outdir, name), "w") as stream:
        stream.write(text)


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--gauss", type=int, default=8,
                        help="NSPGaussLeg (1 = no spin projection)")
    parser.add_argument("--calmode", type=int, choices=(0, 1), default=1)
    parser.add_argument("--samples", type=int, default=256)
    parser.add_argument("--seed", type=int, default=4321)
    parser.add_argument("--warmup", type=int, default=1000)
    parser.add_argument("--boundary", choices=("pbc", "apbc"), default="pbc")
    parser.add_argument("--outdir", default=".")
    args = parser.parse_args()
    if args.gauss < 1:
        parser.error("--gauss must be positive")
    os.makedirs(args.outdir, exist_ok=True)
    optimize = 1 if args.calmode == 0 else 0
    ap = args.boundary == "apbc"

    write(args.outdir, "modpara.def", modpara(args))
    gutz = ["{} 0\n".format(site) for site in range(NSITE)]
    gutz.append("0 {}\n".format(optimize))
    write(args.outdir, "gutzwilleridx.def",
          def_text("NGutzwillerIdx", 1, gutz, complex_type=0))
    jastrow = []
    for first in range(NSITE):
        for second in range(NSITE):
            if first != second:
                distance = abs(second - first) % NSITE
                distance = min(distance, NSITE - distance)
                jastrow.append("{} {} {}\n".format(first, second,
                                                   distance - 1))
    jastrow.extend("{} {}\n".format(index, optimize)
                   for index in range(len(JASTROW)))
    write(args.outdir, "jastrowidx.def",
          def_text("NJastrowIdx", len(JASTROW), jastrow, complex_type=0))
    write(args.outdir, "locspn.def",
          def_text("NlocalSpin", 0,
                   ["{} 0\n".format(site) for site in range(NSITE)]))
    # Identity translation only (abs(NMPTrans)=1), weight 1.
    trans_sym = ["0 {} {}\n".format(fmt(1.0), fmt(0.0))]
    trans_sym.extend("0 {} {} 1\n".format(site, site)
                     for site in range(NSITE))
    write(args.outdir, "qptransidx.def", def_text("NQPTrans", 1, trans_sym))
    # One parameter per (i,j); the fourth column is the boundary sign.
    orbital = ["{} {} {}{}\n".format(i, j, i * NSITE + j, " 1" if ap else "")
               for i in range(NSITE) for j in range(NSITE)]
    orbital.extend("{} {}\n".format(index, optimize)
                   for index in range(NSITE * NSITE))
    write(args.outdir, "orbitalidx.def",
          def_text("NOrbitalIdx", NSITE * NSITE, orbital, complex_type=1))

    # H = -sum t_ij c+_is c_js - mu N + U sum n_up n_dn + V sum n_i n_j
    #     + Delta c+_{3 up} c+_{0 dn} + h.c.; the ring bond 3-0 carries the
    # boundary phase -1 for APBC.
    rows = ["{} {} {} {} {} {}\n".format(o % NSITE, o // NSITE, o % NSITE,
                                         o // NSITE, fmt(MU), fmt(0.0))
            for o in range(NORBITAL)]
    hopping = {}
    for site in range(NSITE):
        target = (site + 1) % NSITE
        phase = -1 if (ap and target == 0) else 1
        hopping[(site, target)] = HOPPING * phase
        hopping[(target, site)] = HOPPING.conjugate() * phase
    for (first, second), value in sorted(hopping.items()):
        for spin in (0, 1):
            rows.append("{} {} {} {} {} {}\n".format(
                first, spin, second, spin, fmt(value.real), fmt(value.imag)))
    write(args.outdir, "trans.def", def_text("NTransfer", len(rows), rows))
    write(args.outdir, "coulombintra.def",
          def_text("NCoulombIntra", NSITE,
                   ["{} {}\n".format(site, fmt(COULOMB_INTRA))
                    for site in range(NSITE)]))
    bonds = sorted((min(site, (site + 1) % NSITE), max(site, (site + 1) %
                                                       NSITE))
                   for site in range(NSITE))
    write(args.outdir, "coulombinter.def",
          def_text("NCoulombInter", len(bonds),
                   ["{} {} {}\n".format(first, second, fmt(COULOMB_INTER))
                    for first, second in bonds]))

    one = ["{}\n".format(factors(((first, second),)))
           for first in range(NORBITAL) for second in range(NORBITAL)]
    write(args.outdir, "greenone.def", def_text("NCisAjs", len(one), one))
    two = two_body_rows()
    write(args.outdir, "greentwo.def",
          def_text("NCisAjsCktAltDC", len(two), two))
    nbody = ["3 {}\n".format(factors(term)) for term in (
        ((0, 4), (5, 1), (2, 2)), ((0, 4), (1, 1), (2, 2)))]
    write(args.outdir, "nbodyg.def", def_text("NNBodyG", len(nbody), nbody))
    anomalous = ["1 3 0 0 1 {} {}\n".format(fmt(DELTA.real), fmt(DELTA.imag)),
                 "0 0 1 3 0 {} {}\n".format(fmt(DELTA.real),
                                            fmt(-DELTA.imag))]
    write(args.outdir, "anomalousterm.def",
          def_text("NAnomalousTerm", len(anomalous), anomalous))
    names = [
        "         ModPara  modpara.def",
        "         LocSpin  locspn.def",
        "      Gutzwiller  gutzwilleridx.def",
        "         Jastrow  jastrowidx.def",
        "OrbitalAntiParallel  orbitalidx.def",
        "        TransSym  qptransidx.def",
        "           Trans  trans.def",
        "    CoulombIntra  coulombintra.def",
        "    CoulombInter  coulombinter.def",
        "        OneBodyG  greenone.def",
        "        TwoBodyG  greentwo.def",
        "          NBodyG  nbodyg.def",
        "   AnomalousTerm  anomalousterm.def",
    ]
    if args.calmode == 1:
        rows = ["{} {} {} {} {}\n".format(kind, *(fields(a) + fields(b)))
                for kind, a, b in ((1, 3, 4), (0, 4, 3), (1, 4, 3), (0, 3, 4),
                                   (1, 0, 1), (0, 1, 0), (1, 5, 6),
                                   (0, 6, 5))]
        write(args.outdir, "anomalousg.def",
              def_text("NAnomalousG", len(rows), rows))
        names.append("      AnomalousG  anomalousg.def")
    write(args.outdir, "namelist.def", "\n".join(names) + "\n")
    values = [0.0] * 6
    for value in (GUTZWILLER,) + JASTROW:
        values.extend((value, 0.0, 0.0))
    for row in pair_amplitudes():
        for value in row:
            values.extend((value.real, value.imag, 0.0))
    write(args.outdir, "initial.def",
          " ".join(fmt(value) for value in values) + "\n")


if __name__ == "__main__":
    main()
