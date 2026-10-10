"""Convert the StdFace output of stan.in into a grand-canonical APBC input.

Usage: python3 make_gc_apbc_input.py <vmcdry.out> [output directory]

StdFace (``phase0 = 180``, ``2Sz = -1``) writes the anti-periodic boundary
phase into trans.def and one parameter plus one boundary sign per pair into
orbitalidxgen.def.  This script keeps those files unchanged and only
 * keeps OrbitalGeneral as the single orbital entry of namelist.def,
 * sets NGrandCanonical=1, NGCInitNelec, and NSPGaussLeg=1 in modpara.def
   (2Sz=-1, NMPTrans=-1 and NExUpdatePath=0 are checked, not changed), and
 * writes explicit complex initial parameters to initial.def.
"""
from __future__ import print_function

import math
import os
import shutil
import subprocess
import sys
import tempfile


INIT_NELEC = 4
UNUSED_ORBITAL_FILES = ("orbitalidx.def", "orbitalidxpara.def")


def header_count(path, keyword):
    with open(path) as stream:
        words = stream.readlines()[1].split()
    if words[0] != keyword:
        raise RuntimeError("{}: expected {}".format(path, keyword))
    return int(words[1])


def convert_namelist(path):
    lines = []
    with open(path) as stream:
        for line in stream:
            words = line.split()
            if words and words[0] in ("Orbital", "OrbitalParallel"):
                continue
            if words[:2] == ["#", "OrbitalGeneral"]:
                line = "  OrbitalGeneral  {}\n".format(words[2])
            lines.append(line)
    if not any(line.split()[:1] == ["OrbitalGeneral"] for line in lines):
        raise RuntimeError("StdFace did not offer an OrbitalGeneral file")
    with open(path, "w") as stream:
        stream.writelines(lines)


def convert_modpara(path):
    required = {"2Sz": "-1", "NMPTrans": "-1", "NExUpdatePath": "0"}
    added = (("NSPGaussLeg", "1"), ("NGrandCanonical", "1"),
             ("NGCInitNelec", str(INIT_NELEC)))
    lines = []
    with open(path) as stream:
        for line in stream:
            words = line.split()
            if words and words[0] in required:
                if words[1] != required.pop(words[0]):
                    raise RuntimeError("unexpected {}".format(line.strip()))
            if words and words[0] in dict(added):
                continue
            lines.append(line)
    if required:
        raise RuntimeError("modpara.def lacks {}".format(sorted(required)))
    lines.extend("{:<15} {}\n".format(key, value) for key, value in added)
    with open(path, "w") as stream:
        stream.writelines(lines)


def write_initial(directory):
    nproj = (header_count(os.path.join(directory, "gutzwilleridx.def"),
                          "NGutzwillerIdx") +
             header_count(os.path.join(directory, "jastrowidx.def"),
                          "NJastrowIdx"))
    norbital = header_count(os.path.join(directory, "orbitalidxgen.def"),
                            "NOrbitalIdx")
    values = [0.0] * 6
    for unused in range(nproj):
        values.extend((0.0, 0.0, 0.0))
    for index in range(norbital):
        # Deterministic complex pair amplitudes of order 0.1-0.3.
        values.extend((0.2 * math.cos(0.7 * index + 0.3),
                       0.2 * math.sin(1.3 * index + 0.5), 0.0))
    with open(os.path.join(directory, "initial.def"), "w") as stream:
        stream.write(" ".join("{:.15e}".format(value) for value in values))
        stream.write("\n")


def main():
    if len(sys.argv) not in (2, 3):
        print(__doc__)
        return 2
    vmcdry = os.path.abspath(sys.argv[1])
    output = os.path.abspath(sys.argv[2] if len(sys.argv) == 3 else "expert")
    here = os.path.dirname(os.path.abspath(__file__))
    work = tempfile.mkdtemp(prefix="gc_apbc_stdface_")
    try:
        shutil.copy(os.path.join(here, "stan.in"), work)
        process = subprocess.run([vmcdry, "stan.in"], cwd=work,
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.STDOUT,
                                 universal_newlines=True)
        if process.returncode != 0:
            print(process.stdout)
            return 1
        if not os.path.isdir(output):
            os.makedirs(output)
        for name in sorted(os.listdir(work)):
            if name.endswith(".def") and name not in UNUSED_ORBITAL_FILES:
                shutil.copy(os.path.join(work, name), output)
        convert_namelist(os.path.join(output, "namelist.def"))
        convert_modpara(os.path.join(output, "modpara.def"))
        write_initial(output)
    finally:
        shutil.rmtree(work)
    print("grand-canonical APBC input written to {}".format(output))
    return 0


if __name__ == "__main__":
    sys.exit(main())
