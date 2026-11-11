"""Compare every binary parameter record with the ordinary text output."""

import math
import os
from pathlib import Path
import struct
import subprocess
import sys
import tempfile


def write_def(directory, name, key, count, rows, complex_type=None):
    header = ["====================", "{} {}".format(key, count),
              "====================", "====================",
              "===================="]
    if complex_type is not None:
        header[2] = "ComplexType {}".format(complex_type)
    (directory / name).write_text("\n".join(header + rows) + "\n")


def prepare_case(directory, complex_type, gutz_count, mode):
    directory.mkdir()
    (directory / "namelist.def").write_text(
        "ModPara modpara.def\nLocSpin locspn.def\nTrans trans.def\n"
        "CoulombIntra coulombintra.def\nGutzwiller gutzwilleridx.def\n"
        "Orbital orbitalidx.def\nTransSym qptransidx.def\n")
    (directory / "modpara.def").write_text(
        "--------------------\nModel_Parameters 0\n--------------------\n"
        "VMC_Cal_Parameters\n--------------------\n"
        "CDataFileHead zvo\nCParaFileHead zqp\n--------------------\n"
        "NVMCCalMode {}\nNLanczosMode 0\nNDataIdxStart 7\n"
        "NDataQtySmp 2\nNsite 4\nNcond 4\n2Sz 0\n"
        "NSPGaussLeg 1\nNSPStot 0\nNMPTrans 1\n"
        "NSROptItrStep 3\nNSROptItrSmp 1\nDSROptRedCut 1e-8\n"
        "DSROptStaDel 0.01\nDSROptStepDt 0.003\nNVMCWarmUp 4\n"
        "NVMCInterval 1\nNVMCSample 32\nNExUpdatePath 0\n"
        "RndSeed 2468\nNSplitSize 1\nNStore 1\nNSRCG 0\n".format(mode))
    write_def(directory, "locspn.def", "NlocalSpin", 0,
              ["{} 0".format(i) for i in range(4)])
    transfers = []
    for i in range(4):
        j = (i + 1) % 4
        for spin in range(2):
            transfers.extend(["{} {} {} {} 1.0 0.0".format(i, spin, j, spin),
                              "{} {} {} {} 1.0 0.0".format(j, spin, i, spin)])
    write_def(directory, "trans.def", "NTransfer", len(transfers), transfers)
    write_def(directory, "coulombintra.def", "NCoulombIntra", 4,
              ["{} 4.0".format(i) for i in range(4)])
    write_def(directory, "gutzwilleridx.def", "NGutzwillerIdx", gutz_count,
              ["{} {}".format(i, i % gutz_count) for i in range(4)] +
              ["{} 1".format(i) for i in range(gutz_count)], complex_type)
    write_def(directory, "orbitalidx.def", "NOrbitalIdx", 16,
              ["{} {} {} 1".format(i, j, 4*i+j)
               for i in range(4) for j in range(4)] +
              ["{} 1".format(i) for i in range(16)], complex_type)
    write_def(directory, "qptransidx.def", "NQPTrans", 1,
              ["0 1.0 0.0"] + ["0 {0} {0} 1".format(i) for i in range(4)])
    # Use nonzero, distinct components throughout, including beyond the old
    # truncated prefix. The final orbital matrix is not a separable rank-one form.
    parameters = [(0.05 + 0.03*i, complex_type*(0.02 + 0.01*i))
                  for i in range(gutz_count)]
    parameters += [(0.5 + 0.13*i + 0.17*j + 0.09*((i*j+2) % 5),
                    complex_type*0.03*(i+2*j+1))
                   for i in range(4) for j in range(4)]
    fields = ["0"] * 6
    for real, imag in parameters:
        fields += ["{:.17e}".format(real), "{:.17e}".format(imag), "0"]
    (directory / "initial.def").write_text(" ".join(fields) + "\n")


def run_vmc(binary, directory, binary_output):
    command = [str(binary)] + (["-b"] if binary_output else [])
    command += ["-e", "namelist.def", "initial.def"]
    if os.environ.get("MVMC_MPI_PROCS"):
        command = [os.environ.get("MVMC_MPIEXEC", "mpiexec"),
                   os.environ.get("MVMC_MPIEXEC_NUMPROC_FLAG", "-n"),
                   os.environ["MVMC_MPI_PROCS"]] + command
    environment = os.environ.copy()
    environment.setdefault("OMP_NUM_THREADS", "1")
    result = subprocess.run(command, cwd=directory, env=environment,
                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                            text=True, timeout=120)
    (directory / "run.log").write_text(result.stdout)
    if result.returncode:
        raise AssertionError("vmc.out failed ({}):\n{}".format(
            result.returncode, result.stdout))


def check_output(text_dir, binary_dir, npara, mode, complex_type):
    files = [7] if mode == 0 else [7, 8]
    steps = 3 if mode == 0 else 1
    expected_names = {"zvo_varbin_{:03d}.dat".format(i) for i in files}
    actual_names = {p.name for p in binary_dir.glob("zvo_varbin_*.dat")}
    if actual_names != expected_names:
        raise AssertionError("binary file indices: {} != {}".format(
            actual_names, expected_names))
    for index in files:
        suffix = "_{:03d}.dat".format(index)
        rows = [[float(s) for s in line.split()]
                for line in (text_dir / ("zvo_var" + suffix)).read_text().splitlines()]
        if len(rows) != steps or any(len(row) != 6+3*npara for row in rows):
            raise AssertionError("unexpected text parameter record dimensions")
        values = []
        for row in rows:
            for parameter in range(npara):
                values.extend(row[6+3*parameter:8+3*parameter])
        if not all(math.isfinite(value) for value in values):
            raise AssertionError("non-finite reference parameters")
        if complex_type and not any(value != 0.0 for value in values[1::2]):
            raise AssertionError("complex fixture has no imaginary parameters")
        if mode == 0 and rows[0][6:] == rows[-1][6:]:
            raise AssertionError("optimization fixture did not update parameters")
        raw = (binary_dir / ("zvo_varbin" + suffix)).read_bytes()
        header_size = struct.calcsize("@ii")
        header = struct.unpack("@ii", raw[:header_size])
        if header != (npara, steps):
            raise AssertionError("binary header {} != {}".format(header, (npara, steps)))
        expected_payload = struct.pack("={}d".format(len(values)), *values)
        if len(raw) != header_size + len(expected_payload):
            raise AssertionError("{}: {} bytes, expected {} (NPara={}, records={})".format(
                suffix, len(raw), header_size+len(expected_payload), npara, steps))
        if raw[header_size:] != expected_payload:
            raise AssertionError("binary parameter values differ from text output")
        if (text_dir / ("zvo_out" + suffix)).read_bytes() != (
                binary_dir / ("zvo_out" + suffix)).read_bytes():
            raise AssertionError("-b changed physical output")
    if mode == 0 and (text_dir / "zqp_opt.dat").read_bytes() != (
            binary_dir / "zqp_opt.dat").read_bytes():
        raise AssertionError("-b changed optimized parameters")


def main():
    if len(sys.argv) != 2:
        raise SystemExit("usage: runtest_varbin.py <vmc.out>")
    binary = Path(sys.argv[1]).resolve()
    work = Path.cwd() / "work"
    work.mkdir(exist_ok=True)
    root = Path(tempfile.mkdtemp(prefix="VarbinOutput_", dir=work))
    print("Varbin evidence: {}".format(root), flush=True)
    failures = []
    for complex_type in (0, 1):
        for gutz_count in (1, 2):
            for mode in (0, 1):
                name = "{}_{}_{}".format(
                    "complex" if complex_type else "real",
                    "odd" if gutz_count == 1 else "even",
                    "opt" if mode == 0 else "physcal")
                try:
                    directories = []
                    for binary_output in (False, True):
                        directory = root / (name + ("_binary" if binary_output else "_text"))
                        prepare_case(directory, complex_type, gutz_count, mode)
                        run_vmc(binary, directory, binary_output)
                        directories.append(directory / "output")
                    check_output(*directories, 16+gutz_count, mode, complex_type)
                    print("{}: PASS".format(name), flush=True)
                except (AssertionError, OSError, ValueError, struct.error,
                        subprocess.TimeoutExpired) as error:
                    failures.append(name)
                    print("{}: FAIL: {}".format(name, error), flush=True)
    print("Varbin output: {}/8 cases passed".format(8-len(failures)))
    return int(bool(failures))


if __name__ == "__main__":
    sys.exit(main())
