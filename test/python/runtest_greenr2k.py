"""Tests for the post-processing tool greenr2k (tool/greenr2k.F90).

usage: runtest_greenr2k.py [<case> ...]

vmc.out is run once (physics mode, a 4-site Hubbard chain from StdFace),
and greenr2k is run on copies of that directory. All cases are run when
none is given; "normal" is always run first because it is the reference.

cases:
  normal         greenr2k exits 0 and writes output/zvo_corr.dat
  slash_paths    the namelist names files with '/' ("sub/greenone.def",
                 "./modpara.def", one separated by a tab); the result is
                 byte for byte the same as "normal"
  quoted_long_paths  quoted paths, embedded spaces and quotes, long records
  path_too_long  reject a name exceeding the tool's 256-character field
  invalid_dimension  reject invalid FullDiag headers before writing output
  missing_index  greentwo.def lacks required entries -> nonzero exit
                 status and the message "Missing indices ..."
  missing_file   the namelist names a file that does not exist -> nonzero
                 exit status, an error message, and no file is created
  fulldiag       an HPhi FullDiag layout (CalcMod with CalcType 2,
                 output/CHECK_Memory.dat, *_eigen*.dat) -> greenr2k exits 0
                 and writes correlations for two eigenstates
"""

from __future__ import print_function

import os
import shutil
import subprocess
import sys

TIMEOUT = 300
HEAD = "zvo"
CASES = ["normal", "slash_paths", "quoted_long_paths", "path_too_long",
         "missing_index", "missing_file", "fulldiag", "invalid_dimension"]

STDFACE = """\
L = 4
Lsub = 2
model = "Hubbard"
lattice = "chain"
U = 4.0
t = 1.0
nelec = 4
2Sz = 0
NVMCCalMode = 1
NVMCSample = 10
"""
KPATH = "2 4\nG 0 0 0\nX 0.5 0 0\n4 1 1\n"
# The first line is what HPhi writes to output/CHECK_Memory.dat.
CHECK_MEMORY = ("  MAX DIMENSION idim_max=2 \n"
                "  APPROXIMATE REQUIRED MEMORY  max_mem=0.000001 GB \n")


def run_command(cmd, cwd):
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = "1"
    proc = subprocess.Popen(cmd, cwd=cwd, env=env,
                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    try:
        out, _ = proc.communicate(timeout=TIMEOUT)
    except subprocess.TimeoutExpired:
        proc.kill()
        out, _ = proc.communicate()
        print(out.decode("utf-8", "replace"))
        print("FAIL: {} did not finish within {} s".format(cmd[0], TIMEOUT))
        return None, ""
    out = out.decode("utf-8", "replace")
    print(out)
    return proc.returncode, out


def prepare_base(build, basedir):
    """Run vmcdry.out and vmc.out once; the cases copy this directory."""
    if os.path.exists(basedir):
        shutil.rmtree(basedir)
    os.makedirs(basedir)
    with open(os.path.join(basedir, "StdFace.def"), "w") as stream:
        stream.write(STDFACE)
    code, _ = run_command(
        [os.path.join(build, "src", "mVMC", "vmcdry.out"), "StdFace.def"],
        basedir)
    if code != 0:
        print("FAIL: vmcdry.out exited with {}".format(code))
        return False
    with open(os.path.join(basedir, "geometry.dat"), "a") as stream:
        stream.write(KPATH)
    code, _ = run_command(
        [os.path.join(build, "src", "mVMC", "vmc.out"), "-e", "namelist.def"],
        basedir)
    if code != 0:
        print("FAIL: vmc.out exited with {}".format(code))
        return False
    for name in ["_cisajs_001.dat", "_cisajscktalt_001.dat"]:
        if not os.path.isfile(os.path.join(basedir, "output", HEAD + name)):
            print("FAIL: vmc.out did not write {}".format(HEAD + name))
            return False
    return True


def prepare_case(basedir, casedir):
    if os.path.exists(casedir):
        shutil.rmtree(casedir)
    shutil.copytree(basedir, casedir)


def run_greenr2k(build, casedir, namelist="namelist.def"):
    return run_command(
        [os.path.join(build, "tool", "greenr2k"), namelist, "geometry.dat"],
        casedir)


def read_bytes(path):
    with open(path, "rb") as stream:
        return stream.read()


def rewrite_namelist(casedir, replace, separator="  ", indent=""):
    """Replace the file names of the given keywords in namelist.def."""
    path = os.path.join(casedir, "namelist.def")
    with open(path) as stream:
        lines = stream.readlines()
    result = []
    for line in lines:
        words = line.split()
        if len(words) == 2 and words[0] in replace:
            line = "{}{}{}{}\n".format(indent, words[0], separator, replace[words[0]])
        result.append(line)
    with open(path, "w") as stream:
        stream.writelines(result)


def read_index_file(path):
    """Header lines and the index rows (lists of int) of a Green-function
    index file."""
    with open(path) as stream:
        lines = stream.readlines()
    header, rows = lines[:5], []
    for line in lines[5:]:
        words = line.split()
        if words:
            rows.append([int(word) for word in words])
    return header, rows


def write_index_file(path, header, rows, keyword):
    header = list(header)
    header[1] = "{} {}\n".format(keyword, len(rows))
    with open(path, "w") as stream:
        stream.writelines(header)
        for row in rows:
            stream.write("".join("{:6d}".format(i) for i in row) + "\n")


def hphi_spin_rows(rows):
    """S+S- and S-S+ entries in the convention of HPhi.

    mVMC writes <c+_{i s} c_{j s} c+_{j t} c_{i t}> (i s j s j t i t);
    HPhi writes  <c+_{i s} c_{i t} c+_{j t} c_{j s}> (i s i t j t j s),
    which is what greenr2k expects when CalcType is not mVMC.
    """
    result = []
    for row in rows:
        i, s, j, t = row[0], row[1], row[2], row[5]
        if (row[6], row[4], row[3], row[7]) == (i, j, s, t) and s != t:
            new = [i, s, i, t, j, t, j, s]
            if new not in result:
                result.append(new)
    return result


def case_normal(build, basedir, casedir):
    prepare_case(basedir, casedir)
    code, _ = run_greenr2k(build, casedir)
    if code != 0:
        print("FAIL: greenr2k exited with {}".format(code))
        return 1
    for name in [os.path.join("output", HEAD + "_corr.dat"),
                 os.path.join("output", HEAD + "_corr_001.dat.frmsf"),
                 "kpath.gp"]:
        if not os.path.isfile(os.path.join(casedir, name)):
            print("FAIL: greenr2k did not write {}".format(name))
            return 1
    return 0


def case_slash_paths(build, basedir, casedir, refdir):
    prepare_case(basedir, casedir)
    os.mkdir(os.path.join(casedir, "sub"))
    for name in ["greenone.def", "greentwo.def"]:
        shutil.move(os.path.join(casedir, name),
                    os.path.join(casedir, "sub", name))
    rewrite_namelist(casedir, {"ModPara": "./modpara.def",
                               "OneBodyG": "sub/greenone.def"})
    rewrite_namelist(casedir, {"TwoBodyG": "sub/greentwo.def"},
                     separator="\t")
    code, out = run_greenr2k(build, casedir)
    if code != 0:
        print("FAIL: greenr2k exited with {}".format(code))
        return 1
    for name in ["./modpara.def", "sub/greenone.def", "sub/greentwo.def"]:
        if "Read from " + name not in out:
            print("FAIL: greenr2k did not read {}".format(name))
            return 1
    corr = os.path.join("output", HEAD + "_corr.dat")
    if not os.path.isfile(os.path.join(casedir, corr)):
        print("FAIL: greenr2k did not write {}".format(corr))
        return 1
    if read_bytes(os.path.join(casedir, corr)) != \
            read_bytes(os.path.join(refdir, corr)):
        print("FAIL: {} differs from the normal case".format(corr))
        return 1
    return 0


def case_quoted_long_paths(build, basedir, casedir, refdir):
    prepare_case(basedir, casedir)
    # A 256-character relative path split over two filesystem components.
    os.mkdir(os.path.join(casedir, "subdir"))
    long_name = "subdir/" + "f" * 245 + ".def"
    shutil.move(os.path.join(casedir, "greenone.def"), os.path.join(casedir, long_name))
    shutil.move(os.path.join(casedir, "greentwo.def"), os.path.join(casedir, "two's body.def"))
    rewrite_namelist(casedir, {"OneBodyG": long_name})
    rewrite_namelist(casedir, {"ModPara": '"./modpara.def"', "TwoBodyG": "'two''s body.def'"},
                     separator="\t", indent=" " * 600)
    code, _ = run_greenr2k(build, casedir)
    if code != 0:
        return 1
    for name in [HEAD + "_corr.dat", HEAD + "_corr_001.dat.frmsf"]:
        if read_bytes(os.path.join(casedir, "output", name)) != read_bytes(os.path.join(refdir, "output", name)):
            print("FAIL: {} differs for quoted or long paths".format(name))
            return 1
    return 0


def case_path_too_long(build, basedir, casedir):
    prepare_case(basedir, casedir)
    name = "subdir/" + "x" * 250
    for value in [name, '"' + name + '"']:
        rewrite_namelist(casedir, {"OneBodyG": value})
        code, out = run_greenr2k(build, casedir)
        if code in (None, 0) or "File name too long." not in out:
            print("FAIL: expected a length error for a 257-character file name")
            return 1
    return 0


def case_missing_index(build, basedir, casedir):
    prepare_case(basedir, casedir)
    path = os.path.join(casedir, "greentwo.def")
    header, rows = read_index_file(path)
    write_index_file(path, header, rows[:1], "NCisAjsCktAltDC")
    code, out = run_greenr2k(build, casedir)
    if code is None:
        return 1
    if code == 0:
        print("FAIL: greenr2k exited with 0 although indices are missing")
        return 1
    if "Missing indices for the Green function." not in out:
        print("FAIL: the message about missing indices was not printed")
        return 1
    if os.path.exists(os.path.join(casedir, "output", HEAD + "_corr.dat")):
        print("FAIL: greenr2k wrote the correlation file")
        return 1
    return 0


def case_missing_file(build, basedir, casedir):
    prepare_case(basedir, casedir)
    os.remove(os.path.join(casedir, "greenone.def"))
    code, out = run_greenr2k(build, casedir)
    if code is None:
        return 1
    if code == 0:
        print("FAIL: greenr2k exited with 0 although greenone.def is missing")
        return 1
    if "Cannot open the OneBodyG file" not in out:
        print("FAIL: the message about the missing file was not printed")
        return 1
    if os.path.exists(os.path.join(casedir, "greenone.def")):
        print("FAIL: greenr2k created the missing input file")
        return 1
    return 0


def case_fulldiag(build, basedir, casedir):
    prepare_case(basedir, casedir)
    with open(os.path.join(casedir, "calcmod.def"), "w") as stream:
        stream.write("CalcType   2\nCalcModel  0\n")
    with open(os.path.join(casedir, "namelist.def"), "a") as stream:
        stream.write("CalcMod  calcmod.def\n")
    output = os.path.join(casedir, "output")
    with open(os.path.join(output, "CHECK_Memory.dat"), "w") as stream:
        stream.write(CHECK_MEMORY)
    # Two-body indices: keep the mVMC entries and add the HPhi-style
    # spin-flip entries, with zero values in the correlation file.
    path = os.path.join(casedir, "greentwo.def")
    header, rows = read_index_file(path)
    extra = hphi_spin_rows(rows)
    if not extra:
        print("FAIL: no S+S- entries were found in greentwo.def")
        return 1
    write_index_file(path, header, rows + extra, "NCisAjsCktAltDC")
    shutil.copy(os.path.join(output, HEAD + "_cisajs_001.dat"),
                os.path.join(output, HEAD + "_cisajs_eigen0.dat"))
    shutil.copy(os.path.join(output, HEAD + "_cisajscktalt_001.dat"),
                os.path.join(output, HEAD + "_cisajscktalt_eigen0.dat"))
    with open(os.path.join(output, HEAD + "_cisajscktalt_eigen0.dat"),
              "a") as stream:
        for row in extra:
            stream.write(" ".join(str(i) for i in row) + "  0.0  0.0\n")
    for prefix in [HEAD + "_cisajs", HEAD + "_cisajscktalt"]:
        shutil.copyfile(os.path.join(output, prefix + "_eigen0.dat"),
                        os.path.join(output, prefix + "_eigen1.dat"))
    code, out = run_greenr2k(build, casedir)
    if code != 0:
        print("FAIL: greenr2k exited with {}".format(code))
        return 1
    for text in ["Method : Full Diagonalization", "Number of States :"]:
        if text not in out:
            print("FAIL: '{}' was not printed".format(text))
            return 1
    for state in (0, 1):
        for suffix in (".dat", ".dat.frmsf"):
            name = HEAD + "_corr_eigen" + str(state) + suffix
            if not os.path.isfile(os.path.join(output, name)):
                print("FAIL: greenr2k did not write {}".format(name))
                return 1
    return 0


def case_invalid_dimension(build, basedir, casedir):
    prepare_case(basedir, casedir)
    with open(os.path.join(casedir, "namelist.def"), "a") as stream:
        stream.write("CalcMod calcmod.def\n")
    with open(os.path.join(casedir, "calcmod.def"), "w") as stream:
        stream.write("CalcType 2\n")
    output = os.path.join(casedir, "output")
    records = ["", "other=1", "MAX DIMENSION idim_max"] + [
        "MAX DIMENSION idim_max=" + value
        for value in ("", "0", "-1", "/", ",", "+", "abc", "1.5", "1 2",
                      "9" * 40, "1" + " " * 512 + "junk")]
    for record in records:
        with open(os.path.join(output, "CHECK_Memory.dat"), "w") as stream:
            stream.write(record + "\n" if record else "")
        code, out = run_greenr2k(build, casedir)
        if code in (None, 0) or "Cannot read idim_max" not in out:
            print("FAIL: expected a dimension error for {!r}".format(record))
            return 1
        if any(name.startswith(HEAD + "_corr") for name in os.listdir(output)):
            print("FAIL: greenr2k wrote correlation output for an invalid dimension")
            return 1
    return 0


def main():
    cases = sys.argv[1:] or CASES
    for case in cases:
        if case not in CASES:
            print(__doc__)
            print("unknown case: {}".format(case))
            return 2
    if "normal" not in cases:
        cases = ["normal"] + cases
    rootdir = os.getcwd()
    build = os.path.join(rootdir, "..", "..")
    workdir = os.path.join(rootdir, "work", "Greenr2k_Tool")
    basedir = os.path.join(workdir, "base")
    if not prepare_base(build, basedir):
        return 1
    refdir = os.path.join(workdir, "normal")
    failed = []
    for case in cases:
        print("===== case: {} =====".format(case))
        casedir = os.path.join(workdir, case)
        if case == "normal":
            code = case_normal(build, basedir, casedir)
        elif case == "slash_paths":
            code = case_slash_paths(build, basedir, casedir, refdir)
        elif case == "quoted_long_paths":
            code = case_quoted_long_paths(build, basedir, casedir, refdir)
        elif case == "path_too_long":
            code = case_path_too_long(build, basedir, casedir)
        elif case == "invalid_dimension":
            code = case_invalid_dimension(build, basedir, casedir)
        elif case == "missing_index":
            code = case_missing_index(build, basedir, casedir)
        elif case == "missing_file":
            code = case_missing_file(build, basedir, casedir)
        else:
            code = case_fulldiag(build, basedir, casedir)
        print("===== case: {} : {} =====".format(
            case, "PASS" if code == 0 else "FAIL"))
        if code != 0:
            failed.append(case)
    if failed:
        print("FAILED cases: {}".format(", ".join(failed)))
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
