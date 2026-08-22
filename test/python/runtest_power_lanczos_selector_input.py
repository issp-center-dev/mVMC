from __future__ import print_function

import os
import shutil
import subprocess
import sys


BASE_CASE = "HubbardChain_SpinJastrow"
P6_KEYWORDS = (
    "NLanczosEstimatorMode",
    "NLanczosCoeffWarmUp",
    "NLanczosCoeffSample",
    "NLanczosCoeffInterval",
    "NLanczosFinalWarmUp",
    "NLanczosFinalSample",
    "NLanczosFinalInterval",
    "NLanczosGuideMode",
    "NLanczosStatMode",
)
LEGACY_WARNING = (
    "WARNING: explicit legacy base-support power-Lanczos estimator; "
    "biased diagnostic only; not a corrected release result."
)
LEGACY_GUIDANCE = (
    'P6 LEGACY OPT-IN: add "NLanczosEstimatorMode 1" to ModPara '
    "to reproduce the previous biased base-support estimator; legacy "
    "output is diagnostic-only and is not a corrected release result."
)


def replace_modpara(path, updates, appended=()):
    found = set()
    lines = []
    with open(path) as source:
        for line in source:
            words = line.split()
            if words and words[0] in updates:
                lines.append("{:<27} {}\n".format(words[0], updates[words[0]]))
                found.add(words[0])
            else:
                lines.append(line)
    for key, value in updates.items():
        if key not in found:
            lines.append("{:<27} {}\n".format(key, value))
    for keyword, value in appended:
        lines.append("{:<27} {}\n".format(keyword, value))
    with open(path, "w") as destination:
        destination.writelines(lines)


def write_observable(path, label, rows):
    with open(path, "w") as destination:
        destination.write("====================\n")
        destination.write("{} {}\n".format(label, len(rows)))
        destination.write("====================\n")
        destination.write("observable rows\n")
        destination.write("====================\n")
        for row in rows:
            destination.write(row + "\n")


def prepare_case(rootdir, label):
    source = os.path.join(rootdir, "data", BASE_CASE)
    workdir = os.path.join(rootdir, "work", label)
    if os.path.exists(workdir):
        shutil.rmtree(workdir)
    shutil.copytree(source, workdir)
    output = os.path.join(workdir, "output")
    if os.path.exists(output):
        shutil.rmtree(output)
    return workdir


def run_case(rootdir, workdir, mpi_procs):
    executable = os.path.join(rootdir, "..", "..", "src", "mVMC", "vmc.out")
    command = [executable, "-e", "namelist.def"]
    if mpi_procs:
        command = ["mpirun", "-np", mpi_procs] + command
    return subprocess.run(
        command,
        cwd=workdir,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    )


def assert_rejected(process, expected, label, expect_legacy_guidance=False):
    if process.returncode != 20:
        raise AssertionError(
            "{} exit={} expected=20\n{}".format(
                label, process.returncode, process.stdout
            )
        )
    if process.stdout.count(expected) != 1:
        raise AssertionError(
            "{} expected one {!r}\n{}".format(label, expected, process.stdout)
        )
    expected_guidance_count = 1 if expect_legacy_guidance else 0
    if process.stdout.count(LEGACY_GUIDANCE) != expected_guidance_count:
        raise AssertionError(
            "{} legacy guidance count={} expected={}\n{}".format(
                label,
                process.stdout.count(LEGACY_GUIDANCE),
                expected_guidance_count,
                process.stdout,
            )
        )
    if (
        "Start: Set memories." in process.stdout
        or "Start: Sampling." in process.stdout
    ):
        raise AssertionError("{} allocated or sampled before reject".format(label))


def test_corrected_unavailable(rootdir, mpi_procs):
    workdir = prepare_case(rootdir, "P6Runtime_CorrectedUnavailable")
    replace_modpara(
        os.path.join(workdir, "modpara.def"),
        {"NLanczosMode": "1", "NVMCCalMode": "1"},
    )
    process = run_case(rootdir, workdir, mpi_procs)
    assert_rejected(
        process,
        "P6 INPUT REJECTED: CORRECTED_PIPELINE_UNAVAILABLE",
        "corrected default",
        expect_legacy_guidance=True,
    )


def test_observable_certificate(rootdir, mpi_procs):
    workdir = prepare_case(rootdir, "P6Runtime_ObservableCertificate")
    replace_modpara(
        os.path.join(workdir, "modpara.def"),
        {"NLanczosMode": "2", "NVMCCalMode": "1"},
    )
    write_observable(
        os.path.join(workdir, "greenone.def"), "NCisAjs", ["0 0 0 0"]
    )
    write_observable(
        os.path.join(workdir, "greentwoex.def"),
        "NCisAjsCktAlt",
        ["0 0 0 0 1 1 1 1"],
    )
    write_observable(
        os.path.join(workdir, "greentwo.def"),
        "NCisAjsCktAltDC",
        ["0 0 0 0 1 1 1 1"],
    )
    with open(os.path.join(workdir, "namelist.def"), "a") as destination:
        destination.write("    TwoBodyGEx  greentwoex.def\n")
    process = run_case(rootdir, workdir, mpi_procs)
    assert_rejected(
        process,
        "P6 INPUT REJECTED: OBSERVABLE_CERTIFICATE_UNAVAILABLE",
        "corrected observable certificate",
        expect_legacy_guidance=True,
    )
    if "OBSERVABLE_CENSUS_" in process.stdout:
        raise AssertionError(
            "corrected observable ran census before unavailable gate\n{}".format(
                process.stdout
            )
        )


def test_observable_unavailable_precedes_census_error(rootdir, mpi_procs):
    workdir = prepare_case(rootdir, "P6Runtime_ObservableBeforeCensusError")
    replace_modpara(
        os.path.join(workdir, "modpara.def"),
        {"NLanczosMode": "2", "NVMCCalMode": "1"},
    )
    write_observable(
        os.path.join(workdir, "greenone.def"), "NCisAjs", ["malformed"]
    )
    process = run_case(rootdir, workdir, mpi_procs)
    assert_rejected(
        process,
        "P6 INPUT REJECTED: OBSERVABLE_CERTIFICATE_UNAVAILABLE",
        "corrected observable before census error",
        expect_legacy_guidance=True,
    )
    if "OBSERVABLE_CENSUS_" in process.stdout:
        raise AssertionError(
            "malformed observable reached census before unavailable gate\n{}".format(
                process.stdout
            )
        )


def test_strict_grammar(rootdir):
    values = ("01", "+1", "1.0", "2147483648")
    for ordinal, keyword in enumerate(P6_KEYWORDS):
        value = values[ordinal % len(values)]
        workdir = prepare_case(
            rootdir, "P6Runtime_Grammar_{:02d}".format(ordinal)
        )
        replace_modpara(os.path.join(workdir, "modpara.def"), {keyword: value})
        process = run_case(rootdir, workdir, None)
        assert_rejected(
            process,
            "P6 INPUT REJECTED: INVALID_INTEGER_TOKEN {}".format(keyword),
            "strict grammar {}={}".format(keyword, value),
        )


def test_duplicate(rootdir):
    for ordinal, keyword in enumerate(P6_KEYWORDS):
        workdir = prepare_case(
            rootdir, "P6Runtime_Duplicate_{:02d}".format(ordinal)
        )
        replace_modpara(
            os.path.join(workdir, "modpara.def"),
            {keyword: "0"},
            appended=((keyword, "0"),),
        )
        process = run_case(rootdir, workdir, None)
        assert_rejected(
            process,
            "P6 INPUT REJECTED: DUPLICATE_KEYWORD {}".format(keyword),
            "duplicate keyword {}".format(keyword),
        )


def test_cross_field_rejection(rootdir, mpi_procs):
    cases = (
        (
            "DisabledEstimator",
            {"NLanczosMode": "0", "NLanczosEstimatorMode": "1"},
            "NLanczosMode=0 requires every P6 estimator",
        ),
        (
            "LegacyChain",
            {
                "NLanczosMode": "1",
                "NLanczosEstimatorMode": "1",
                "NLanczosCoeffSample": "1",
            },
            "explicit legacy estimator mode requires all coefficient",
        ),
        (
            "GuideMode",
            {"NLanczosGuideMode": "1"},
            "NLanczosGuideMode only supports integer 0",
        ),
        (
            "StatMode",
            {"NLanczosStatMode": "1"},
            "NLanczosStatMode only supports integer 0",
        ),
    )
    for label, updates, expected in cases:
        workdir = prepare_case(rootdir, "P6Runtime_{}".format(label))
        replace_modpara(os.path.join(workdir, "modpara.def"), updates)
        process = run_case(rootdir, workdir, mpi_procs)
        assert_rejected(process, expected, label)


def test_explicit_legacy_augmented(rootdir, mpi_procs):
    workdir = prepare_case(rootdir, "P6Runtime_ExplicitLegacyAugmented")
    replace_modpara(
        os.path.join(workdir, "modpara.def"),
        {
            "NLanczosMode": "2",
            "NLanczosEstimatorMode": "1",
            "NVMCCalMode": "1",
        },
    )
    process = run_case(rootdir, workdir, mpi_procs)
    if process.returncode != 0:
        raise AssertionError(
            "explicit legacy augmented exit={}\n{}".format(
                process.returncode, process.stdout
            )
        )
    if process.stdout.count(LEGACY_WARNING) != 1:
        raise AssertionError(
            "explicit legacy expected one warning\n{}".format(process.stdout)
        )
    if "P6 INPUT REJECTED:" in process.stdout:
        raise AssertionError("explicit legacy was rejected")


def main():
    rootdir = os.getcwd()
    mpi_procs = os.environ.get("MVMC_MPI_PROCS")
    test_corrected_unavailable(rootdir, mpi_procs)
    test_observable_certificate(rootdir, mpi_procs)
    test_observable_unavailable_precedes_census_error(rootdir, mpi_procs)
    test_cross_field_rejection(rootdir, mpi_procs)
    test_explicit_legacy_augmented(rootdir, mpi_procs)
    if not mpi_procs:
        test_strict_grammar(rootdir)
        test_duplicate(rootdir)
    print("power-Lanczos runtime input: PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
