from __future__ import print_function

import copy
import os
import sys

import validate_bounded_krylov_p4s9_hardened_reparse_policy as validator


def expect_failure(function, label):
    try:
        function()
    except AssertionError:
        return
    raise AssertionError("negative hardened policy accepted: {}".format(
        label))


def main():
    if len(sys.argv) != 3:
        raise AssertionError(
            "usage: runtest_bounded_krylov_p4s9_hardened_reparse_policy.py "
            "POLICY SOURCE_ROOT")
    policy_path = os.path.abspath(sys.argv[1])
    source_root = os.path.abspath(sys.argv[2])
    policy = validator.p4c.read_json(policy_path)
    validator.validate_policy(policy, source_root)

    changed_execution = copy.deepcopy(policy)
    changed_execution["execution_contract"]["threshold_changes"] = True
    expect_failure(
        lambda: validator.validate_policy(changed_execution, source_root),
        "execution threshold change")

    changed_source = copy.deepcopy(policy)
    changed_source["execution_contract"]["source_commit"] = "f" * 40
    expect_failure(
        lambda: validator.validate_policy(changed_source, source_root),
        "source commit mutation")

    changed_package = copy.deepcopy(policy)
    changed_package["execution_contract"][
        "transfer_package_sha256"] = "f" * 64
    expect_failure(
        lambda: validator.validate_policy(changed_package, source_root),
        "transfer package mutation")

    changed_manifest = copy.deepcopy(policy)
    parser_path = next(iter(changed_manifest["parser_manifest_sha256"]))
    changed_manifest["parser_manifest_sha256"][parser_path] = "0" * 64
    expect_failure(
        lambda: validator.validate_policy(changed_manifest, source_root),
        "parser manifest mutation")

    changed_promotion = copy.deepcopy(policy)
    changed_promotion["promotion"][
        "hardened_v13_reparse_closure_required"] = "STOP"
    expect_failure(
        lambda: validator.validate_policy(changed_promotion, source_root),
        "promotion mutation")
    print("bounded Krylov P4-S9 hardened reparse policy: PASS")


if __name__ == "__main__":
    main()
