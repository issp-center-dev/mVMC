"""Parse the stable end-of-run grand-canonical component statistics line."""
import re


FIELDS = (
    "thread_capacity", "rebuild", "to_pfonly", "recover",
    "probe_rebuild", "debug_near_singular", "direct_value",
    "cofactor_gradient", "composite_zero",
)
PATTERN = re.compile(
    r"^GC component stats:" +
    "".join(r" {}=(\d+)".format(field) for field in FIELDS) + r"$")


def parse_component_stats(output):
    lines = [line for line in output.splitlines()
             if line.startswith("GC component stats:")]
    if len(lines) != 1:
        raise AssertionError(
            "expected exactly one GC component stats line, found {}".format(
                len(lines)))
    match = PATTERN.match(lines[0])
    if match is None:
        raise AssertionError("invalid GC component stats line: {}".format(
            lines[0]))
    return dict(zip(FIELDS, (int(value) for value in match.groups())))
