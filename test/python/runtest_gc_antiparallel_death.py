"""Require the intended anti-parallel GC sampler abort; signals do not pass."""
import subprocess
import sys
from gc_antiparallel_failure import reject_signal_reports, require_expected_failure

EXPECTED = {
    '--antiparallel-death-init': ['makeInitialSampleGC exceeded 100 attempts',
                                  'anti-parallel', 'Ncur=2'],
    '--antiparallel-death-nonfinite': ['GC anti-parallel sampler',
                                       'nonfinite'],
    '--antiparallel-death-oldlog': ['GC anti-parallel sampler', 'nonfinite'],
    '--antiparallel-death-burn': ['GC anti-parallel', 'restored',
                                  'Sz=0'],
}

binary, mode = sys.argv[1:]
try:
    proc = subprocess.run([binary, mode], stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE, universal_newlines=True,
                          timeout=30)
except subprocess.TimeoutExpired:
    raise SystemExit('{} timed out'.format(mode))
print(proc.stdout, end='')
print(proc.stderr, end='', file=sys.stderr)
reject_signal_reports(proc.stdout)
require_expected_failure(proc.returncode, proc.stderr, EXPECTED[mode])
print('expected anti-parallel GC abort {} PASS'.format(mode))
