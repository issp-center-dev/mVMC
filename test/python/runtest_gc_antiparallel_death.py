"""Require the intended anti-parallel GC sampler abort; signals do not pass."""
import subprocess
import sys
from runtest_gc import mpi_command
from gc_antiparallel_failure import reject_signal_reports, require_expected_failure

EXPECTED = {
    '--spin-sector-init': ['NGCInitNelec=2', 'exactly antisymmetric',
                           'multiple of 4'],
    '--spin-sector-rank': ['rank-inconsistent GC spin support'],
    '--antiparallel-death-init': ['zero projected amplitude',
                                  'anti-parallel mode', 'Ncur=2'],
    '--antiparallel-death-nonfinite': ['GC sampler met a nonfinite',
                                       'anti-parallel mode'],
    '--antiparallel-death-oldlog': ['GC sampler met a nonfinite',
                                    'anti-parallel mode'],
    '--antiparallel-death-burn': ['GC anti-parallel', 'restored',
                                  'Sz=0'],
    '--death-vacuum_k': ['zero projected amplitude', 'after 1 attempt',
                         'general mode', 'Ncur=0'],
    '--death-general_nonfinite': ['GC sampler met a nonfinite',
                                  'general mode'],
    '--death-cache_uninitialized': ['component cache is not initialized'],
    '--death-rebuild_zero': ['zero projected amplitude',
                             'after periodic rebuild'],
}

binary, mode = sys.argv[1:3]
command = [binary, mode]
if len(sys.argv) > 3:
    command = mpi_command(int(sys.argv[3]), command)
try:
    proc = subprocess.run(command, stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE, universal_newlines=True,
                          timeout=30)
except subprocess.TimeoutExpired:
    raise SystemExit('{} timed out'.format(mode))
print(proc.stdout, end='')
print(proc.stderr, end='', file=sys.stderr)
reject_signal_reports(proc.stdout)
require_expected_failure(proc.returncode, proc.stderr, EXPECTED[mode])
print('expected anti-parallel GC abort {} PASS'.format(mode))
