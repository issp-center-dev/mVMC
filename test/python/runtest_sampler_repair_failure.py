"""Require the intended failure; unrelated errors and signals must not pass."""
import subprocess
import sys

expected = {
    'invalid': ['accept only 0 or 1'],
    'unsupported': ['enabled outside'],
    'factor-failure': ['component factorization failed', 'sampler repair failure'],
    'rollback-fails': ['action=abort-rollback'],
    'rebuild-fails': ['action=abort-proposal-rebuild'],
    'rollback-factor-fails': ['action=abort-rollback', 'info=9'],
    'rebuild-nan': ['action=abort-proposal-rebuild'],
    'current-zero': ['action=abort-initial'],
    'current-factor-fails': ['action=abort-full-recompute', 'info=8'],
}
binary, mode = sys.argv[1:]
proc = subprocess.run([binary, mode], stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                      universal_newlines=True, timeout=30)
print(proc.stdout, end='')
print(proc.stderr, end='', file=sys.stderr)
if proc.returncode != 1 or not all(message in proc.stderr for message in expected[mode]):
    raise SystemExit('unexpected failure for {}: exit={}'.format(mode, proc.returncode))
print('expected sampler failure {} PASS'.format(mode))
