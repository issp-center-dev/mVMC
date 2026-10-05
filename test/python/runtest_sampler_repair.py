"""Check the real solver's repair logging without a long sampling benchmark."""
import collections
import hashlib
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile

binary, dry, fixture, scratch = map(lambda v: Path(v).resolve(), sys.argv[1:5])
launcher = sys.argv[5:]
rank_count = int(launcher[-1]) if launcher else 1
scratch.mkdir(parents=True, exist_ok=True)
env = {k: v for k, v in os.environ.items() if not k.startswith('MVMC_SAMPLER_')}
env['OMP_NUM_THREADS'] = '1'

def run(cmd, cwd, settings, expected_error=None):
    proc = subprocess.run(cmd, cwd=str(cwd), env=dict(env, **settings),
                          stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                          universal_newlines=True, timeout=90)
    (cwd/'run.log').write_text(proc.stdout)
    if expected_error:
        assert proc.returncode != 0 and expected_error in proc.stdout, proc.stdout
    elif proc.returncode:
        raise RuntimeError(proc.stdout)
    return proc.stdout

def physical(d):
    return {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
            for p in (d/'output').glob('*.dat')
            if p.name.startswith(('zvo_out_', 'zvo_var_', 'zvo_ls_'))}

with tempfile.TemporaryDirectory(prefix='sampler-repair-', dir=str(scratch)) as tmp:
    work = Path(tmp)
    generated = work/'input'
    generated.mkdir()
    shutil.copy2(str(fixture/'StdFace.def'), str(generated/'StdFace.def'))
    shutil.copy2(str(fixture/'zqp_opt.dat'), str(generated/'zqp_opt.dat'))
    run([str(dry), 'StdFace.def'], generated, {})
    p = generated/'modpara.def'
    text = p.read_text()
    for key, value in [('NVMCCalMode', 1), ('NLanczosMode', 1), ('NLanczosStep', 1),
                       ('NLanczosEstimatorMode', 0), ('NExUpdatePath', 1),
                       ('NDataQtySmp', 2), ('NVMCSample', 64), ('NVMCWarmUp', 8),
                       ('NSplitSize', 1)]:
        text, count = re.subn(r'^'+key+r'\s+\S+', key+' '+str(value), text, flags=re.M)
        if not count:
            text += '\n'+key+' '+str(value)+'\n'
    p.write_text(text)
    for name, settings in [('logged', {}), ('quiet', {'MVMC_SAMPLER_DRIFT_LOG': '0'})]:
        d = work/name
        shutil.copytree(str(generated), str(d))
        out = run(launcher+[str(binary), '-e', 'namelist.def', 'zqp_opt.dat'], d, settings)
        assert 'Sampler component repair: on;' in out
    assert physical(work/'logged') == physical(work/'quiet')
    assert len(physical(work/'logged')) >= 10
    assert not list((work/'quiet').glob('sampler_drift_r*.dat'))
    logs = list((work/'logged').glob('sampler_drift_r*.dat'))
    assert len(logs) == rank_count
    for p in logs:
        rows = [l.split() for l in p.read_text().splitlines() if l and l[0] != '#']
        q = [r for r in rows if r[0] == 'Q']
        z = [r for r in rows if r[0] == 'Z']
        recovery = [r for r in rows if r[0] == 'R']
        assert len(recovery) == 2 and all(len(r) == 6 for r in recovery)
        assert all(all(int(v) == 0 for v in r[2:]) for r in recovery)
        events = [r for r in rows if r[0] == 'E']
        assert len(q) == len(z) == 2
        assert [int(r[1]) for r in q] == [0, 1]
        assert all(len(r) == 11 and int(r[2]) == 64 and int(r[3]) == 0 for r in q)
        assert all(len(r) == 16 for r in z)
        assert len(events) <= 16
        assert all(n <= 4 for n in collections.Counter(r[5] for r in events).values())
        assert p.stat().st_size < 16384
    # A positive guide epsilon overlaps the supported real/split1/hop-exchange
    # conditions above, but uses a different sampling weight and must stay off.
    guide_input = work/'guide-input'
    shutil.copytree(str(generated), str(guide_input))
    p = guide_input/'modpara.def'
    text, count = re.subn(r'^DLanczosGuideEps\s+\S+', 'DLanczosGuideEps 1.0',
                          p.read_text(), flags=re.M)
    if not count:
        text += '\nDLanczosGuideEps 1.0\n'
    p.write_text(text)
    for name, settings in [('guide-auto', {}), ('guide-off', {'MVMC_SAMPLER_REPAIR': '0'}),
                           ('guide-force', {'MVMC_SAMPLER_REPAIR': '1'}),
                           ('guide-log', {'MVMC_SAMPLER_DRIFT_LOG': '1'})]:
        d = work/name
        shutil.copytree(str(guide_input), str(d))
        error = 'enabled outside' if name == 'guide-force' else (
            'drift logging requested outside' if name == 'guide-log' else None)
        out = run(launcher+[str(binary), '-e', 'namelist.def', 'zqp_opt.dat'],
                  d, settings, error)
        assert not list(d.glob('sampler_drift_r*.dat'))
        if error:
            assert not physical(d)
        else:
            assert 'Sampler component repair: off;' in out
    assert physical(work/'guide-auto') == physical(work/'guide-off')
    assert len(physical(work/'guide-auto')) >= 10
    guide_stats = list((work/'guide-auto'/'output').glob('zvo_nzguide_*.dat'))
    assert len(guide_stats) == 2
    for p in guide_stats:
        assert p.read_bytes() == (work/'guide-off'/'output'/p.name).read_bytes()
print('sampler repair solver/logging smoke PASS')
