"""Narrow read-only identity proof for the pre-reviewed timedatex system service.

No process-name-only exemption, arbitrary unit enrollment, or network calls.
The service must match its preflight-pinned executable and systemd unit bytes.
"""
import hashlib
import os
from pathlib import Path
import stat
import subprocess

UNIT = 'timedatex.service'
COMM = 'timedatex'
EXE = '/usr/sbin/timedatex'
CGROUP = '/system.slice/timedatex.service'


def secure_file(path):
    p = Path(path)
    if not p.is_absolute() or p.resolve() != p:
        raise RuntimeError('System service file is redirected')
    before = p.stat()
    if before.st_uid != 0 or before.st_gid != 0 or before.st_mode & 0o022 or not stat.S_ISREG(before.st_mode):
        raise RuntimeError('System service file ownership or permissions unsafe')
    digest = hashlib.sha256(p.read_bytes()).hexdigest()
    after = p.stat()
    if (before.st_ino, before.st_dev, before.st_size, before.st_mtime_ns) != (after.st_ino, after.st_dev, after.st_size, after.st_mtime_ns):
        raise RuntimeError('System service file changed while hashing')
    return dict(path=str(p), sha256=digest, bytes=after.st_size, root_owned_nonwritable=True)


def unit_properties(run_fn=subprocess.run):
    cmd = ['systemctl', 'show', UNIT, '--property=Id', '--property=MainPID',
           '--property=FragmentPath', '--property=ControlGroup', '--property=LoadState']
    res = run_fn(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=3)
    if res.returncode:
        raise RuntimeError('Cannot verify timedatex service unit')
    values = dict(line.split('=', 1) for line in res.stdout.decode('utf-8', 'strict').splitlines() if '=' in line)
    if values.get('Id') != UNIT or values.get('LoadState') != 'loaded':
        raise RuntimeError('Unexpected timedatex systemd identity')
    fragment = Path(values['FragmentPath']).resolve()
    if not any(str(fragment).startswith(prefix) for prefix in ('/usr/lib/systemd/system/', '/etc/systemd/system/')):
        raise RuntimeError('Unexpected timedatex unit fragment location')
    values['fragment'] = secure_file(fragment)
    return values


def enroll(run_fn=subprocess.run):
    """Called during operator-reviewed preflight, then sealed in the baseline."""
    unit = unit_properties(run_fn)
    return dict(role='pinned_restarting_system_service', comm=COMM, unit=UNIT,
                executable=secure_file(EXE), fragment=unit['fragment'], cgroup=CGROUP,
                reason='timedatex was the sole previous false blocker after PID turnover; require root-owned executable/unit hashes and live systemd MainPID/cgroup proof')


def process_proof(entry, row, run_fn=subprocess.run):
    """Collect only this fixed service's safe identity fields; fail closed."""
    if row.get('comm') != COMM:
        return None
    try:
        uids = None
        for line in (entry / 'status').read_text().splitlines():
            if line.startswith('Uid:'):
                uids = [int(x) for x in line.split()[1:]]
        exe = os.readlink(str(entry / 'exe'))
        if exe != EXE:
            raise RuntimeError('Service executable mismatch')
        groups = [line.split(':', 2)[2] for line in (entry / 'cgroup').read_text().splitlines()]
        unit = unit_properties(run_fn)
        if int(unit['MainPID']) != row['pid'] or unit.get('ControlGroup') != CGROUP or CGROUP not in groups:
            raise RuntimeError('Service process is not systemd main process')
        if row['ppid'] != 1 or uids != [0, 0, 0, 0]:
            raise RuntimeError('Service process ownership or ancestry mismatch')
        return dict(status='VERIFIED', role='pinned_restarting_system_service', comm=COMM, unit=UNIT,
                    executable=secure_file(EXE), fragment=unit['fragment'], cgroup=CGROUP,
                    pid=row['pid'], start_ticks=row['start_ticks'], main_pid=int(unit['MainPID']), uids=uids)
    except (OSError, ValueError, KeyError, RuntimeError, subprocess.SubprocessError) as error:
        return dict(status='UNVERIFIED', error_type=type(error).__name__)


def valid_profile(profile):
    if not isinstance(profile, dict):
        return False
    if (profile.get('role'), profile.get('comm'), profile.get('unit'), profile.get('cgroup')) != (
            'pinned_restarting_system_service', COMM, UNIT, CGROUP):
        return False
    if not profile.get('reason'):
        return False
    exe, unit = profile.get('executable', {}), profile.get('fragment', {})
    if exe.get('path') != EXE or not any(unit.get('path', '').startswith(p) for p in ('/usr/lib/systemd/system/', '/etc/systemd/system/')):
        return False
    for f in (exe, unit):
        if not f.get('root_owned_nonwritable') or not isinstance(f.get('bytes'), int) or f['bytes'] <= 0:
            return False
        digest=f.get('sha256','')
        if len(digest)!=64 or any(c not in '0123456789abcdef' for c in digest):
            return False
    return True


def matches(row, profiles):
    proof = row.get('pinned_system_service_proof', {}) or {}
    if proof.get('status') != 'VERIFIED' or row.get('comm') != COMM or row.get('ppid') != 1:
        return False
    if proof.get('pid') != row.get('pid') or proof.get('start_ticks') != row.get('start_ticks') or proof.get('main_pid') != row.get('pid') or proof.get('uids') != [0,0,0,0]:
        return False
    keys=('role','comm','unit','cgroup','executable','fragment')
    return any(valid_profile(p) and all(proof.get(k)==p.get(k) for k in keys) for p in profiles)
