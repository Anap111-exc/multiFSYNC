#!/usr/bin/env python3
"""Completion-driven normal StopCharging, independent of a local PC.

No API call occurs on import or --validate-only. All actual API calls use the
existing ECS RAM-role CLI profile, never credential files or printed secrets.
Successful API acceptance is NOT a claim of Stopped or stopped billing.
"""
import argparse
import datetime
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import time

try:
    from . import system_service
except ImportError:
    # Support the documented direct script entrypoint as well as module imports.
    import system_service

EXPERIMENT = 'cavi_true222_fit444_m4_20261117_v1'
INSTANCE = 'i-bp1egmsairexm9w1rnks'
REGION = 'cn-hangzhou'
EX = Path('/root/cavi-true222-fit444-m4-exec-20261117-v1')
OUT = Path('/root/multifsync-results/cavi_true222_fit444_m4_20261117_v1')
STAGE = Path('/root/cavi-true222-fit444-m4-stage-20261117-v1')
PROFILE = 'multifsync-self-stop'
NUMERICAL = re.compile(r'^(?:R|Rscript|julia|matlab|torchrun)$', re.I)
RUNTIME = re.compile(r'^(?:python[0-9.]*|pypy[0-9.]*|node(?:js)?|java)$', re.I)


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def now():
    return datetime.datetime.utcnow().isoformat() + 'Z'


def strict_path(path, root, require_file=True):
    p, base = Path(path), Path(root)
    require(p.is_absolute() and p.resolve() == p and base.resolve() == base,
            'Nonabsolute or redirected lifecycle path')
    try:
        p.relative_to(base)
    except ValueError:
        raise RuntimeError('Lifecycle path outside its designated root')
    require(p != base and (p.is_file() if require_file else True),
            'Missing lifecycle file or root used as a file')
    return p


def read_bound(path, expected_hash, root):
    p = strict_path(path, root)
    require(re.fullmatch('[0-9a-f]{64}', expected_hash or '') is not None,
            'Invalid lifecycle SHA256')
    require(sha(p) == expected_hash, 'Lifecycle input hash mismatch')
    return json.loads(p.read_text(encoding='utf-8-sig'))


def persist(path, value, exclusive=False):
    """fsync both file and directory; never overwrite an initial stop request."""
    path = Path(path)
    temp = path if exclusive else path.with_name(path.name + '.part')
    with temp.open('x', encoding='utf-8') as stream:
        json.dump(value, stream, ensure_ascii=False, indent=2)
        stream.write('\n')
        stream.flush()
        os.fsync(stream.fileno())
    if not exclusive:
        os.replace(str(temp), str(path))
    if hasattr(os, 'O_DIRECTORY'):
        fd = os.open(str(path.parent), os.O_DIRECTORY)
        try:
            os.fsync(fd)
        finally:
            os.close(fd)


def validate_contract(completion_path, completion_hash, baseline_path, baseline_hash):
    complete = read_bound(completion_path, completion_hash, OUT)
    baseline = read_bound(baseline_path, baseline_hash, STAGE)
    require(complete.get('schema') == 'multifsync_lifecycle_completion_v1',
            'Wrong completion schema')
    require(baseline.get('schema') == 'multifsync_process_baseline_v1',
            'Wrong process baseline schema')
    for obj in (complete, baseline):
        require(obj.get('experiment_id') == EXPERIMENT and
                obj.get('instance_id') == INSTANCE and obj.get('region') == REGION,
                'Wrong experiment, instance, or region')
    require(complete.get('status') in ('COMPLETE', 'FAILED_SAFE_EXIT'),
            'Completion must be COMPLETE or FAILED_SAFE_EXIT')
    for field in ('fitting_terminal', 'evaluation_terminal', 'owned_processes_exited'):
        require(complete.get(field) is True, 'Completion gate not satisfied: ' + field)
    for field in ('all_processes_reviewed', 'no_other_tasks_require_instance',
                  'billing_mode_supports_stopcharging'):
        require(baseline.get(field) is True, 'Preflight gate not satisfied: ' + field)
    require(isinstance(baseline.get('reviewed_unix'), (int, float)) and
            0 < baseline['reviewed_unix'] <= time.time() + 60,
            'Invalid preflight review timestamp')
    require(isinstance(baseline.get('processes'), list) and baseline['processes'],
            'Reviewed process inventory is absent')
    required = ('pid', 'ppid', 'pgrp', 'start_ticks', 'comm', 'state')
    for row in baseline['processes']:
        require(all(k in row for k in required) and row['start_ticks'] is not None,
                'Incomplete reviewed process identity')
    identities = [(r['pid'], r['start_ticks']) for r in baseline['processes']]
    require(len(set(identities)) == len(identities), 'Duplicate process baseline identity')
    for row in baseline.get('approved_system_runtime_processes', []):
        require(row.get('role') == 'system_service' and bool(row.get('reason')) and
                (row.get('pid'), row.get('start_ticks')) in identities,
                'Unbound or unexplained system runtime exception')
    profiles = baseline.get('approved_restarting_system_services', [])
    require(isinstance(profiles, list) and len(profiles) <= 1 and
            all(system_service.valid_profile(p) for p in profiles),
            'Invalid narrow restarting-service fingerprint')
    groups = complete.get('owned_process_groups')
    require(isinstance(groups, list), 'Owned worker group identities must be recorded')
    for group in groups:
        require(all(k in group for k in ('pid', 'start_ticks', 'pgrp')) and
                group['pid'] == group['pgrp'], 'Invalid owned worker session identity')
    files = complete.get('files')
    require(isinstance(files, list) and files, 'No terminal result artifact bindings')
    roles, seen = set(), set()
    for record in files:
        p = strict_path(record.get('path', ''), OUT)
        require(str(p) not in seen, 'Duplicate lifecycle result artifact')
        seen.add(str(p))
        require(isinstance(record.get('bytes'), int) and record['bytes'] >= 0 and
                p.stat().st_size == record['bytes'], 'Terminal artifact size mismatch')
        require(re.fullmatch('[0-9a-f]{64}', record.get('sha256', '')) is not None and
                sha(p) == record['sha256'], 'Terminal artifact hash mismatch')
        require(isinstance(record.get('role'), str), 'Terminal artifact role absent')
        roles.add(record['role'])
    require('terminal_status' in roles, 'Terminal status artifact is not preserved')
    need = 'scientific_results' if complete['status'] == 'COMPLETE' else 'failure_evidence'
    require(need in roles, 'Required terminal evidence absent: ' + need)
    return complete, baseline


def process_inventory(proc_root=Path('/proc')):
    rows = []
    for entry in proc_root.iterdir():
        if not entry.name.isdigit():
            continue
        try:
            raw = (entry / 'stat').read_text()
            fields = raw[raw.rfind(')') + 2:].split()
            row = dict(pid=int(entry.name), ppid=int(fields[1]),
                             pgrp=int(fields[2]), start_ticks=int(fields[19]),
                             comm=raw[raw.find('(') + 1:raw.rfind(')')], state=fields[0])
            proof = system_service.process_proof(entry, row)
            if proof is not None:
                row['pinned_system_service_proof'] = proof
            rows.append(row)
        except (FileNotFoundError, ProcessLookupError):
            continue  # A process exited while the inventory was being read.
    return rows


def process_check(rows, baseline, complete, self_pid):
    """Ancestry exemption is exact; PID 1 descendants are NOT all exempted."""
    by_pid = {r['pid']: r for r in rows}
    require(self_pid in by_pid, 'Cannot bind finalizer process identity')
    own_lineage = {self_pid}
    p = self_pid
    while p in by_pid:
        parent = by_pid[p]['ppid']
        if parent <= 0 or parent in own_lineage:
            break
        own_lineage.add(parent)
        p = parent
    # Only descendants of this finalizer, never descendants of its ancestors.
    own_children = set()
    changed = True
    while changed:
        changed = False
        for row in rows:
            if row['ppid'] in ({self_pid} | own_children) and row['pid'] not in own_children:
                own_children.add(row['pid'])
                changed = True
    reviewed = {(r['pid'], r['start_ticks']) for r in baseline['processes']}
    services = {(r['pid'], r['start_ticks']) for r in
                baseline.get('approved_system_runtime_processes', [])}
    blockers = []
    for row in rows:
        if row['state'].startswith('Z') or row['pid'] == 2 or row['ppid'] == 2:
            continue
        if row['pid'] in own_lineage or row['pid'] in own_children:
            continue
        identity = (row['pid'], row['start_ticks'])
        reasons = []
        proven_service = system_service.matches(row, baseline.get('approved_restarting_system_services', []))
        if identity not in reviewed and not proven_service:
            reasons.append('new_or_unreviewed_process')
        if NUMERICAL.fullmatch(row['comm']):
            reasons.append('numerical_process_still_running')
        elif RUNTIME.fullmatch(row['comm']) and identity not in services:
            reasons.append('runtime_not_approved_as_system_service')
        if reasons:
            blockers.append(dict(row, reasons=reasons))
    for group in complete['owned_process_groups']:
        live = [r for r in rows if r['pgrp'] == group['pgrp'] and not r['state'].startswith('Z')]
        if live:
            blockers.append(dict(reason='owned_worker_process_group_still_live', group=group,
                                 members=live))
    return dict(observed_unix=time.time(), processes=rows, blocking_processes=blockers,
                finalizer_identity=by_pid[self_pid], own_lineage=sorted(own_lineage),
                own_children=sorted(own_children))


def safe_response(result):
    out = dict(cli_exit_code=result.returncode)
    for name in ('stdout', 'stderr'):
        raw = getattr(result, name, b'') or b''
        if isinstance(raw, str):
            raw = raw.encode()
        out[name + '_sha256'] = hashlib.sha256(raw).hexdigest()
        try:
            obj = json.loads(raw or b'{}')
            if isinstance(obj, dict):
                out[name + '_safe_fields'] = {k: obj[k] for k in ('RequestId', 'Code') if k in obj}
        except (ValueError, TypeError):
            pass
    return out


def stop_command():
    return ['aliyun', 'ecs', 'StopInstance', '--profile', PROFILE,
            '--RegionId', REGION, '--InstanceId', INSTANCE,
            '--StoppedMode', 'StopCharging', '--ForceStop', 'false']


def finalize(completion_path, completion_hash, baseline_path, baseline_hash,
             validate_only=False, inventory_fn=process_inventory, run_fn=subprocess.run,
             self_pid=None, sleep_fn=time.sleep):
    complete, baseline = validate_contract(completion_path, completion_hash,
                                          baseline_path, baseline_hash)
    pid = os.getpid() if self_pid is None else self_pid
    before = process_check(inventory_fn(), baseline, complete, pid)
    if validate_only:
        return dict(status='VALIDATED_NOT_STOPPED' if not before['blocking_processes'] else 'BLOCKED',
                    automatic_api_calls=0, process_audit=before,
                    final_stopped_verified=False, stopped_billing_verified=False)
    initial_process_audit = before
    # A bounded completion-time recheck only; no process is ignored or killed.
    # Unknown processes that remain alive still prohibit the stop request.
    if before['blocking_processes']:
        sleep_fn(2)
        before = process_check(inventory_fn(), baseline, complete, pid)
    directory = OUT / 'lifecycle'
    require(directory.resolve() == directory, 'Redirected lifecycle output directory')
    directory.mkdir(exist_ok=True)
    request_path = directory / 'STOP_REQUEST.json'
    # At-most-one actual request. An uncertain prior request requires observation,
    # not another StopInstance call.
    require(not os.path.lexists(str(request_path)),
            'STOP_REQUEST already exists; query externally and never automatically repeat')
    attempt = directory / ('attempt_' + time.strftime('%Y%m%dT%H%M%SZ', time.gmtime()) + '_' + str(pid))
    attempt.mkdir(exist_ok=False)
    audit = dict(experiment_id=EXPERIMENT, instance_id=INSTANCE, region=REGION,
                 completion_path=str(completion_path), completion_sha256=completion_hash,
                 baseline_path=str(baseline_path), baseline_sha256=baseline_hash,
                 completion_status=complete['status'], core_results_verified=True,
                 result_bindings=complete['files'], initial_process_audit=initial_process_audit,
                 bounded_recheck_process_audit=before if before is not initial_process_audit else None,
                 status='CHECKING', time_utc=now(), local_pc_or_transfer_required=False,
                 final_stopped_verified=False, stopped_billing_verified=False)
    persist(attempt / 'AUDIT.json', audit, exclusive=True)
    if before['blocking_processes']:
        audit.update(status='BLOCKED_OTHER_OR_UNRESOLVED_PROCESSES', actual_stop_requests=0)
        persist(attempt / 'RESULT.json', audit, exclusive=True)
        return audit
    cmd = stop_command()
    try:
        dry = run_fn(cmd + ['--DryRun', 'true'], stdout=subprocess.PIPE,
                     stderr=subprocess.PIPE, timeout=25)
        summary = safe_response(dry)
        codes = [summary.get(k + '_safe_fields', {}).get('Code') for k in ('stdout', 'stderr')]
        dry_ok = 'DryRunOperation' in codes
        # Older aliyun CLI may print the structured error code in a text envelope.
        if not dry_ok:
            dry_ok = b'DryRunOperation' in ((dry.stdout or b'') + (dry.stderr or b''))
        audit['dry_run'] = summary
        require(dry_ok, 'Normal StopCharging dry-run did not confirm permission')
        # Revalidate all critical results and process state immediately before API.
        validate_contract(completion_path, completion_hash, baseline_path, baseline_hash)
        latest = process_check(inventory_fn(), baseline, complete, pid)
        audit['final_process_audit'] = latest
        require(not latest['blocking_processes'], 'A new or unresolved process blocks normal stop')
    except BaseException as error:
        audit.update(status='BLOCKED_PREFLIGHT', error_type=type(error).__name__,
                     reason=str(error) if isinstance(error, RuntimeError) else 'Dry-run or verification failed',
                     actual_stop_requests=0)
        persist(attempt / 'RESULT.json', audit, exclusive=True)
        return audit
    proof = dict(schema='multifsync_stop_request_v1', experiment_id=EXPERIMENT,
                 instance_id=INSTANCE, region=REGION, requested_mode='StopCharging',
                 force_stop=False, request_state='PENDING', request_accepted=None,
                 automatic_retry_allowed=False, actual_stop_requests=1,
                 completion_sha256=completion_hash, baseline_sha256=baseline_hash,
                 audit_path=str(attempt / 'AUDIT.json'), audit_sha256=sha(attempt / 'AUDIT.json'),
                 process_audit=latest, dry_run=summary, requested_utc=now(),
                 external_state_confirmation_pending=True, final_stopped_verified=False,
                 stopped_billing_verified=False)
    persist(request_path, proof, exclusive=True)
    try:
        result = run_fn(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=25)
        response = safe_response(result)
        accepted = result.returncode == 0 and bool(response.get('stdout_safe_fields', {}).get('RequestId'))
        proof.update(response=response, request_state='ACCEPTED' if accepted else 'UNKNOWN',
                     request_accepted=True if accepted else None)
    except BaseException as error:
        proof.update(request_state='UNKNOWN', request_accepted=None,
                     error_type=type(error).__name__)
    finally:
        proof['last_recorded_utc'] = now()
        persist(request_path, proof)
        audit.update(status='STOP_REQUEST_' + proof['request_state'],
                     request_path=str(request_path), request_sha256=sha(request_path),
                     request_accepted=proof['request_accepted'], actual_stop_requests=1,
                     external_state_confirmation_pending=True)
        persist(attempt / 'RESULT.json', audit, exclusive=True)
    return audit


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--completion', required=True)
    parser.add_argument('--completion-sha256', required=True)
    parser.add_argument('--baseline', required=True)
    parser.add_argument('--baseline-sha256', required=True)
    parser.add_argument('--validate-only', action='store_true')
    args = parser.parse_args()
    require(sys.platform.startswith('linux'), 'ECS finalizer is Linux only')
    result = finalize(args.completion, args.completion_sha256,
                      args.baseline, args.baseline_sha256, args.validate_only)
    brief = {k: result[k] for k in ('status', 'request_path', 'request_sha256', 'request_accepted',
             'actual_stop_requests', 'external_state_confirmation_pending',
             'final_stopped_verified', 'stopped_billing_verified') if k in result}
    print(json.dumps(brief), flush=True)
    return 0 if result['status'] in ('VALIDATED_NOT_STOPPED', 'STOP_REQUEST_ACCEPTED') else 2


if __name__ == '__main__':
    try:
        sys.exit(main())
    except Exception as error:
        # Safe errors contain no raw CLI response, secret, environment, or args.
        print(json.dumps(dict(status='FINALIZER_BLOCKED', error_type=type(error).__name__,
              reason=str(error) if isinstance(error, RuntimeError) else 'Inspect persisted lifecycle evidence',
              final_stopped_verified=False, stopped_billing_verified=False)), flush=True)
        sys.exit(2)
