#!/usr/bin/env python3
"""Thread-compatible observation report; no independent process or scheduling changes."""
from pathlib import Path
import argparse, csv, ctypes, hashlib, json, math, os, re, select, shutil, sys, time, threading

EXPERIMENT = 'cavi_true222_fit444_m4_20261117_v1'
GiB, MiB = 1024 ** 3, 1024 ** 2
HEX = re.compile(r'^[0-9a-f]{64}$')


class AuditError(RuntimeError):
    pass


class Terminal(RuntimeError):
    pass


def require(ok, message):
    if not ok:
        raise AuditError(message)


def strict(path, root):
    path, root = Path(path), Path(root)
    require(path.is_absolute() and path.resolve() == path and root.resolve() == root,
            'Redirected or nonabsolute observation path')
    try:
        path.relative_to(root)
    except ValueError:
        raise AuditError('Observation path escaped its registered root')
    return path


def sha(path, check=lambda: None):
    before = Path(path).stat()
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        while True:
            check()
            chunk = stream.read(MiB)
            if not chunk:
                break
            digest.update(chunk)
    after = Path(path).stat()
    require((before.st_dev, before.st_ino, before.st_size, before.st_mtime_ns) ==
            (after.st_dev, after.st_ino, after.st_size, after.st_mtime_ns), 'File changed during hashing')
    return digest.hexdigest()


def read(path):
    return json.loads(Path(path).read_text(encoding='utf-8-sig'))


def identity(pid, proc_root=Path('/proc')):
    try:
        fields = (Path(proc_root) / str(pid) / 'stat').read_text().rsplit(')', 1)[1].split()
        return dict(pid=int(pid), state=fields[0], pgrp=int(fields[2]), start_ticks=int(fields[19]))
    except (FileNotFoundError, ProcessLookupError):
        return None


def terminal_reason(out, launch, identity_fn=identity):
    if (Path(out) / 'FINAL_STATUS.json').exists():
        return 'FINAL_STATUS_EXISTS'
    current = identity_fn(launch['pid'])
    if not current or current['state'].startswith('Z'):
        return 'REGISTERED_PIPELINE_EXITED'
    if any(current[k] != launch[k] for k in ('pid', 'pgrp', 'start_ticks')):
        return 'REGISTERED_PIPELINE_IDENTITY_CHANGED'
    return None


def verify_frozen(ex, expected, check=lambda: None):
    manifest = strict(Path(ex) / 'CODE_SHA256SUMS.txt', ex)
    require(HEX.fullmatch(expected or '') is not None and sha(manifest, check) == expected,
            'Frozen source manifest hash mismatch')
    lines = manifest.read_text(encoding='utf-8').splitlines()
    require(len(lines) > 35, 'Incomplete frozen source manifest')
    for line in lines:
        check()
        require(re.match(r'^[0-9a-f]{64}  .+', line) is not None, 'Malformed frozen manifest')
        rel = Path(line[66:])
        require(not rel.is_absolute() and '..' not in rel.parts, 'Unsafe manifest relative path')
        require(sha(strict(Path(ex) / rel, ex), check) == line[:64], 'Frozen source changed: ' + str(rel))
    return dict(status='PASS', manifest_sha256=expected, files_checked=len(lines))


def read_plan(ex, out, check):
    source, copied = Path(ex) / 'config/PLAN.csv', Path(out) / 'PLAN.csv'
    require(sha(source, check) == sha(copied, check), 'Output PLAN differs from frozen PLAN')
    with copied.open(encoding='utf-8-sig', newline='') as stream:
        rows = list(csv.DictReader(stream))
    require(len(rows) == 120 and len({r['fit_id'] for r in rows}) == 120 and
            len({r['dataset_id'] for r in rows}) == 10 and all(r['arm'] == 'original' for r in rows),
            'Wrong 10-dataset / 120-fit single-arm registry')
    require(all(re.fullmatch(r'[A-Za-z0-9_.-]+', r['fit_id']) for r in rows), 'Unsafe fit ID')
    return rows


def tree_size(folder, root, check=lambda: None):
    """Size closed files only; partial/temp files are listed but excluded from model size."""
    records, excluded = [], []
    for path in sorted(Path(folder).rglob('*')):
        check(); strict(path, root)
        if not path.is_file():
            continue
        rel = path.relative_to(folder).as_posix()
        record = dict(relative_path=rel, bytes=path.stat().st_size)
        parts = path.relative_to(folder).parts
        temporary = path.name.endswith(('.part', '.partial', '.tmp', '~')) or any(
            x.lower() in ('tmp', 'temp', 'partial') for x in parts)
        (excluded if temporary else records).append(record)
    return dict(bytes=sum(r['bytes'] for r in records), files=records, excluded_partial_files=excluded,
                excluded_partial_bytes=sum(r['bytes'] for r in excluded))


def verify_storage(folder, receipt_name, artifact_name, root, check):
    receipt_path = strict(folder / receipt_name, root)
    receipt = read(receipt_path)
    artifact = strict(folder / artifact_name, root)
    require(Path(receipt.get('path', '')) == artifact, 'Storage receipt artifact path differs')
    require(receipt.get('identical_roundtrip') is True and receipt.get('compression') == 'xz' and
            receipt.get('serialization_version') == 3 and HEX.fullmatch(receipt.get('object_sha256', '')),
            'Producer storage roundtrip attestation missing')
    digest = sha(artifact, check)
    require(digest == receipt.get('sha256') and artifact.stat().st_size == receipt.get('bytes'),
            'Full artifact hash or size differs from roundtrip receipt')
    return dict(path=str(artifact), bytes=artifact.stat().st_size, sha256=digest,
                receipt_sha256=sha(receipt_path, check), producer_identical_roundtrip=True,
                producer_object_sha256=receipt['object_sha256'], independent_R_object_recheck=False)


def endpoint(out, fit_id, expected_manifest, check):
    folder = strict(Path(out) / 'fits' / fit_id, out)
    worker = strict(Path(out) / 'scheduler' / fit_id, out)
    exit_path, status_path = worker / 'EXIT.json', folder / 'STATUS.csv'
    if not exit_path.is_file() or not status_path.is_file():
        return None
    ended = read(exit_path)
    if ended.get('status') != 'SUCCESS':
        return None
    with status_path.open(encoding='utf-8-sig', newline='') as stream:
        status = list(csv.DictReader(stream))
    require(len(status) == 1 and status[0].get('fit_id') == fit_id and
            status[0].get('status') == 'complete', 'SUCCESS worker lacks one complete fit STATUS')
    status = status[0]
    require(int(status['t1_sweeps']) == 400 and int(status['annealing_sweeps']) == 99,
            'Successful endpoint did not complete frozen update budget')
    require(ended.get('frozen_bindings_pass') is True and ended.get('child_still_running') is False,
            'Successful worker did not attest closed child and frozen bindings')
    binding = read(strict(folder / 'SOURCE_BINDING.json', out))
    require(binding.get('fit_id') == fit_id and binding.get('source_manifest_sha256') == expected_manifest,
            'Endpoint source binding differs from frozen source')
    full = verify_storage(folder, 'STORAGE_ROUNDTRIP.json', 'fit.rds', out, check)
    compact = verify_storage(folder, 'COMPACT_STORAGE_ROUNDTRIP.json', 'COMPACT_DIAGNOSTICS.rds', out, check)
    require(status['fit_sha256'] == full['sha256'] and int(float(status['fit_bytes'])) == full['bytes'],
            'STATUS fit hash/size differs from full storage receipt')
    sizes = tree_size(folder, out, check)
    logs = [dict(name=name, bytes=(worker / name).stat().st_size) for name in
            ('stdout.log', 'stderr.log', 'worker_stdout.log', 'worker_stderr.log') if (worker / name).is_file()]
    peak = ended.get('peak_rss_bytes')
    require(isinstance(peak, int) and peak >= 0, 'Worker peak RSS missing')
    return dict(fit_id=fit_id, finished_unix=ended.get('finished_unix'), full_storage=full,
                compact_storage=compact, directory=sizes, worker_logs=logs,
                candidate_and_log_bytes=sizes['bytes'] + sum(r['bytes'] for r in logs),
                worker_peak_rss_bytes=peak, worker_resource_sample_count=ended.get('sample_count'),
                worker_elapsed_seconds=ended.get('elapsed_seconds'),
                worker_exit_sha256=sha(exit_path, check), status_csv_sha256=sha(status_path, check),
                resource_sampling_limitation='Two-second RSS samples can miss peaks; no full-runtime memory bound.')


def forecast(config, probe, endpoints, free_bytes):
    completed = len(endpoints)
    largest = max(row['candidate_and_log_bytes'] for row in endpoints)
    allowance = max(int(probe['per_candidate_pair_allowance_bytes']), int(math.ceil(1.25 * largest)))
    unfinished = config['new_fit_count'] - completed
    science = max(512 * MiB, int(probe.get('science_and_controller_allowance_bytes', 0)))
    # Extra slot-sized allowance for every possible concurrent atomic write.
    temporary = config['max_concurrency'] * allowance
    reserve = config['minimum_free_disk_bytes']
    required = unfinished * allowance + science + temporary + reserve
    return dict(status='WITHIN_PADDED_FORECAST' if free_bytes >= required else 'INSUFFICIENT_PADDED_FORECAST',
                completed_verified_count=completed, remaining_conservative_count=unfinished,
                per_remaining_candidate_allowance_bytes=allowance, largest_complete_candidate_bytes=largest,
                remaining_candidate_bytes=unfinished * allowance, science_controller_allowance_bytes=science,
                concurrent_temporary_write_allowance_bytes=temporary, disk_reserve_bytes=reserve,
                current_free_disk_bytes=free_bytes, required_remaining_free_bytes=required,
                headroom_after_forecast_bytes=free_bytes - required,
                partial_files_excluded_from_candidate_model_bytes=sum(e['directory']['excluded_partial_bytes'] for e in endpoints),
                existing_and_partial_files_already_reflected_in_current_free_space=True,
                failed_or_unattempted_fits_counted_conservatively_without_authorizing_retry=True,
                full_storage_upper_bound_proven=False,
                limitation='Padding is an estimate from completed endpoints, not a bound on later file sizes or unrelated disk use.')


class WakeOnOutput:
    """Linux completion-event wakeup; no persistent observer without this facility."""
    def __init__(self, out):
        self.lib = ctypes.CDLL(None, use_errno=True)
        self.lib.inotify_init1.argtypes = [ctypes.c_int]
        self.lib.inotify_add_watch.argtypes = [ctypes.c_int, ctypes.c_char_p, ctypes.c_uint32]
        self.fd = self.lib.inotify_init1(os.O_NONBLOCK | os.O_CLOEXEC)
        require(self.fd >= 0, 'inotify unavailable; refuse long-lived observer')
        mask = 0x100 | 0x80 | 0x8 | 0x400 | 0x800  # CREATE, MOVED_TO, CLOSE_WRITE, DELETE_SELF, MOVE_SELF
        if self.lib.inotify_add_watch(self.fd, os.fsencode(str(out)), mask) < 0:
            os.close(self.fd)
            raise AuditError('Cannot watch completion marker; observer will not wait')

    def wait(self, seconds, check):
        deadline = time.monotonic() + seconds
        while time.monotonic() < deadline:
            check()
            ready, _, _ = select.select([self.fd], [], [], max(0, deadline - time.monotonic()))
            if ready:
                try:
                    os.read(self.fd, 65536)
                except BlockingIOError:
                    pass
            check()

    def close(self):
        os.close(self.fd)


def write_report(path, report):
    path.parent.mkdir(exist_ok=True)
    require(path.parent.resolve() == path.parent and not path.exists(), 'Existing or redirected observer report')
    with path.open('x', encoding='utf-8') as stream:
        json.dump(report, stream, indent=2); stream.write('\n'); stream.flush(); os.fsync(stream.fileno())


def await_started_launch(path, max_wait=5.0, clock=time.monotonic, sleeper=time.sleep):
    """Allow the launcher a bounded interval to persist the pipeline identity."""
    deadline = clock() + max_wait
    while True:
        try:
            value = read(path)
        except (FileNotFoundError, json.JSONDecodeError):
            value = None
        if isinstance(value, dict) and value.get('status') == 'STARTED':
            return value
        left = deadline - clock()
        if left <= 0:
            raise AuditError('PIPELINE_LAUNCH did not become STARTED within five seconds')
        sleeper(min(0.05, left))


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source-manifest-sha256', required=True)
    parser.add_argument('--observer-source-sha256', required=True)
    args = parser.parse_args(argv)
    require(sys.platform.startswith('linux'), 'Observer execution is Linux only')
    require(threading.current_thread() is not threading.main_thread(),
            'Thread observer must run within the registered pipeline thread; standalone execution is forbidden')
    ex = Path('/root/cavi-true222-fit444-m4-exec-20261117-v1')
    out = Path('/root/multifsync-results/' + EXPERIMENT)
    stage = Path('/root/cavi-true222-fit444-m4-stage-20261117-v1')
    destination = out / 'resource_monitor/FIRST_ENDPOINT_STORAGE_REASSESSMENT.json'
    for p in (ex, out, stage):
        strict(p, p)
    require(not destination.exists(), 'Observer report already exists; no repeat')
    own_sha = sha(Path(__file__).resolve())
    require(own_sha == args.observer_source_sha256, 'Observer source hash differs')
    report = dict(schema='first_endpoint_resource_observation_v1', experiment_id=EXPERIMENT,
                  observer_pid=os.getpid(), observer_source_sha256=own_sha,
                  expected_frozen_source_manifest_sha256=args.source_manifest_sha256,
                  observed_unix=time.time(), cloud_api_calls=0, process_signals_sent=0,
                  scheduler_state_modified=False, fitting_modified=False, runtime_poll_seconds=30,
                  completion_event_wakeup=True, hard_resource_guarantee=False,
                  observer_mode='THREAD_IN_REGISTERED_PIPELINE', standalone_process_forbidden=True)
    wake = None
    try:
        launch_path = strict(stage / 'PIPELINE_LAUNCH.json', stage)
        launch = await_started_launch(launch_path)
        require(launch.get('status') == 'STARTED' and launch.get('experiment_id') == EXPERIMENT and
                launch.get('registered_fits') == 120 and launch.get('registered_datasets') == 10 and
                launch.get('pid') == launch.get('pgrp') and isinstance(launch.get('start_ticks'), int),
                'No valid registered pipeline identity')
        report['pipeline_launch_sha256'] = sha(launch_path)
        report['pipeline_identity'] = {key: launch[key] for key in ('pid', 'pgrp', 'start_ticks')}
        def check():
            why = terminal_reason(out, launch)
            if why:
                raise Terminal(why)
        check()
        wake = WakeOnOutput(out)
        report['initial_frozen_source_verification'] = verify_frozen(ex, args.source_manifest_sha256, check)
        plan = read_plan(ex, out, check)
        config_path = stage / 'amendments/ops_v3/EXECUTION_CONFIG.json'
        probe_path = stage / 'RESOURCE_PROBE_ACCEPTANCE_V3.json'
        config, probe = read(config_path), read(probe_path)
        require(config.get('status') == 'FROZEN' and config.get('new_fit_count') == 120 and
                config.get('max_concurrency') == 8 and config.get('minimum_free_disk_bytes') == 3 * GiB,
                'Execution resource contract changed')
        require(config.get('planned_new_artifact_bytes') == 27 * GiB // 4,
                'OPS_V3 must use the separate 6.75 GiB artifact allowance')
        require(probe.get('status') == 'PASS' and probe.get('registered_fit_count') == 120,
                'Registered resource probe acceptance absent')
        report['resource_probe_acceptance_sha256'] = sha(probe_path, check)
        report['ops_v3_execution_config_sha256'] = sha(config_path, check)
        while True:
            check()
            complete = []
            for row in plan:
                check()
                value = endpoint(out, row['fit_id'], args.source_manifest_sha256, check)
                if value is not None:
                    complete.append(value)
            if complete:
                complete.sort(key=lambda x: (x['finished_unix'] if isinstance(x['finished_unix'], (int, float)) else float('inf'), x['fit_id']))
                report['first_verified_endpoint'] = complete[0]
                report['verified_completed_endpoints'] = complete
                report['forecast'] = forecast(config, probe, complete, shutil.disk_usage(str(out)).free)
                report['final_frozen_source_verification'] = verify_frozen(ex, args.source_manifest_sha256, check)
                require(sha(config_path, check) == report['ops_v3_execution_config_sha256'] and
                        sha(probe_path, check) == report['resource_probe_acceptance_sha256'],
                        'OPS_V3 resource inputs changed during observation')
                check()
                report['status'] = 'REASSESSED'
                report['decision_authority'] = 'Observation only; no automatic admission, cancellation, cleanup or shutdown.'
                break
            wake.wait(30, check)
    except Terminal as error:
        for key in ('first_verified_endpoint', 'verified_completed_endpoints', 'forecast'):
            report.pop(key, None)
        report.update(status='TERMINATED_WITHOUT_COMPLETED_REASSESSMENT', termination_reason=str(error),
                      final_source_verification_not_claimed=True,
                      no_complete_endpoint_claimed=True)
    except Exception as error:
        report.update(status='OBSERVATION_FAILED', error_type=type(error).__name__,
                      reason=str(error) if isinstance(error, AuditError) else 'Read-only observer failed; no state was changed.')
    finally:
        if wake is not None:
            wake.close()
    report['finished_unix'] = time.time()
    write_report(destination, report)
    print(json.dumps(dict(status=report['status'], report=str(destination), shutdown_requested=False)))
    return 0 if report['status'] in ('REASSESSED', 'TERMINATED_WITHOUT_COMPLETED_REASSESSMENT') else 2


if __name__ == '__main__':
    raise SystemExit(main())
