#!/usr/bin/env python3
"""One Linux-only preparation attempt; no registered fitting or actual ECS stop."""
import csv
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

EX = Path('/root/cavi-true222-fit444-m4-exec-20261117-v1')
OUT = Path('/root/multifsync-results/cavi_true222_fit444_m4_20261117_v1')
STAGE = Path('/root/cavi-true222-fit444-m4-stage-20261117-v1')
THREADS = ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
           'BLIS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'NUMEXPR_NUM_THREADS')


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(1048576), b''):
            digest.update(chunk)
    return digest.hexdigest()


def read(path):
    return json.loads(Path(path).read_text(encoding='utf-8'))


def save(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('x', encoding='utf-8') as stream:
        json.dump(value, stream, indent=2)
        stream.write('\n')
        stream.flush()
        os.fsync(stream.fileno())


def verify_manifest(path, root, exact_names=None):
    names = []
    for line in Path(path).read_text().splitlines():
        require(len(line) > 66 and line[64:66] == '  ', 'Malformed manifest')
        rel = Path(line[66:])
        require(not rel.is_absolute() and '..' not in rel.parts, 'Unsafe manifest path')
        item = root / rel
        require(item.resolve() == item and sha(item) == line[:64], 'Manifest mismatch: ' + str(rel))
        names.append(rel.as_posix())
    require(len(names) > 0 and len(names) == len(set(names)), 'Empty or duplicate manifest')
    if exact_names is not None:
        require(set(names) == set(exact_names), 'Unexpected projected input whitelist')
    return names


def safe_cli(result):
    answer = {'cli_exit_code': result.returncode}
    for label in ('stdout', 'stderr'):
        value = getattr(result, label, '') or ''
        answer[label + '_sha256'] = hashlib.sha256(value.encode('utf-8')).hexdigest()
        try:
            obj = json.loads(value or '{}')
            if isinstance(obj, dict):
                answer[label + '_safe_fields'] = {key: obj[key] for key in ('Code', 'RequestId') if key in obj}
        except ValueError:
            pass
    return answer


def main():
    require(sys.platform.startswith('linux'), 'Linux execution only')
    require(all(path.resolve() == path for path in (EX, OUT, STAGE)), 'Redirected roots')
    require(STAGE.is_dir() and EX.is_dir(), 'Upload stage and execution package must exist')
    logs = STAGE / 'prepare_logs'
    require(not logs.exists(), 'Preparation was already attempted; no automatic retry')
    logs.mkdir()
    save(logs / 'START.json', {'started_unix': time.time(), 'pid': os.getpid(),
                              'registered_fit_attempts': 0, 'automatic_retries': 0})
    results = []
    prepared = None
    try:
        config = read(EX / 'config/EXECUTION_CONFIG.json')
        require(config['status'] == 'FROZEN' and config['new_fit_count'] == 120, 'Configuration not frozen')
        require(config['arms'] == ['original'] and config['fitted_factor_capacities'] == [4, 4, 4]
                and config['fitted_fpca_capacity'] == 4, 'Unexpected model configuration')
        require(config['planned_new_artifact_bytes'] == 6 * 1024**3
                and config['minimum_free_disk_bytes'] == 3 * 1024**3, 'Unexpected frozen storage budget')
        frozen = read(EX / 'config/PROTOCOL_FROZEN.json')
        require(frozen['status'] == 'FROZEN' and sha(EX / 'config/PROTOCOL.json') == frozen['protocol_sha256'],
                'Protocol binding differs')
        verify_manifest(EX / 'CODE_SHA256SUMS.txt', EX)
        require(shutil.disk_usage(str(EX)).free >= config['required_initial_free_disk_bytes'],
                'Initial free disk is below the registered artifact allowance plus reserve')
        dep_source = EX / 'RUNTIME_LIBRARY_BINDING.json'
        dep = read(dep_source)
        for item in dep['package_files']:
            path = Path(dep['package_root']) / item['relative_path']
            require(path.is_file() and sha(path) == item['sha256'], 'Bound runtime dependency changed')
        dep_dest = STAGE / 'RUNTIME_LIBRARY_BINDING.json'
        if dep_dest.exists():
            require(sha(dep_dest) == sha(dep_source), 'Stage runtime binding differs')
        else:
            with dep_dest.open('xb') as stream:
                stream.write(dep_source.read_bytes())
        require(sha(dep_dest) == sha(dep_source), 'Runtime binding copy failed its hash check')
        env = dict(os.environ, PYTHONDONTWRITEBYTECODE='1', FSYNC_EXEC_ROOT=str(EX),
                   FSYNC_OUTPUT_ROOT=str(OUT), R_LIBS=dep['R_LIBS'], R_LIBS_USER=dep['R_LIBS_USER'])
        for key in THREADS:
            env[key] = '1'

        def run(name, command):
            verify_manifest(EX / 'CODE_SHA256SUMS.txt', EX)
            began = time.monotonic()
            save(logs / (name + '.START.json'), {'command': command, 'started_unix': time.time()})
            with (logs / (name + '.log')).open('x', encoding='utf-8') as stream:
                child = subprocess.Popen(command, cwd=str(EX), env=env, stdin=subprocess.DEVNULL,
                                         stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                         universal_newlines=True, bufsize=1)
                try:
                    for line in child.stdout:
                        stream.write(line)
                        stream.flush()
                    code = child.wait()
                finally:
                    if child.poll() is None:
                        child.wait()
                    stream.flush()
                    os.fsync(stream.fileno())
            result = {'stage': name, 'exit_code': code, 'wall_seconds': time.monotonic() - began}
            save(logs / (name + '.EXIT.json'), result)
            results.append(result)
            require(code == 0, 'Preparation stage failed: ' + name + '; no automatic retry')
            verify_manifest(EX / 'CODE_SHA256SUMS.txt', EX)

        # Only permission checking: --DryRun true can never stop the instance.
        dry_command = ['aliyun', 'ecs', 'StopInstance', '--profile', 'multifsync-self-stop',
                       '--RegionId', 'cn-hangzhou', '--InstanceId', 'i-bp1egmsairexm9w1rnks',
                       '--StoppedMode', 'StopCharging', '--ForceStop', 'false', '--DryRun', 'true']
        dry = subprocess.run(dry_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                             universal_newlines=True, timeout=30, env=env)
        dry_pass = 'DryRunOperation' in ((dry.stdout or '') + (dry.stderr or ''))
        dry_receipt = safe_cli(dry)
        dry_receipt.update({'status': 'PASS' if dry_pass else 'FAIL', 'dry_run': True,
                            'actual_stop_requests': 0, 'final_stopped_verified': False,
                            'stopped_billing_verified': False})
        save(STAGE / 'STOPCHARGING_DRYRUN.json', dry_receipt)
        require(dry_pass, 'Normal StopCharging permission dry-run did not pass')

        run('generation', ['Rscript', '--vanilla', str(EX / 'generation/generate_registered.R'),
                           str(EX), str(OUT), config['baseline_library']])
        generation = read(OUT / 'GENERATION_STATUS.json')
        require(generation['status'] == 'PASS' and generation['new_datasets'] == 10
                and generation['reused_datasets'] == 0 and generation['registered_fits'] == 120,
                'Generation completion differs')
        require(sha(EX / 'config/PLAN.csv') == sha(OUT / 'PLAN.csv'), 'Copied plan differs')
        with (OUT / 'PLAN.csv').open(encoding='utf-8-sig', newline='') as stream:
            plan = list(csv.DictReader(stream))
        ids = sorted(set(row['dataset_id'] for row in plan))
        require(ids == ['NC01_TRUE222_20261117__rep%02d' % index for index in range(1, 11)]
                and len(plan) == 120 and len(set(row['fit_id'] for row in plan)) == 120,
                'Dataset/fit registry differs')
        for data_id in ids:
            subset = [row for row in plan if row['dataset_id'] == data_id]
            require(sorted(int(row['start']) for row in subset) == list(range(1, 13))
                    and all(row['arm'] == 'original' for row in subset), 'Start registry differs')
        projected = 0
        for data_id in ids:
            folder = OUT / 'fit_inputs' / data_id
            projected += len(verify_manifest(folder / 'INPUT_SHA256SUMS.txt', folder,
                                              ('observation_bundle.rds', 'model_metadata.rds')))
        require(projected == 20 and not (OUT / 'fits').exists(), 'Unexpected projected inputs or formal fit attempt')
        prepared = {'status': 'PASS', 'new_datasets': 10, 'registered_fits': 120,
                    'projected_input_files_verified': projected, 'formal_fits': 0,
                    'generation_status_sha256': sha(OUT / 'GENERATION_STATUS.json'),
                    'source_manifest_sha256': sha(EX / 'CODE_SHA256SUMS.txt'),
                    'registration_sha256': sha(EX / 'config/DATASET_REGISTRATION.csv')}
        run('startup', ['Rscript', '--vanilla', str(EX / 'scripts/startup_acceptance.R')])
        require(read(OUT / 'startup/RUNNER_ACCEPTANCE.json')['status'] == 'PASS', 'Original startup gate failed')
        oracle_output = OUT / 'science_tests/oracle'
        summary_output = OUT / 'science_tests/summary'
        run('science_oracle', ['Rscript', '--vanilla', str(EX / 'tests/science_oracle.R'), str(EX), str(oracle_output)])
        run('science_summary', ['python3', '-B', str(EX / 'tests/science_summary.py'), str(EX), str(summary_output)])
        for folder, name in ((oracle_output, 'science_ORACLE_PASS.json'), (summary_output, 'science_SUMMARY_PASS.json')):
            require(read(folder / name)['status'] == 'PASS', 'Science test failed: ' + name)
        run('resource_probe', ['Rscript', '--vanilla', str(EX / 'scripts/resource_probe.R')])
        run('science_probe', ['Rscript', '--vanilla', str(EX / 'science/evaluate_timing_endpoint.R'),
                             '--exec-root', str(EX), '--output-root', str(OUT), '--data-id', ids[0],
                             '--fit-path', str(OUT / 'resource_probe/original_L4M4_2sweep/PROBE_FIT.rds')])
        require(read(OUT / 'resource_probe/SCIENCE_SMOKE/PASS.json')['status'] == 'PASS', 'Actual science smoke failed')
        run('resource_assessment', ['python3', '-B', str(EX / 'control/assess_resource_probe.py')])
        require(read(STAGE / 'RESOURCE_PROBE_ACCEPTANCE.json')['status'] == 'PASS', 'Resource assessment failed')
        require(not (OUT / 'fits').exists(), 'Preparation unexpectedly created formal fit candidates')
        save(STAGE / 'PREPARE_ACCEPTANCE.json', prepared)
        runtime = {'status': 'PASS', 'stopcharging_permission_dryrun': True,
                   'startup_original_gate': True, 'science_tests': True, 'registration_checked': True,
                   'runtime_library_binding_sha256': sha(dep_dest),
                   'runtime_library_binding_execution_sha256': sha(dep_source),
                   'runtime_library_binding_stage_execution_identical': sha(dep_dest) == sha(dep_source),
                   'actual_science_smoke': True, 'resource_probe_acceptance': True,
                   'formal_fit_attempts': 0, 'actual_stop_requests': 0, 'stages': results}
        save(STAGE / 'RUNTIME_ACCEPTANCE.json', runtime)
        print(json.dumps({'status': 'PASS', 'new_datasets': 10, 'formal_fit_attempts': 0,
                          'actual_stop_requests': 0, 'launcher_invoked': False}))
        return 0
    except BaseException as error:
        failure = {'status': 'FAIL', 'error_type': type(error).__name__, 'error': str(error),
                   'stages': results, 'automatic_retries': 0, 'actual_stop_requests': 0,
                   'formal_fit_attempts': 0}
        save(logs / 'FAILURE.json', failure)
        if not (STAGE / 'RUNTIME_ACCEPTANCE.json').exists():
            save(STAGE / 'RUNTIME_ACCEPTANCE.json', failure)
        if not (STAGE / 'PREPARE_ACCEPTANCE.json').exists():
            save(STAGE / 'PREPARE_ACCEPTANCE.json', prepared or failure)
        print(json.dumps({'status': 'FAIL', 'error_type': type(error).__name__, 'error': str(error)}))
        return 2


if __name__ == '__main__':
    sys.exit(main())
