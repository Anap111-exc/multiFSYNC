#!/usr/bin/env python3
"""Conservative launch screen from a two-update probe, not a full-fit bound."""
import csv
import hashlib
import json
import math
from pathlib import Path
import shutil
import sys
import time

EX = Path('/root/cavi-true222-fit444-m4-exec-20261117-v1')
OUT = Path('/root/multifsync-results/cavi_true222_fit444_m4_20261117_v1')
STAGE = Path('/root/cavi-true222-fit444-m4-stage-20261117-v1')
GiB = 1024**3
MiB = 1024**2


def read(path):
    return json.loads(Path(path).read_text(encoding='utf-8'))


def sha(path):
    result = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(MiB), b''):
            result.update(chunk)
    return result.hexdigest()


def tree_bytes(root):
    total = 0
    for path in root.rglob('*'):
        if path.is_symlink():
            raise RuntimeError('Redirected artifact path')
        if path.is_file():
            total += path.stat().st_size
    return total


def estimate(config, probe, fields, existing_output_bytes, package_bytes, free_bytes, available_bytes):
    """Pure estimate. Explicit padding admits unmeasured full endpoint growth."""
    # Formal endpoints contain two checkpoints; probe contains one. Snapshot history
    # is internal memory, not 61 serialized copies. Diagnostics grow to 400 rows.
    checkpoint = fields.get('practical_checkpoints', 0)
    compact_now = sum(fields.get(key, 0) for key in
                      ('initialization_diagnostics', 'driver_trace', 'practical_final', 'practical_checkpoints'))
    small_diagnostic_padding = 4 * MiB
    uncompressed_full_pair = probe['fit_object_bytes'] + compact_now + 2 * checkpoint + small_diagnostic_padding
    measured_pair = probe['fit_xz_bytes'] + probe['compact_xz_bytes']
    # An allowance, not a theorem about 499-step compression or posterior growth.
    per_candidate = int(math.ceil(max(uncompressed_full_pair * 1.25, measured_pair * 3.0)))
    fitting_artifacts = config['new_fit_count'] * per_candidate
    science_and_controller_allowance = 512 * MiB
    prediction = fitting_artifacts + science_and_controller_allowance + existing_output_bytes + package_bytes
    remaining_prediction = fitting_artifacts + science_and_controller_allowance
    snapshot_history = int(probe.get('full_history_fitted_vectors_uncompressed_estimate_bytes', 0))
    per_process = max(GiB, int(math.ceil(probe['peak_rss_mib'] * MiB * 2)) + snapshot_history)
    memory_reserve = 3 * GiB
    max_supported = max(0, min(config['max_concurrency'], (available_bytes - memory_reserve) // per_process))
    required_disk = max(config['planned_new_artifact_bytes'], remaining_prediction) + config['minimum_free_disk_bytes']
    failures = []
    if prediction > config['planned_new_artifact_bytes']:
        failures.append('Conservative artifact estimate exceeds the frozen 6 GiB allowance')
    if free_bytes < required_disk:
        failures.append('Free disk cannot support the frozen allowance and 3 GiB reserve')
    if max_supported < config['initial_concurrency']:
        failures.append('Available memory does not support the first two formal workers plus reserve')
    return {'status': 'PASS' if not failures else 'BLOCKED', 'blockers': failures,
            'predicted_new_artifact_bytes': prediction, 'predicted_remaining_artifact_bytes': remaining_prediction,
            'per_candidate_pair_allowance_bytes': per_candidate, 'candidate_artifact_bytes': fitting_artifacts,
            'science_and_controller_allowance_bytes': science_and_controller_allowance,
            'existing_output_bytes': existing_output_bytes, 'uploaded_package_bytes': package_bytes,
            'free_disk_bytes': free_bytes, 'required_free_disk_bytes': required_disk,
            'disk_reserve_bytes': config['minimum_free_disk_bytes'], 'artifact_budget_bytes': config['planned_new_artifact_bytes'],
            'measured_pair_xz_bytes': measured_pair, 'uncompressed_full_pair_proxy_bytes': uncompressed_full_pair,
            'conservative_process_memory_bytes': per_process, 'mem_available_bytes': available_bytes,
            'system_memory_reserve_bytes': memory_reserve, 'initial_concurrency': config['initial_concurrency'],
            'max_concurrency_from_short_probe': int(max_supported),
            'scheduler_must_measure_active_rss_before_growth': True,
            'full_state_size_upper_bound_proven': False, 'full_fit_time_upper_bound_proven': False,
            'first_formal_endpoint_storage_reassessment_required': True,
            'limitation': 'Two ordinary updates do not bound 499-step runtime, memory or compression. These padded allowances are a launch screen; the first full endpoint must re-estimate all remaining storage.'}


def main():
    if not sys.platform.startswith('linux'):
        raise RuntimeError('Linux execution only')
    if not all(path.resolve() == path for path in (EX, OUT, STAGE)):
        raise RuntimeError('Redirected resource roots')
    destination = STAGE / 'RESOURCE_PROBE_ACCEPTANCE.json'
    if destination.exists():
        raise RuntimeError('Resource assessment already exists; no automatic repeat')
    folder = OUT / 'resource_probe/original_L4M4_2sweep'
    probe_path = folder / 'RESOURCE_PROBE_COMPLETE.json'
    probe = read(probe_path)
    config = read(EX / 'config/EXECUTION_CONFIG.json')
    if not (config['status'] == 'FROZEN' and config['new_fit_count'] == 120
            and config['planned_new_artifact_bytes'] == 6 * GiB and config['minimum_free_disk_bytes'] == 3 * GiB
            and probe['status'] == 'PASS' and probe['seed'] == 200611171
            and probe['ordinary_sweeps'] == 2 and probe['annealing_sweeps'] == 0
            and probe['formal_fit_budget_used'] == 0 and probe['endpoint_is_candidate'] is False):
        raise RuntimeError('Probe or configuration contract differs')
    with (OUT / 'PLAN.csv').open(encoding='utf-8-sig', newline='') as stream:
        plan = list(csv.DictReader(stream))
    if not (len(plan) == 120 and all(int(row['seed']) != probe['seed'] for row in plan)
            and not (OUT / 'fits').exists()):
        raise RuntimeError('Resource probe is not independent of an unstarted formal pool')
    with (folder / 'FIELD_OBJECT_SIZES.csv').open(newline='') as stream:
        fields = {row['field']: int(row['object_bytes']) for row in csv.DictReader(stream)}
    memory = {}
    for line in Path('/proc/meminfo').read_text().splitlines():
        if line.startswith(('MemAvailable:', 'MemTotal:')):
            memory[line.split(':')[0]] = int(line.split()[1]) * 1024
    result = estimate(config, probe, fields, tree_bytes(OUT), tree_bytes(EX),
                      shutil.disk_usage(str(OUT)).free, memory['MemAvailable'])
    result.update({'experiment_id': OUT.name, 'registered_fit_attempts': 0,
                   'independent_of_registered_pool': True, 'registered_fit_count': 120,
                   'predicted_full_endpoint_file_count': 120, 'predicted_compact_diagnostic_file_count': 120,
                   'probe_report_sha256': sha(probe_path), 'config_sha256': sha(EX / 'config/EXECUTION_CONFIG.json'),
                   'observed_unix': time.time(), 'numerical_budget_changed': False,
                   'automatic_retries': 0, 'actual_stop_requests': 0})
    with destination.open('x', encoding='utf-8') as stream:
        json.dump(result, stream, indent=2)
        stream.write('\n')
    print(json.dumps(result))
    return 0 if result['status'] == 'PASS' else 2


if __name__ == '__main__':
    sys.exit(main())
