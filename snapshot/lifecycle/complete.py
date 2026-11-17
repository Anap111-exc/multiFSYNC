"""Seal explicitly closed scientific artifacts, then request normal StopCharging.

No whole-output directory walk, no mutable log binding, no network call on import.
The caller must wait for fitting and evaluation processes before calling this helper.
"""
import json
import os
from pathlib import Path
import time

from . import finalize as F


def artifact_record(item):
    role = item.get('role')
    F.require(role in ('terminal_status', 'scientific_results', 'failure_evidence',
                       'selection_seal', 'source_manifest'), 'Unsupported terminal artifact role')
    path = F.strict_path(item.get('path', ''), F.OUT)
    # A live controller can append to stdout after a completion manifest is made.
    F.require(path.suffix.lower() != '.log' and
              not any(x in path.name.lower() for x in ('stdout', 'stderr')),
              'Mutable stdout/stderr/log cannot be a terminal artifact')
    F.require(path.name != 'RETURN_MANIFEST.json',
              'Bind closed scientific artifacts directly, not a broad return manifest')
    before = path.stat()
    with path.open('rb') as stream:
        os.fsync(stream.fileno())
    digest = F.sha(path)
    after = path.stat()
    F.require((before.st_size, before.st_mtime_ns) == (after.st_size, after.st_mtime_ns),
              'Artifact changed while being sealed')
    return dict(path=str(path), bytes=after.st_size, sha256=digest, role=role)


def seal_and_finalize(artifacts, status, owned_process_groups,
                      process_baseline_sha256, validate_only=False,
                      inventory_fn=None, run_fn=None):
    """Return the real finalizer receipt; injected functions are for local mocks only."""
    F.require(status in ('COMPLETE', 'FAILED_SAFE_EXIT'), 'Invalid terminal status')
    F.require(isinstance(artifacts, list) and artifacts, 'Explicit artifact list required')
    F.require(isinstance(owned_process_groups, list), 'Owned worker identity list required')
    baseline_path = F.STAGE / 'PROCESS_BASELINE.json'
    baseline = F.read_bound(baseline_path, process_baseline_sha256, F.STAGE)
    F.require(baseline.get('experiment_id') == F.EXPERIMENT, 'Wrong process baseline')
    rows = (inventory_fn or F.process_inventory)()
    for group in owned_process_groups:
        F.require(set(('pid', 'start_ticks', 'pgrp')) <= set(group) and
                  group['pid'] == group['pgrp'], 'Invalid worker session identity')
        F.require(not any(r['pgrp'] == group['pgrp'] and not r['state'].startswith('Z')
                          for r in rows), 'An owned worker process group is still live')
    records = [artifact_record(item) for item in artifacts]
    roles = {record['role'] for record in records}
    F.require('terminal_status' in roles, 'Terminal status must be bound')
    need = 'scientific_results' if status == 'COMPLETE' else 'failure_evidence'
    F.require(need in roles, 'Required evidence is absent: ' + need)
    F.require(len({record['path'] for record in records}) == len(records),
              'Duplicate terminal artifact')
    receipt = F.OUT / 'LIFECYCLE_READY.json'
    F.require(not receipt.exists(), 'Lifecycle completion already sealed; no implicit retry')
    complete = dict(schema='multifsync_lifecycle_completion_v1',
        experiment_id=F.EXPERIMENT, instance_id=F.INSTANCE, region=F.REGION,
        status=status, fitting_terminal=True, evaluation_terminal=True,
        evaluation_successful=status == 'COMPLETE',
        failed_safe_exit_may_have_incomplete_science=status == 'FAILED_SAFE_EXIT',
        owned_processes_exited=True, files=records,
        owned_process_groups=owned_process_groups, created_unix=time.time(),
        local_return_required=False, active_logs_excluded=True)
    F.persist(receipt, complete, exclusive=True)
    kwargs = dict(validate_only=validate_only)
    if inventory_fn is not None:
        kwargs['inventory_fn'] = inventory_fn
    if run_fn is not None:
        kwargs['run_fn'] = run_fn
    return F.finalize(str(receipt), F.sha(receipt), str(baseline_path),
                      process_baseline_sha256, **kwargs)
