"""Local, synthetic tests only; no fit, cloud API, process signals or real /proc read."""
from pathlib import Path
import csv, hashlib, importlib.util, json, tempfile, unittest
from unittest import mock

ROOT = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('resource_observer', ROOT / 'observe_first_endpoint.py')
ob = importlib.util.module_from_spec(spec); spec.loader.exec_module(ob)


class ObserverTest(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name).resolve()
        self.out = self.root / 'out'; self.out.mkdir()
        self.ex = self.root / 'exec'; self.ex.mkdir()
        self.fit_id = 'DATA01_original_start01'
        self.folder = self.out / 'fits' / self.fit_id; self.folder.mkdir(parents=True)
        self.worker = self.out / 'scheduler' / self.fit_id; self.worker.mkdir(parents=True)
        self.manifest = 'a' * 64
        self.launch = dict(pid=123, pgrp=123, start_ticks=456)

    def tearDown(self):
        self.temp.cleanup()

    def put(self, path, value):
        path.write_text(json.dumps(value), encoding='utf-8')

    def fixture(self):
        for artifact, receipt in [('fit.rds', 'STORAGE_ROUNDTRIP.json'),
                                  ('COMPACT_DIAGNOSTICS.rds', 'COMPACT_STORAGE_ROUNDTRIP.json')]:
            target = self.folder / artifact; target.write_bytes((artifact * 80).encode())
            self.put(self.folder / receipt, dict(path=str(target), bytes=target.stat().st_size,
                sha256=ob.sha(target), object_sha256='b' * 64, identical_roundtrip=True,
                compression='xz', serialization_version=3))
        self.put(self.folder / 'SOURCE_BINDING.json', dict(fit_id=self.fit_id, source_manifest_sha256=self.manifest))
        status = dict(fit_id=self.fit_id, status='complete', t1_sweeps=400, annealing_sweeps=99,
                      fit_sha256=ob.sha(self.folder / 'fit.rds'), fit_bytes=(self.folder / 'fit.rds').stat().st_size)
        with (self.folder / 'STATUS.csv').open('w', newline='') as stream:
            writer = csv.DictWriter(stream, status.keys()); writer.writeheader(); writer.writerow(status)
        self.put(self.worker / 'EXIT.json', dict(status='SUCCESS', frozen_bindings_pass=True,
            child_still_running=False, peak_rss_bytes=800 * ob.MiB, sample_count=100,
            elapsed_seconds=300, finished_unix=1000))
        (self.worker / 'stdout.log').write_text('completed\n')

    def test_complete_endpoint_hash_receipts_and_read_only(self):
        self.fixture()
        before = {str(p): ob.sha(p) for p in self.out.rglob('*') if p.is_file()}
        value = ob.endpoint(self.out, self.fit_id, self.manifest, lambda: None)
        after = {str(p): ob.sha(p) for p in self.out.rglob('*') if p.is_file()}
        self.assertEqual(before, after)
        self.assertEqual(value['worker_peak_rss_bytes'], 800 * ob.MiB)
        self.assertFalse(value['full_storage']['independent_R_object_recheck'])

    def test_hash_mismatch_rejected(self):
        self.fixture(); (self.folder / 'fit.rds').write_bytes(b'corruption')
        with self.assertRaises(ob.AuditError):
            ob.endpoint(self.out, self.fit_id, self.manifest, lambda: None)

    def test_partial_missing_and_failed_not_success(self):
        self.assertIsNone(ob.endpoint(self.out, self.fit_id, self.manifest, lambda: None))
        self.fixture(); ended=ob.read(self.worker / 'EXIT.json'); ended['status']='NUMERICAL_FAILED_NO_RETRY'
        self.put(self.worker / 'EXIT.json', ended)
        self.assertIsNone(ob.endpoint(self.out, self.fit_id, self.manifest, lambda: None))

    def test_success_needs_complete_status_and_full_budget(self):
        self.fixture()
        path = self.folder / 'STATUS.csv'; path.write_text(path.read_text().replace('complete', 'running'))
        with self.assertRaises(ob.AuditError):
            ob.endpoint(self.out, self.fit_id, self.manifest, lambda: None)

    def test_source_binding_rejected(self):
        self.fixture()
        with self.assertRaises(ob.AuditError):
            ob.endpoint(self.out, self.fit_id, 'c' * 64, lambda: None)

    def test_partial_files_explicitly_excluded(self):
        self.fixture(); (self.folder / 'incomplete.part').write_bytes(b'x' * 99)
        result = ob.endpoint(self.out, self.fit_id, self.manifest, lambda: None)
        self.assertEqual(result['directory']['excluded_partial_bytes'], 99)
        self.assertEqual(len(result['directory']['excluded_partial_files']), 1)

    def test_remaining_forecast_includes_reserve_science_and_write_slots(self):
        self.fixture(); result = ob.endpoint(self.out, self.fit_id, self.manifest, lambda: None)
        config = dict(new_fit_count=120, max_concurrency=8, minimum_free_disk_bytes=3 * ob.GiB)
        probe = dict(per_candidate_pair_allowance_bytes=10 * ob.MiB)
        value = ob.forecast(config, probe, [result], 20 * ob.GiB)
        self.assertEqual(value['remaining_conservative_count'], 119)
        self.assertEqual(value['required_remaining_free_bytes'], (119 + 8) * 10 * ob.MiB + 512 * ob.MiB + 3 * ob.GiB)
        self.assertEqual(value['status'], 'WITHIN_PADDED_FORECAST')
        self.assertFalse(value['full_storage_upper_bound_proven'])
        self.assertEqual(ob.forecast(config, probe, [result], 1)['status'], 'INSUFFICIENT_PADDED_FORECAST')

    def test_final_status_wakes_terminal_before_identity_read(self):
        (self.out / 'FINAL_STATUS.json').write_text('{')
        def forbidden(pid): raise AssertionError('Do not need /proc after completion marker')
        self.assertEqual(ob.terminal_reason(self.out, self.launch, forbidden), 'FINAL_STATUS_EXISTS')

    def test_reused_pipeline_pid_or_exit_is_terminal(self):
        same = dict(self.launch, state='R')
        self.assertIsNone(ob.terminal_reason(self.out, self.launch, lambda pid: same))
        changed = dict(same, start_ticks=999)
        self.assertEqual(ob.terminal_reason(self.out, self.launch, lambda pid: changed), 'REGISTERED_PIPELINE_IDENTITY_CHANGED')
        self.assertEqual(ob.terminal_reason(self.out, self.launch, lambda pid: None), 'REGISTERED_PIPELINE_EXITED')

    def test_hashing_honors_completion_abort(self):
        path = self.root / 'large'; path.write_bytes(b'x' * 3 * ob.MiB)
        calls = []
        def check():
            calls.append(True)
            if len(calls) == 2: raise ob.Terminal('FINAL_STATUS_EXISTS')
        with self.assertRaises(ob.Terminal): ob.sha(path, check)
        self.assertEqual(len(calls), 2)

    def test_frozen_manifest_rechecked_and_source_change_detected(self):
        lines = []
        for i in range(36):
            path = self.ex / ('source_%02d.txt' % i); path.write_text(str(i))
            lines.append(ob.sha(path) + '  ' + path.name)
        manifest = self.ex / 'CODE_SHA256SUMS.txt'; manifest.write_text('\n'.join(lines) + '\n')
        expected = ob.sha(manifest)
        self.assertEqual(ob.verify_frozen(self.ex, expected)['files_checked'], 36)
        (self.ex / 'source_35.txt').write_text('changed')
        with self.assertRaises(ob.AuditError): ob.verify_frozen(self.ex, expected)

    def test_completion_event_interrupts_thirty_second_wait(self):
        watcher = ob.WakeOnOutput.__new__(ob.WakeOnOutput); watcher.fd = 999
        events = []
        def check():
            if events: raise ob.Terminal('FINAL_STATUS_EXISTS')
        def ready(*args):
            events.append(True); return ([999], [], [])
        with mock.patch.object(ob.select, 'select', side_effect=ready), mock.patch.object(ob.os, 'read', return_value=b'event'):
            with self.assertRaises(ob.Terminal): watcher.wait(30, check)

    def test_report_never_overwrites(self):
        path = self.out / 'resource_monitor/FIRST_ENDPOINT_STORAGE_REASSESSMENT.json'
        ob.write_report(path, {'status':'TEST'})
        with self.assertRaises(ob.AuditError): ob.write_report(path, {'status':'OVERWRITE'})


if __name__ == '__main__':
    unittest.main()
