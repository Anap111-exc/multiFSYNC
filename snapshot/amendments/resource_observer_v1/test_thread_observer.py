"""Inherited resource checks plus bounded launch wait; no pipeline or cloud execution."""
from pathlib import Path
import importlib.util, inspect, json, unittest
from unittest import mock

ROOT = Path(__file__).resolve().parent
def load(name, file):
    spec=importlib.util.spec_from_file_location(name, ROOT/file)
    module=importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    return module

ob=load('thread_observer', 'observe_first_endpoint_thread.py')
base=load('observer_shared_tests', 'test_observer.py')
base.ob=ob


class SharedThreadObserverTest(base.ObserverTest):
    pass


class ThreadEntryTest(unittest.TestCase):
    def test_signature_accepts_explicit_argument_vector(self):
        signature=inspect.signature(ob.main)
        self.assertEqual(list(signature.parameters), ['argv'])
        self.assertIsNone(signature.parameters['argv'].default)

    def test_already_started_no_sleep(self):
        value={'status':'STARTED','pid':123}
        with mock.patch.object(ob, 'read', return_value=value):
            result=ob.await_started_launch('unused', sleeper=lambda seconds:self.fail('Unnecessary sleep'))
        self.assertEqual(result, value)

    def test_waits_for_missing_partial_and_pending_receipt(self):
        now=[0.0]
        def sleep(seconds): now[0]+=seconds
        value={'status':'STARTED','pid':123}
        failures=[FileNotFoundError(),json.JSONDecodeError('partial','{',1),{'status':'PENDING'},value]
        with mock.patch.object(ob,'read',side_effect=failures):
            result=ob.await_started_launch('unused', clock=lambda:now[0], sleeper=sleep)
        self.assertEqual(result,value)
        self.assertAlmostEqual(now[0],.15)

    def test_missing_launch_exits_within_five_seconds(self):
        now=[0.0]
        def sleep(seconds): now[0]+=seconds
        with mock.patch.object(ob,'read',side_effect=FileNotFoundError()):
            with self.assertRaises(ob.AuditError):
                ob.await_started_launch('unused', clock=lambda:now[0], sleeper=sleep)
        self.assertLessEqual(now[0],5.0)
        self.assertGreaterEqual(now[0],4.99)

    def test_main_thread_is_rejected_without_any_cloud_or_filesystem_action(self):
        with mock.patch.object(ob.sys,'platform','linux'):
            with self.assertRaises(ob.AuditError):
                ob.main(['--observer-source-sha256','a'*64,'--source-manifest-sha256','b'*64])


if __name__=='__main__':
    unittest.main()
