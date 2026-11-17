from pathlib import Path
import unittest,sys,json,tempfile,csv,hashlib,copy,importlib.util
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'cloud_scheduler'))
import schedule,resource_policy as rp

class RegistrationTest(unittest.TestCase):
 def setUp(self):
  self.temp=tempfile.TemporaryDirectory();self.base=Path(self.temp.name)
  self.code=self.base/'code.txt';self.code.write_text('frozen-test-only')
  self.protocol=self.base/'PROTOCOL.json';self.protocol.write_text('{}')
  self.config=self.base/'EXECUTION_CONFIG.json';self.config.write_text(json.dumps({'status':'FROZEN','new_fit_count':120,'minimum_free_disk_bytes':3*rp.GiB}))
  with (ROOT/'config/PLAN.csv').open(encoding='utf-8') as f:registry=list(csv.DictReader(f))
  self.tasks=[{'id':r['fit_id'],'dataset_id':r['dataset_id'],'arm':r['arm'],'seed':int(r['seed']),'start':int(r['start']),'cwd':str(self.base),'output_dir':str(self.base/r['fit_id']),'command':['Rscript','registered-test-only']} for r in registry]
  self.bindings=[{'path':str(self.code),'sha256':schedule.sha(self.code)}]
  self.manifest=self.base/'BINDINGS.json';self.manifest.write_text(json.dumps(self.bindings))
  self.plan={'experiment_id':ROOT.name,'execution_root':str(self.base),'output_root':str(self.base),'tasks':self.tasks,'bindings':self.bindings}
  self.release={'status':'RELEASED','experiment_id':ROOT.name,'source_frozen':True,'generation_pass':True,'protocol_path':str(self.protocol),'protocol_sha256':schedule.sha(self.protocol),'binding_manifest_path':str(self.manifest),'binding_manifest_sha256':schedule.sha(self.manifest),'approved_task_count':120,'approved_datasets':10,'approved_starts_per_dataset':12,'approved_arms':['original'],'approved_independent_populations':10,'initial_concurrency':2,'max_concurrency':8,'pilot_observation_seconds':300,'system_reserve_bytes':3*rp.GiB,'execution_config_path':str(self.config),'execution_config_sha256':schedule.sha(self.config),'min_disk_free_bytes':3*rp.GiB,'numeric_failure_policy':schedule.POLICY,'pilot_policy':'first_two_registered_original_starts_300s_v1','fit_timeout_seconds':None}
 def tearDown(self):self.temp.cleanup()
 def validate(self,p=None):
  p=self.plan if p is None else p;pp=self.base/'PLAN.json';pp.write_text(json.dumps(p));r=copy.deepcopy(self.release);r['task_plan_sha256']=schedule.sha(pp);rr=self.base/'RELEASE.json';rr.write_text(json.dumps(r));return schedule.validate(pp,rr,schedule.sha(rr))
 def test_actual_registry_accepts_exactly_120_original(self):self.assertEqual(len(self.validate()[0]['tasks']),120)
 def test_rejects_extra_mgp(self):
  p=copy.deepcopy(self.plan);p['tasks'][5]['arm']='mgp'
  with self.assertRaises(AssertionError):self.validate(p)
 def test_rejects_missing_or_duplicate_start(self):
  p=copy.deepcopy(self.plan);p['tasks'][3]['start']=2
  with self.assertRaises(AssertionError):self.validate(p)
 def test_rejects_modified_frozen_source(self):
  self.code.write_text('changed')
  with self.assertRaises(AssertionError):self.validate()
 def test_resource_caps_and_reserve(self):
  args={'available_bytes':15*rp.GiB,'owned_rss_bytes':0,'running_budgets':[rp.GiB]*2,'candidate_peak_bytes':100*1024**2,'pilot_ready':False,'running_count':2}
  self.assertFalse(rp.can_submit(**args));args['pilot_ready']=True;self.assertTrue(rp.can_submit(**args));args['running_count']=8;self.assertFalse(rp.can_submit(**args));args['running_count']=2;args['available_bytes']=2*rp.GiB;self.assertFalse(rp.can_submit(**args))
 def test_single_arm_pilot_needs_both_observed(self):
  ids=[t['id'] for t in self.tasks[:2]];state={i:{'elapsed_seconds':301,'sample_count':150,'monitoring_errors':[]} for i in ids}
  self.assertTrue(rp.pilot_ready(state,self.tasks));state.pop(ids[1]);self.assertFalse(rp.pilot_ready(state,self.tasks))

if __name__=='__main__':unittest.main()
