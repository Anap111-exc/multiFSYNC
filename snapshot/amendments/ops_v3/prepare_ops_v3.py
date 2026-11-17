import json,hashlib,shutil,time
from pathlib import Path
EX=Path('/root/cavi-true222-fit444-m4-exec-20261117-v1')
OUT=Path('/root/multifsync-results/cavi_true222_fit444_m4_20261117_v1')
S=Path('/root/cavi-true222-fit444-m4-stage-20261117-v1');A=S/'amendments/ops_v3'
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
read=lambda p:json.loads(p.read_text())
def save(p,x):
    with p.open('x') as f:json.dump(x,f,indent=2);f.write('\n')
assert not (OUT/'fits').exists() and not (S/'RUNTIME_ACCEPTANCE_V3.json').exists()
for row in read(A/'AMENDMENT.json')['files']:assert sha(A/row['relative_path'])==row['sha256']
for line in (EX/'CODE_SHA256SUMS.txt').read_text().splitlines():assert sha(EX/line[66:])==line[:64]
cfg=read(A/'EXECUTION_CONFIG.json');original=read(EX/'config/EXECUTION_CONFIG.json')
keys=['planned_new_artifact_bytes','required_initial_free_disk_bytes','storage_estimate_status','operational_amendment','original_config_sha256']
assert {k:v for k,v in cfg.items() if k not in keys}=={k:v for k,v in original.items() if k not in keys}
assert cfg['minimum_free_disk_bytes']==3*1024**3 and cfg['planned_new_artifact_bytes']==27*1024**3//4
assert cfg['original_config_sha256']==sha(EX/'config/EXECUTION_CONFIG.json')
assert read(OUT/'resource_probe/original_L4M4_2sweep_v2/RESOURCE_PROBE_COMPLETE.json')['status']=='PASS'
assert read(OUT/'resource_probe/SCIENCE_SMOKE/PASS.json')['status']=='PASS'
assert read(OUT/'startup/RUNNER_ACCEPTANCE.json')['status']=='PASS'
assert read(S/'PREPARE_ACCEPTANCE.json')['status']=='PASS'
assert read(S/'STOPCHARGING_DRYRUN.json')['status']=='PASS'
assert sha(EX/'RUNTIME_LIBRARY_BINDING.json')==sha(S/'RUNTIME_LIBRARY_BINDING.json')
cleanup=[json.loads(x) for x in (S/'SAFE_RECLAIM.jsonl').read_text().splitlines()]
assert cleanup[-1]['status']=='COMPLETE' and cleanup[-1]['installed_packages_changed'] is False
old=read(S/'RESOURCE_PROBE_ACCEPTANCE.json');resource=dict(old)
assert old['status']=='BLOCKED' and old['registered_fit_attempts']==0
resource.update(status='PASS',blockers=[],artifact_budget_bytes=cfg['planned_new_artifact_bytes'],free_disk_bytes=shutil.disk_usage(str(OUT)).free,required_free_disk_bytes=cfg['required_initial_free_disk_bytes'],config_sha256=sha(A/'EXECUTION_CONFIG.json'),original_resource_assessment_sha256=sha(S/'RESOURCE_PROBE_ACCEPTANCE.json'),operational_config_path=str(A/'EXECUTION_CONFIG.json'),first_endpoint_observer_inside_pipeline_thread=True,observed_unix=time.time())
assert resource['predicted_new_artifact_bytes']<=resource['artifact_budget_bytes']
assert resource['free_disk_bytes']>=resource['required_free_disk_bytes']
assert resource['max_concurrency_from_short_probe']>=2
save(S/'RESOURCE_PROBE_ACCEPTANCE_V3.json',resource)
runtime={'status':'PASS','stopcharging_permission_dryrun':True,'startup_original_gate':True,'science_tests':True,'registration_checked':True,'runtime_library_binding_sha256':sha(S/'RUNTIME_LIBRARY_BINDING.json'),'actual_science_smoke':True,'resource_probe_acceptance':True,'formal_fit_attempts':0,'actual_stop_requests':0,'original_failure_retained':True,'inference_and_registered_budget_unchanged':True,'probe_call_correction_sha256':sha(S/'amendments/probe_call_v2/AMENDMENT.json'),'operational_amendment_sha256':sha(A/'AMENDMENT.json'),'prepared_data_sha256':sha(OUT/'GENERATION_STATUS.json')}
save(S/'RUNTIME_ACCEPTANCE_V3.json',runtime)
print(json.dumps({'status':'ALL_PRELAUNCH_ACCEPTANCE_PASS','formal_fits':0,'free_disk_bytes':resource['free_disk_bytes'],'required_free_bytes':resource['required_free_disk_bytes'],'max_concurrency_from_probe':resource['max_concurrency_from_short_probe']}))
