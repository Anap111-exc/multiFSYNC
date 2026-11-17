import sys,threading
from pathlib import Path
OPS=Path(__file__).resolve().parent
sys.path.insert(0,'/root/cavi-true222-fit444-m4-exec-20261117-v1/cloud_scheduler')
"""Fit 120 original-CAVI endpoints with true222 and fitted444/M4 on ten new datasets."""
from pathlib import Path
import argparse, collections, csv, hashlib, json, os, shutil, subprocess, sys, time, traceback
from failure_policy import POLICY
from drain_owned import drain as drain_owned_groups

EXPERIMENT='cavi_true222_fit444_m4_20261117_v1'
EX=Path(os.environ.get('FSYNC_EXEC_ROOT','/root/cavi-true222-fit444-m4-exec-20261117-v1'))
OUT=Path(os.environ.get('FSYNC_OUTPUT_ROOT','/root/multifsync-results/'+EXPERIMENT))

def require(condition, message):
    if not condition:raise RuntimeError(message)

def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()

def read(path):return json.loads(Path(path).read_text(encoding='utf-8'))

def save(path,value):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('x',encoding='utf-8') as f:
        json.dump(value,f,indent=2);f.write('\n');f.flush();os.fsync(f.fileno())

def identity(pid):
    a=Path('/proc/%s/stat'%pid).read_text().rsplit(')',1)[1].split()
    return dict(pid=pid,start_ticks=int(a[19]),pgrp=int(a[2]))

def source_bindings():
    result=[]
    for line in (EX/'CODE_SHA256SUMS.txt').read_text().splitlines():
        require(len(line)>66 and line[64:66]=='  ','Bad code manifest')
        rel=Path(line[66:]);require(not rel.is_absolute() and '..' not in rel.parts,'Unsafe source path')
        path=EX/rel;require(path.resolve()==path and sha(path)==line[:64],'Source changed: '+str(rel))
        result.append(dict(path=str(path),sha256=line[:64]))
    operations=read(OPS/'AMENDMENT.json')
    for row in operations['files']:
        p=OPS/row['relative_path']
        require(p.resolve()==p and sha(p)==row['sha256'],'Operational amendment changed')
        result.append(dict(path=str(p),sha256=row['sha256']))
    result.append(dict(path=str(OPS/'AMENDMENT.json'),sha256=sha(OPS/'AMENDMENT.json')))
    require(len(result)>35,'Incomplete source binding')
    return result

def registered_plan():
    with (EX/'config/PLAN.csv').open(encoding='utf-8-sig',newline='') as f:rows=list(csv.DictReader(f))
    require(len(rows)==120 and len({r['fit_id'] for r in rows})==120,'Expected 120 registered fits')
    require({r['arm'] for r in rows}=={'original'},'Expected the registered original CAVI method')
    require(sorted(collections.Counter(r['dataset_id'] for r in rows).values())==[12]*10,'Expected ten datasets')
    for dataset in {r['dataset_id'] for r in rows}:
        paired={}
        for arm in ('original',):
            arm_rows=[r for r in rows if r['dataset_id']==dataset and r['arm']==arm]
            require(sorted(int(r['start']) for r in arm_rows)==list(range(1,13)),'Start registry differs')
            paired[arm]={int(r['start']):int(r['seed']) for r in arm_rows}
    require([r['arm'] for r in rows[:2]]==['original','original'],'Resource pilot requires original CAVI starts')
    require(rows[0]['dataset_id']==rows[1]['dataset_id'] and [int(r['start']) for r in rows[:2]]==[1,2],'Resource pilot requires first two registered starts')
    return rows

def execution_config():
    config=read(OPS/'EXECUTION_CONFIG.json')
    require(config['status']=='FROZEN' and config['new_fit_count']==120,'Execution configuration is not frozen')
    require(config['fitted_factor_capacities']==[4,4,4] and config['fitted_fpca_capacity']==4,'Unexpected capacity')
    require(config['annealing_sweeps']==99 and config['ordinary_t1_sweeps']==400,'Update budget changed')
    for key in ('minimum_free_disk_bytes','planned_new_artifact_bytes'):
        require(type(config[key]) is int and config[key]>0,'Invalid disk budget: '+key)
    return config

def runtime_env():
    env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1',FSYNC_EXEC_ROOT=str(EX),FSYNC_OUTPUT_ROOT=str(OUT))
    dep=read(EX/'RUNTIME_LIBRARY_BINDING.json')
    env['R_LIBS_USER']=dep['R_LIBS_USER']
    env['R_LIBS']=dep['R_LIBS']
    for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):env[key]='1'
    return env

def freeze_scheduler(rows):
    g=read(OUT/'GENERATION_STATUS.json')
    require(g['status']=='PASS' and g['new_datasets']==10 and g['reused_datasets']==0 and g['registered_fits']==120,'Data preparation failed')
    require((OUT/'DATA_READY.txt').is_file() and sha(OUT/'PLAN.csv')==sha(EX/'config/PLAN.csv'),'PLAN changed')
    folder=OUT/'fit_plan_v1';folder.mkdir(exist_ok=False)
    bindings=source_bindings()+[dict(path=str(OUT/'PLAN.csv'),sha256=sha(OUT/'PLAN.csv'))]
    for dataset in sorted({r['dataset_id'] for r in rows}):
        directory=OUT/'fit_inputs'/dataset;manifest=directory/'INPUT_SHA256SUMS.txt'
        bindings.append(dict(path=str(manifest),sha256=sha(manifest)))
        for line in manifest.read_text().splitlines():
            require(len(line)>66 and line[64:66]=='  ','Bad input manifest')
            rel=Path(line[66:]);require(not rel.is_absolute() and '..' not in rel.parts and 'truth' not in str(rel).lower(),'Forbidden fit input')
            path=directory/rel;require(path.resolve()==path and sha(path)==line[:64],'Input changed')
            bindings.append(dict(path=str(path),sha256=line[:64]))
    save(folder/'BINDINGS.json',bindings)
    tasks=[dict(id=r['fit_id'],dataset_id=r['dataset_id'],start=int(r['start']),arm=r['arm'],seed=int(r['seed']),
                command=['Rscript','--vanilla',str(EX/'scripts/run_cloud_fit.R'),r['fit_id']],cwd=str(EX),output_dir=str(OUT/'fits'/r['fit_id'])) for r in rows]
    save(folder/'PLAN.json',dict(schema=1,experiment_id=EXPERIMENT,execution_root=str(EX),output_root=str(OUT),tasks=tasks,bindings=bindings))
    config=execution_config()
    release=dict(schema=1,status='RELEASED',experiment_id=EXPERIMENT,source_frozen=True,generation_pass=True,
        protocol_path=str(EX/'config/PROTOCOL.json'),protocol_sha256=sha(EX/'config/PROTOCOL.json'),
        execution_config_path=str(OPS/'EXECUTION_CONFIG.json'),execution_config_sha256=sha(OPS/'EXECUTION_CONFIG.json'),
        task_plan_sha256=sha(folder/'PLAN.json'),binding_manifest_path=str(folder/'BINDINGS.json'),binding_manifest_sha256=sha(folder/'BINDINGS.json'),
        approved_task_count=120,approved_datasets=10,approved_independent_populations=10,approved_starts_per_dataset=12,
        approved_arms=['original'],initial_concurrency=2,max_concurrency=8,pilot_observation_seconds=300,
        system_reserve_bytes=3*1024**3,min_disk_free_bytes=config['minimum_free_disk_bytes'],
        numeric_failure_policy=POLICY,pilot_policy='first_two_registered_original_starts_300s_v1',fit_timeout_seconds=None,
        automatic_retries=0,automatic_shutdown_deployed=True,timed_shutdown=False)
    save(folder/'RELEASE.json',release)
    return folder,sha(folder/'RELEASE.json')

def fit_phase_manifest(rows,completion):
    records=[]
    for row in rows:
        for name in ('STATUS.csv','fit.rds'):
            p=OUT/'fits'/row['fit_id']/name
            if p.is_file():records.append(dict(path=str(p),bytes=p.stat().st_size,sha256=sha(p)))
    save(OUT/'FIT_PHASE_COMPLETE.json',dict(status='ALL_REGISTERED_ATTEMPTS_TERMINAL',planned=120,
        plan_sha256=sha(OUT/'PLAN.csv'),scheduler_completion=completion,artifacts=records))

def scientific_manifest():
    records=[]
    # Closed scientific outputs only; controller logs and lifecycle receipts remain mutable.
    for name in ('selection','science','aggregate'):
        for p in sorted((OUT/name).rglob('*')):
            if p.is_file():
                require(p.resolve()==p,'Redirected scientific output')
                records.append(dict(path=p.relative_to(OUT).as_posix(),bytes=p.stat().st_size,sha256=sha(p)))
    require(records,'No scientific output to seal')
    save(OUT/'SCIENTIFIC_MANIFEST.json',dict(experiment_id=EXPERIMENT,files=records,live_logs_excluded=True))

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--validate-only',action='store_true');args=parser.parse_args()
    require(sys.platform.startswith('linux'),'Execution requires Linux')
    require(EX.resolve()==EX and OUT.resolve()==OUT,'Redirected root')
    source_bindings();rows=registered_plan();config=execution_config()
    for rel in ('generation/generate_registered.R','scripts/run_cloud_fit.R','science/evaluate_dataset.R',
                'science/summarize.py','lifecycle/complete.py','lifecycle/finalize.py','config/PROTOCOL_FROZEN.json'):
        require((EX/rel).is_file(),'Missing entry '+rel)
    require(read(OUT/'startup/RUNNER_ACCEPTANCE.json').get('status')=='PASS','Original-CAVI startup equivalence acceptance absent')
    if args.validate_only:
        print(json.dumps(dict(status='VALIDATED_NOT_STARTED',datasets=10,fits=120,arms=['original'],automatic_retries=0,fit_timeout_seconds=None)));return 0
    log=OUT/'pipeline';log.mkdir(parents=True,exist_ok=False)
    save(log/'START.json',dict(identity=identity(os.getpid()),started_unix=time.time(),global_wall_clock_deadline=None))
    env=runtime_env();stages=[];failure=None;resource_thread=None
    def run(name,command):
        source_bindings();start=time.monotonic()
        with (log/(name+'.stdout.log')).open('xb') as so,(log/(name+'.stderr.log')).open('xb') as se:
            child=subprocess.Popen(command,cwd=str(EX),env=env,stdin=subprocess.DEVNULL,stdout=so,stderr=se)
            try:
                save(log/(name+'.START.json'),dict(identity=identity(child.pid),command=command,started_unix=time.time()))
                code=child.wait()
            finally:
                # Even an identity/receipt failure must not leave a live stage
                # when the completion handler considers normal instance stop.
                if child.poll() is None:
                    child.wait()
        result=dict(name=name,exit_code=code,wall_seconds=time.monotonic()-start)
        save(log/(name+'.EXIT.json'),result);stages.append(result)
        require(code==0,'Stage failed: '+name+'; no automatic retry');source_bindings()
    try:
        free=shutil.disk_usage(OUT).free
        required_free=config['minimum_free_disk_bytes']+config['planned_new_artifact_bytes']
        save(log/'DISK_GATE.json',dict(free_bytes=free,floor_bytes=config['minimum_free_disk_bytes'],
            new_artifact_allowance_bytes=config['planned_new_artifact_bytes'],required_free_bytes=required_free))
        require(free>=required_free,'Insufficient space for the frozen artifact allowance and disk reserve')
        prepared=read(OUT/'GENERATION_STATUS.json')
        require(prepared['status']=='PASS' and prepared['new_datasets']==10 and prepared['reused_datasets']==0 and prepared['registered_fits']==120,'Pre-launch data preparation did not pass')
        save(log/'VERIFIED_PRELAUNCH_GENERATION.json',dict(status='PASS',generation_status_sha256=sha(OUT/'GENERATION_STATUS.json'),new_data_generation=10))
        import observe_first_endpoint_thread as resource_monitor
        resource_thread=threading.Thread(target=resource_monitor.main,args=([
            '--source-manifest-sha256',sha(EX/'CODE_SHA256SUMS.txt'),
            '--observer-source-sha256',sha(OPS/'observe_first_endpoint_thread.py')],),daemon=True)
        resource_thread.start()
        plan,release_sha=freeze_scheduler(rows)
        run('fits',['python3','-B',str(EX/'cloud_scheduler/schedule.py'),'--plan',str(plan/'PLAN.json'),
                    '--release',str(plan/'RELEASE.json'),'--release-sha256',release_sha,'--execute'])
        completion=read(OUT/'scheduler/COMPLETION.json')
        require(completion['status'] in ('SUCCESS','COMPLETE_WITH_FAILURES') and completion['all_registered_attempted']
                and completion['total_registered']==120 and completion['all_known_owned_groups_exited'] and completion['safe_completion_claim'],
                'Fitting is not safely terminal')
        fit_phase_manifest(rows,completion)
        science_failures=[]
        for dataset in sorted({r['dataset_id'] for r in rows}):
            try:
                run(dataset+'_science',['Rscript','--vanilla',str(EX/'science/evaluate_dataset.R'),'--exec-root',str(EX),'--output-root',str(OUT),'--data-id',dataset])
                marker=read(OUT/'science'/dataset/'SCIENCE_COMPLETE.json')
                require(marker.get('status')=='SCIENCE_COMPLETE','Missing science completion marker: '+dataset)
            except Exception as error:
                science_failures.append(dict(dataset_id=dataset,error_type=type(error).__name__,message=str(error)))
        if science_failures:
            save(OUT/'SCIENCE_FAILURES.json',dict(failed_datasets=science_failures,automatic_retries=0))
            raise RuntimeError('Scientific evaluation failures; other pools evaluated; no automatic retries')
        run('aggregate',['python3','-B',str(EX/'science/summarize.py'),str(OUT)])
        require((OUT/'aggregate/SUMMARY.json').is_file(),'Science summary absent')
        scientific_manifest()
    except BaseException as error:
        failure=dict(type=type(error).__name__,message=str(error))
        with (log/'FAILURE.txt').open('x',encoding='utf-8') as f:f.write(traceback.format_exc());f.flush();os.fsync(f.fileno())
    drain_result=None
    if failure is not None:
        drain_result=drain_owned_groups(OUT,[r['fit_id'] for r in rows],reason=failure)
    final=dict(status='PASS' if failure is None else 'FAIL',experiment_id=EXPERIMENT,stages=stages,failure=failure,
        registered_fits=120,automatic_retries=0,fit_timeout_seconds=None,automatic_shutdown_deployed=True,
        timed_shutdown=False,local_return_required_before_stop=False,finished_unix=time.time(),
        scientific_acceptability='DESCRIPTIVE_PILOT_NOT_AUTOMATIC_ACCEPTANCE')
    save(OUT/'FINAL_STATUS.json',final)
    if resource_thread is not None:resource_thread.join()
    groups=[]
    for p in sorted((OUT/'scheduler').glob('*/WORKER_IDENTITY.json')):
        v=read(p);groups.append({k:v[k] for k in ('pid','start_ticks','pgrp')})
    artifacts=[dict(role='terminal_status',path=str(OUT/'FINAL_STATUS.json'))]
    if failure is None:
        artifacts.extend([dict(role='scientific_results',path=str(OUT/'aggregate/SUMMARY.json')),
                          dict(role='scientific_results',path=str(OUT/'SCIENTIFIC_MANIFEST.json'))])
    else:artifacts.append(dict(role='failure_evidence',path=str(log/'FAILURE.txt')))
    print(json.dumps(final,indent=2),flush=True)
    sys.path.insert(0,str(EX))
    from lifecycle.complete import seal_and_finalize
    if drain_result is not None and drain_result['status']!='DRAINED':
        save(OUT/'STOP_BLOCKED_UNRESOLVED_OWNED_GROUPS.json',drain_result)
        return 3
    seal_and_finalize(artifacts=artifacts,status='COMPLETE' if failure is None else 'FAILED_SAFE_EXIT',
                      owned_process_groups=groups,process_baseline_sha256=os.environ.get('FSYNC_PROCESS_BASELINE_SHA256',''))
    return 0 if failure is None else 2

if __name__=='__main__':raise SystemExit(main())
