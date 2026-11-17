"""Linux durable 120-fit scheduler. Validation is default; --execute is explicit.

The supervisor never compares objectives, retries a task, or stops unrelated work.
Workers persist their own exit records and survive supervisor disconnection.
"""
from pathlib import Path
import argparse,collections,csv,hashlib,json,os,shutil,signal,subprocess,sys,time
from resource_policy import GiB,can_submit,pilot_witnesses,pilot_ready,candidate_peak,memory_budget
from failure_policy import POLICY,classify_failure

def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    return h.hexdigest()
def read(path):return json.loads(Path(path).read_text(encoding='utf-8'))
def write(path,obj):
    path=Path(path);tmp=path.with_suffix(path.suffix+'.tmp')
    with tmp.open('w',encoding='utf-8') as f:
        json.dump(obj,f,indent=2);f.write('\n');f.flush();os.fsync(f.fileno())
    os.replace(tmp,path)
def within(path,root):
    try:Path(path).resolve().relative_to(Path(root).resolve());return True
    except ValueError:return False
def proc_identity(pid):
    try:
        txt=Path(f'/proc/{pid}/stat').read_text();parts=txt[txt.rfind(')')+2:].split()
        return dict(pid=int(pid),state=parts[0],pgrp=int(parts[2]),start_ticks=int(parts[19]),rss_bytes=int(parts[21])*os.sysconf('SC_PAGE_SIZE'))
    except (FileNotFoundError,ProcessLookupError):return None
def group_rss(pgrp):
    total=0
    for p in Path('/proc').iterdir():
        if p.name.isdigit():
            row=proc_identity(int(p.name))
            if row and row['pgrp']==pgrp:total+=row['rss_bytes']
    return total
def group_live(pgrp):
    for p in Path('/proc').iterdir():
        if p.name.isdigit():
            row=proc_identity(int(p.name))
            if row and row['pgrp']==pgrp and row['state']!='Z':return True
    return False
def available_memory():
    for line in Path('/proc/meminfo').read_text().splitlines():
        if line.startswith('MemAvailable:'):return int(line.split()[1])*1024
    raise RuntimeError('Linux MemAvailable is required')
def owned_alive(identity):
    now=proc_identity(identity['pid'])
    return bool(now and now['state']!='Z' and now['start_ticks']==identity['start_ticks'] and now['pgrp']==identity['pgrp'])

def validate(planpath,releasepath,release_sha):
    p=read(planpath);r=read(releasepath)
    assert sha(releasepath)==release_sha,'Independent release hash mismatch'
    assert r['status']=='RELEASED' and r['task_plan_sha256']==sha(planpath)
    assert p['experiment_id']==r['experiment_id']=='cavi_true222_fit444_m4_20261117_v1'
    assert r['source_frozen'] is True and r['generation_pass'] is True
    assert sha(r['protocol_path'])==r['protocol_sha256']
    assert sha(r['binding_manifest_path'])==r['binding_manifest_sha256']
    assert (r['approved_task_count'],r['approved_datasets'],r['approved_starts_per_dataset'])==(120,10,12)
    assert r['approved_arms']==['original']
    assert r['approved_independent_populations']==10
    assert (r['initial_concurrency'],r['max_concurrency'],r['pilot_observation_seconds'],r['system_reserve_bytes'])==(2,8,300,3*GiB)
    assert sha(r['execution_config_path'])==r['execution_config_sha256']
    config=read(r['execution_config_path'])
    assert config['status']=='FROZEN' and config['new_fit_count']==120
    assert type(config['minimum_free_disk_bytes']) is int and config['minimum_free_disk_bytes']>0
    assert r['min_disk_free_bytes']==config['minimum_free_disk_bytes']
    assert r['numeric_failure_policy']==POLICY
    assert r['pilot_policy']=='first_two_registered_original_starts_300s_v1'
    assert r['fit_timeout_seconds'] is None # user explicitly removed wall-clock cutoffs
    ts=p['tasks'];assert len(ts)==120 and len({t['id'] for t in ts})==120
    assert {t['arm'] for t in ts}=={'original'}
    assert sorted(collections.Counter(t['dataset_id'] for t in ts).values())==[12]*10
    for dataset in {t['dataset_id'] for t in ts}:
        seeds={}
        for arm in ('original',):
            arm_tasks=[t for t in ts if t['dataset_id']==dataset and t['arm']==arm]
            assert sorted(t['start'] for t in arm_tasks)==list(range(1,13))
            seeds[arm]={t['start']:t['seed'] for t in arm_tasks}
    assert [t['arm'] for t in ts[:2]]==['original','original']
    assert ts[0]['dataset_id']==ts[1]['dataset_id'] and [t['start'] for t in ts[:2]]==[1,2]
    assert len({t['output_dir'] for t in ts})==120
    for t in ts:
        assert t['id'] and all(c.isalnum() or c in '_-.' for c in t['id'])
        assert isinstance(t['command'],list) and t['command'] and all(isinstance(s,str) for s in t['command'])
        assert within(t['cwd'],p['execution_root']) and within(t['output_dir'],p['output_root'])
    assert p['bindings'] and all(sha(b['path'])==b['sha256'] for b in p['bindings'])
    manifest=read(r['binding_manifest_path'])
    assert manifest==p['bindings'],'Release manifest must equal the exact plan binding list'
    return p,r

def worker(args):
    p,r=validate(args.plan,args.release,args.release_sha256)
    task=next(t for t in p['tasks'] if t['id']==args.worker);folder=Path(args.worker_root)
    assert folder.exists() and not (folder/'EXIT.json').exists()
    start=time.monotonic();ident=proc_identity(os.getpid());assert ident['pgrp']==os.getpid()
    write(folder/'WORKER_IDENTITY.json',dict(**ident,started_unix=time.time(),plan_sha256=sha(args.plan),release_sha256=args.release_sha256))
    env=os.environ.copy()
    env.update(FSYNC_EXEC_ROOT=p['execution_root'],FSYNC_OUTPUT_ROOT=p['output_root'])
    for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):env[k]='1'
    tmp=Path(p['execution_root'])/'tmp'/task['id'];tmp.mkdir(parents=True,exist_ok=False)
    for k in ('TMPDIR','TMP','TEMP'):env[k]=str(tmp)
    peak=ident['rss_bytes'];samples=1;code=None;failure=None;child=None;monitoring_errors=[];rss=ident['rss_bytes']
    try:
        with (folder/'stdout.log').open('wb') as so,(folder/'stderr.log').open('wb') as se:
            child=subprocess.Popen(task['command'],cwd=task['cwd'],env=env,stdout=so,stderr=se)
            write(folder/'FIT_IDENTITY.json',dict(pid=child.pid,proc=proc_identity(child.pid),command=task['command'],started_unix=time.time()))
            while child.poll() is None:
                try:rss=group_rss(ident['pgrp']);peak=max(peak,rss);samples+=1;available=available_memory()
                except Exception as e:
                    monitoring_errors.append(f'{type(e).__name__}: {e}');available=None
                    failure='Resource monitoring failed; waited for the child to exit normally'
                write(folder/'RESOURCE.json',dict(status='RUNNING',elapsed_seconds=time.monotonic()-start,
                    sample_count=samples,current_rss_bytes=rss,peak_rss_bytes=peak,available_memory_bytes=available,
                    monitoring_errors=monitoring_errors))
                time.sleep(2)
            code=child.wait()
        if not all(sha(b['path'])==b['sha256'] for b in p['bindings']):failure='Frozen binding changed'
    except Exception as e:
        failure=f'{type(e).__name__}: {e}'
        # Do not orphan a still-running R child after a monitoring/write error.
        # Default no-timeout semantics preserve its normal completion.
        if child is not None:code=child.wait()
    binding_ok=all(sha(b['path'])==b['sha256'] for b in p['bindings'])
    status_rows=[];status_sha=None;status_error=None;status_file=Path(task['output_dir'])/'STATUS.csv'
    if code==2 and status_file.exists():
        try:
            with status_file.open(encoding='utf-8-sig',newline='') as f:status_rows=list(csv.DictReader(f))
            status_sha=sha(status_file)
        except Exception as e:status_error=f'{type(e).__name__}: {e}'
    classification='SUCCESS' if code==0 and failure is None and binding_ok else classify_failure(
        exit_code=code,fit_id=task['id'],status_rows=status_rows,bindings_ok=binding_ok,
        monitoring_errors=monitoring_errors,worker_error=failure or status_error)
    result=dict(status=classification,exit_code=code,error=failure,
        elapsed_seconds=time.monotonic()-start,sample_count=samples,current_rss_bytes=0,peak_rss_bytes=peak,
        retries=0,monitoring_errors=monitoring_errors,child_still_running=child is not None and child.poll() is None,
        numeric_failure_policy=POLICY,status_csv_sha256=status_sha,status_rows=status_rows,status_read_error=status_error,
        frozen_bindings_pass=binding_ok,
        plan_sha256=sha(args.plan),release_sha256=args.release_sha256,finished_unix=time.time())
    write(folder/'RESOURCE.json',result);write(folder/'EXIT.json',result)
    return 0 if result['status']=='SUCCESS' else 2

def supervise(args,p,r):
    if sys.platform!='linux':raise RuntimeError('Execution is Linux-only')
    root=Path(p['output_root'])/'scheduler';ledger=root/'LEDGER.json'
    if ledger.exists():
        assert args.resume,'Existing ledger; explicit --resume required (never retry)'
        state=read(ledger);assert state['plan_sha256']==sha(args.plan) and state['release_sha256']==args.release_sha256
    else:
        assert not args.resume;root.mkdir(parents=True,exist_ok=False)
        state=dict(plan_sha256=sha(args.plan),release_sha256=args.release_sha256,
            tasks={t['id']:dict(status='PENDING',attempts=0) for t in p['tasks']},stop_submitting=False,created_unix=time.time())
        write(ledger,state)
    pilot_ids=[t['id'] for t in p['tasks'][:2]];task_by_id={t['id']:t for t in p['tasks']}
    # Two pilot slots need the minimum two-task budget plus the system reserve.
    if not any(s['attempts'] for s in state['tasks'].values()) and available_memory()<5*GiB:
        write(root/'RESOURCE_GATE_WAIT.json',dict(available=available_memory(),required=5*GiB));return 3
    while True:
        resources={};running=[]
        for id,s in state['tasks'].items():
            folder=root/id
            if (folder/'RESOURCE.json').exists():
                resources[id]=read(folder/'RESOURCE.json')
                if resources[id].get('monitoring_errors'):state['stop_submitting']=True
            if s['status'] not in ('RUNNING','STARTING','TIMEOUT_REQUESTED','ORPHAN_LIVE'):continue
            if not s.get('identity') and (folder/'WORKER_IDENTITY.json').exists():
                identity=read(folder/'WORKER_IDENTITY.json')
                assert identity['plan_sha256']==sha(args.plan) and identity['release_sha256']==args.release_sha256
                s['identity']=identity
            live_group=bool(s.get('identity') and group_live(s['identity']['pgrp']))
            if (folder/'EXIT.json').exists():
                done=read(folder/'EXIT.json');s['exit']=done
                if not s.get('timeout_seconds') and done['status'] not in ('SUCCESS','NUMERICAL_FAILED_NO_RETRY'):state['stop_submitting']=True
                if live_group:
                    s['exit_waiting_for_group']=True;running.append(id)
                else:s.update(status='TIMED_OUT' if s.get('timeout_seconds') else done['status'],exit_waiting_for_group=False)
            elif live_group:
                running.append(id)
                if not owned_alive(s['identity']):s['status']='ORPHAN_LIVE';state['stop_submitting']=True
            else:
                s['status']='TIMED_OUT' if s.get('timeout_seconds') else 'UNRESOLVED_EXIT'
                s['error']='No durable exit and no matching live owned worker; never retried'
                state['stop_submitting']=True
                if not s.get('identity'):state['unresolved_unowned_launch']=True
        pilot_ok=pilot_ready(resources,p['tasks'],300)
        state['pilot_witnesses']=pilot_witnesses(resources,p['tasks'],300)
        rss=sum(group_rss(state['tasks'][id]['identity']['pgrp']) for id in running)
        state['pilot_ready']=pilot_ok;state['available_memory_bytes']=available_memory();state['owned_rss_bytes']=rss
        # No elapsed-time signal is sent: frozen 99+400 iteration budgets govern each fit.
        if not state['stop_submitting']:
            for task in p['tasks']:
                id=task['id'];s=state['tasks'][id]
                if s['status']!='PENDING':continue
                disk=shutil.disk_usage(p['output_root']).free
                state['free_disk_bytes']=disk
                if disk<r['min_disk_free_bytes']:
                    state['stop_submitting']=True;state['stop_reason']='DISK_RESOURCE_GATE';break
                budgets=[max(memory_budget(resources.get(j,{}).get('peak_rss_bytes',0)),
                    state['tasks'][j].get('allocated_budget_bytes',GiB)) for j in running]
                if not can_submit(available_bytes=available_memory(),owned_rss_bytes=sum(group_rss(state['tasks'][j]['identity']['pgrp']) for j in running),
                    running_budgets=budgets,candidate_peak_bytes=candidate_peak(task,resources,p['tasks']),
                    pilot_ready=pilot_ok,running_count=len(running)):break
                assert s['attempts']==0
                folder=root/id;folder.mkdir(exist_ok=False)
                s.update(status='STARTING',attempts=1,started_unix=time.time(),
                    allocated_budget_bytes=memory_budget(candidate_peak(task,resources,p['tasks'])));write(ledger,state)
                command=[sys.executable,str(Path(__file__).resolve()),'--plan',str(Path(args.plan).resolve()),
                    '--release',str(Path(args.release).resolve()),'--release-sha256',args.release_sha256,
                    '--worker',id,'--worker-root',str(folder)]
                with (folder/'worker_stdout.log').open('wb') as so,(folder/'worker_stderr.log').open('wb') as se:
                    proc=subprocess.Popen(command,stdout=so,stderr=se,start_new_session=True)
                identity=proc_identity(proc.pid)
                if not identity or identity['pgrp']!=proc.pid:raise RuntimeError('Worker identity could not be persisted; no retry')
                s.update(status='RUNNING',identity=identity);running.append(id)
                state['max_concurrent_observed']=max(state.get('max_concurrent_observed',0),len(running));write(ledger,state)
        write(ledger,state)
        if state['stop_submitting'] and not running:
            for s in state['tasks'].values():
                if s['status']=='PENDING':s['status']='UNATTEMPTED_RESOURCE' if state.get('stop_reason')=='DISK_RESOURCE_GATE' else 'UNATTEMPTED_STOPPED'
            state['all_known_owned_groups_exited']=True
            state['safe_completion_claim']=not state.get('unresolved_unowned_launch',False)
            write(ledger,state)
            write(root/'COMPLETION.json',dict(status='FAILED',task_counts=dict(collections.Counter(s['status'] for s in state['tasks'].values())),
                all_known_owned_groups_exited=True,safe_completion_claim=state['safe_completion_claim'],finished_unix=time.time()))
            return 2
        if all(s['status'] in ('SUCCESS','NUMERICAL_FAILED_NO_RETRY','TIMED_OUT') for s in state['tasks'].values()):
            state['all_known_owned_groups_exited']=True;state['safe_completion_claim']=True;write(ledger,state)
            counts=dict(collections.Counter(s['status'] for s in state['tasks'].values()))
            completion='COMPLETE_WITH_FAILURES' if counts.get('NUMERICAL_FAILED_NO_RETRY',0) or counts.get('TIMED_OUT',0) else 'SUCCESS'
            write(root/'COMPLETION.json',dict(status=completion,task_counts=counts,all_registered_attempted=True,
                total_registered=120,automatic_retries=0,all_known_owned_groups_exited=True,
                safe_completion_claim=True,max_concurrent_observed=state.get('max_concurrent_observed'),finished_unix=time.time()))
            return 0
        time.sleep(5)

def main():
    a=argparse.ArgumentParser();a.add_argument('--plan',required=True);a.add_argument('--release',required=True)
    a.add_argument('--release-sha256',required=True);a.add_argument('--execute',action='store_true');a.add_argument('--resume',action='store_true')
    a.add_argument('--worker');a.add_argument('--worker-root');args=a.parse_args()
    if args.worker:return worker(args)
    p,r=validate(args.plan,args.release,args.release_sha256)
    if not args.execute:
        print(json.dumps(dict(status='VALIDATED_NOT_RUN',tasks=len(p['tasks']),timeout_seconds=r.get('fit_timeout_seconds'))));return 0
    return supervise(args,p,r)

if __name__=='__main__':raise SystemExit(main())
