#!/usr/bin/env python3
"""Registered, bounded batch. Fits and selections never deserialize truth."""
import concurrent.futures as futures
import csv
import datetime as dt
import hashlib
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
import tarfile
import threading
import time

BASE = Path(__file__).resolve().parent
RESULT = Path('/root/multifsync-results/stage7l_prospective_20261001_v1')
LIB = '/root/stage7i-exec-20260929-v1/stage7i_lib_ecs_20260929_v1'
SUPPORT = '/root/multiFSYNC-covariate-routeb-pilot-20260927-v1'
ARCHIVE = Path('/root/stage7i-exec-20260929-v1/multiFSYNC-stage7i-20260929-v1.tar.gz')
MAX_WORKERS, FIT_SECONDS, BATCH_SECONDS = 8, 10800, 64800
CRITICAL = threading.Event()
ENV = dict(os.environ, **{k:'1' for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS',
    'MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','VECLIB_MAXIMUM_THREADS','BLIS_NUM_THREADS')})
ENV['LC_ALL'] = 'C.UTF-8'
TERMINAL_COLUMNS = ['fit_id','data_id','branch','start_index','fit_seed','source_U_fit_id',
 'status','error','elapsed_seconds','objective_eligible','invalid_reasons',
 'final_ordinary_elbo','global_iteration','ordinary_t1_sweeps','actual_updates',
 'strict_practical_converged','fit_path','snapshot50_path']

def stamp():
    return dt.datetime.now(dt.timezone.utc).isoformat()

def sha(path):
    h = hashlib.sha256()
    with open(path,'rb') as stream:
        for b in iter(lambda:stream.read(1024*1024),b''): h.update(b)
    return h.hexdigest()

def read_csv(path):
    with open(path,encoding='utf-8-sig',newline='') as stream:
        return list(csv.DictReader(stream))

def write_csv(path,rows,columns=None):
    assert not path.exists(), 'Refusing overwrite: '+str(path)
    path.parent.mkdir(parents=True,exist_ok=True)
    with open(path,'x',encoding='utf-8',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=columns or list(rows[0]),extrasaction='raise')
        writer.writeheader(); writer.writerows(rows)

def write_text(path,value):
    path.parent.mkdir(parents=True,exist_ok=True)
    with open(path,'x',encoding='utf-8') as stream: stream.write(value+'\n')

def register():
    data=[]; fits=[]
    for i in range(1,6):
        logical='covrb7l_%02d'%i; data_id=logical+'__C1'
        row=dict(logical_data_id=logical,data_id=data_id,data_index=i,beta_arm='C1',
          observation_source_data_id=data_id,data_seed=101011000+i,
          truth_basis_seed=101012000+i,truth_loading_seed=101013000+i,
          mean_truth_seed=101014000+i,beta_bspline_seed=101015000+i,
          d_generated=2,d_fit=2,truth_available_to_fit='FALSE',
          truth_available_to_selection='FALSE',truth_unseal_authorized='FALSE')
        data.append(row)
        for j in range(1,11):
            fits.append(dict(fit_id=data_id+'__start%02d'%j,data_id=data_id,
              logical_data_id=logical,beta_arm='C1',start_index=j,
              fit_seed=101016000+i*100+j,n_cpus=1,
              truth_available_to_fit='FALSE',truth_available_to_selection='FALSE'))
    seeds=[r[k] for r in data for k in ('data_seed','truth_basis_seed',
      'truth_loading_seed','mean_truth_seed','beta_bspline_seed')]+[r['fit_seed'] for r in fits]
    assert len(seeds)==75 and len(set(seeds))==75
    write_csv(BASE/'DATA_MANIFEST.csv',data)
    write_csv(BASE/'FIT_MANIFEST.csv',fits)
    freeze={'registered_at':stamp(),'data_manifest_sha256':sha(BASE/'DATA_MANIFEST.csv'),
      'fit_manifest_sha256':sha(BASE/'FIT_MANIFEST.csv'),
      'protocol_sha256':sha(BASE/'PROTOCOL_STAGE7L_PROSPECTIVE_20261001_v1.md'),
      'datasets':5,'U':50,'S_max':5,'seed_collision_check':'PASS_exact_integer_boundaries',
      'generation_started':False,'truth_opened':False}
    write_text(BASE/'REGISTRATION_FROZEN.txt',json.dumps(freeze,indent=2))

def preflight():
    reg=json.loads((BASE/'REGISTRATION_FROZEN.txt').read_text())
    for filename,key in [('DATA_MANIFEST.csv','data_manifest_sha256'),
       ('FIT_MANIFEST.csv','fit_manifest_sha256'),
       ('PROTOCOL_STAGE7L_PROSPECTIVE_20261001_v1.md','protocol_sha256')]:
        assert sha(BASE/filename)==reg[key], 'Registration binding mismatch'
    for line in (BASE/'SOURCE_FILES.sha256').read_text().splitlines():
        expected,filename=line.split('  ',1)
        assert sha(BASE/filename)==expected,'Source binding mismatch: '+filename
    expected={ARCHIVE:'b6c2210089061d99e31825e358d9c03b6d867bcd7fdde21e8432e784e4436178',
      Path(LIB)/'multiFSYNC/DESCRIPTION':'80780d8b34af04fc5519552b96e332342e92938dff33315efb9239805f3953dd',
      Path(LIB)/'multiFSYNC/NAMESPACE':'fd33c1d9e9ae138e3ee1e5ee80dc67c3a9e261d576af09cd8ae19df5d3033ee3',
      Path(LIB)/'multiFSYNC/R/multiFSYNC.rdb':'073cbdb6b53667fe1a1cc5b996cb8b50b52b4e39565dff74a96083bcc9d387a8',
      Path(LIB)/'multiFSYNC/R/multiFSYNC.rdx':'42a6c949608daebdf5da08e2bffaafe0398c12505e927e5c1897a07ef0c34194'}
    for path,value in expected.items(): assert sha(path)==value,'Package binding conflict: '+str(path)
    with tarfile.open(ARCHIVE,'r:gz') as archive:
        members=[m for m in archive.getmembers() if m.name.endswith('/R/generate_data_structured.R')]
        assert len(members)==1
        assert hashlib.sha256(archive.extractfile(members[0]).read()).hexdigest()==\
          '226131949adf52a28fc5e8e91b765a8cdb98a4e2c3a2ffbfa7ee40fc7bc932dd'
    parity=Path('/root/multifsync-results/stage7j_seed02_u_20260930_v1/u_parity_v1/PARITY_PASS.txt')
    assert parity.read_text().splitlines()[0]=='STAGE7J_U_PARITY_PASS'
    for name,count in [('ALL_STATE_PARITY.csv',19307),('ORDINARY_ELBO_PARITY.csv',400)]:
        rows=read_csv(parity.parent/name)
        assert len(rows)==count and all(r['pass']=='TRUE' for r in rows),'Historical parity evidence failed'
        expected[parity.parent/name]=sha(parity.parent/name)
    expected[parity]=sha(parity)
    write_csv(BASE/'PREFLIGHT_EXTERNAL_BINDING.csv',[
        dict(role=str(p.name),path=str(p),sha256=h) for p,h in expected.items()])
    subprocess.run(['Rscript',str(BASE/'preflight_stage7l.R'),str(BASE),LIB,SUPPORT],
                   env=ENV,check=True)
    assert os.statvfs('/root').f_bavail*os.statvfs('/root').f_frsize>3*1024**3,'Insufficient ECS disk'
    write_text(BASE/'PREFLIGHT_PASS.txt','STAGE7L_PREFLIGHT_PASS\n'+stamp())

def r_command(action,*args):
    return ['Rscript',str(BASE/'run_stage7l.R'),action,str(BASE),str(RESULT),LIB]+list(args)

def failure(row,branch,status,error,elapsed=0):
    fit_id=row['fit_id']; folder=RESULT/'fits'/fit_id
    folder.mkdir(parents=True,exist_ok=True)
    if (folder/'terminal.csv').exists(): return
    state=folder/'STATE_SWEEP050.rds'
    saved_iterations=[int(p.stem.replace('STATE_SWEEP','')) for p in folder.glob('STATE_SWEEP*.rds')]
    # Only a lower bound is known after a killed process. Never invent a final iteration.
    reached=max(saved_iterations,default=50 if branch=='S' else 0)
    terminal=dict(zip(TERMINAL_COLUMNS,['']*len(TERMINAL_COLUMNS)))
    terminal.update(fit_id=fit_id,data_id=row['data_id'],branch=branch,
      start_index=row.get('start_index',0),fit_seed=row.get('fit_seed',''),
      source_U_fit_id=row.get('source_U_fit_id',''),status=status,error=error,
      elapsed_seconds=elapsed,objective_eligible='FALSE',invalid_reasons=status,
      actual_updates=max(0,reached-(50 if branch=='S' else 0)),
      strict_practical_converged='FALSE',snapshot50_path=str(state) if state.exists() else '')
    write_csv(folder/'terminal.csv',[terminal],TERMINAL_COLUMNS)
    write_text(folder/'TECHNICAL_FAILURE.txt',status+'\n'+error+'\nupdates_are_recorded_lower_bound=TRUE')

def bounded_task(row,branch,deadline):
    fit_id=row['fit_id']; begin=time.monotonic()
    if CRITICAL.is_set():
        failure(row,branch,'DEPENDENCY_BLOCKED','Critical input/source/interface failure'); return
    if begin>=deadline:
        failure(row,branch,'BATCH_CAP_NOT_STARTED','Batch fitting wall cap reached'); return
    seconds=min(FIT_SECONDS,deadline-begin)
    log=RESULT/'logs'/(fit_id+'.log'); log.parent.mkdir(exist_ok=True)
    with open(log,'x') as stream:
        process=subprocess.Popen(r_command('fit',branch,fit_id),env=ENV,
                                 stdout=stream,stderr=subprocess.STDOUT,start_new_session=True)
        limit=begin+seconds
        while True:
            remaining=limit-time.monotonic()
            if remaining<=0 or CRITICAL.is_set():
                break
            try:
                code=process.wait(timeout=min(15,remaining)); break
            except subprocess.TimeoutExpired:
                pass
        if process.poll() is None:
            os.killpg(process.pid,signal.SIGTERM)
            try: process.wait(timeout=30)
            except subprocess.TimeoutExpired:
                os.killpg(process.pid,signal.SIGKILL); process.wait()
            reason='DEPENDENCY_BLOCKED' if CRITICAL.is_set() else (
                'TIMEOUT_BATCH' if seconds<FIT_SECONDS else 'TIMEOUT_FIT')
            failure(row,branch,reason,'Technical stop; no retry',time.monotonic()-begin)
            return
    text=log.read_text(errors='replace')
    if any(phrase in text for phrase in ('Binding conflict:', 'input binding missing',
        'metadata disagrees with fit static fields', 'Required frozen hooks unavailable')):
        CRITICAL.set()
    if not (RESULT/'fits'/fit_id/'terminal.csv').exists():
        failure(row,branch,'PROCESS_ERROR','R exit code '+str(code),time.monotonic()-begin)

def eligible(row):
    import math
    try: return row['status']=='COMPLETE' and row['objective_eligible']=='TRUE' and math.isfinite(float(row['final_ordinary_elbo']))
    except (ValueError,TypeError): return False

def freeze():
    original=read_csv(BASE/'FIT_MANIFEST.csv')
    data=read_csv(BASE/'DATA_MANIFEST.csv'); endpoints=[]; paired=[]
    u_winners=read_csv(RESULT/'U_WINNERS_FROZEN.csv')
    for d in data:
        rows=[read_csv(RESULT/'fits'/r['fit_id']/'terminal.csv')[0] for r in original if r['data_id']==d['data_id']]
        s=read_csv(RESULT/'fits'/(d['data_id']+'__S')/'terminal.csv')[0]
        endpoints.extend(rows+[s]); pool=[r for r in rows+[s] if eligible(r)]
        chosen=sorted(pool,key=lambda r:(-float(r['final_ordinary_elbo']),r['fit_id']))[0] if pool else None
        u=[r for r in u_winners if r['data_id']==d['data_id']][0]
        paired.append(dict(data_id=d['data_id'],U_winner_fit_id=u['fit_id'],S_fit_id=s['fit_id'],
          enhanced_winner_fit_id=chosen['fit_id'] if chosen else '',
          selected_S='TRUE' if chosen and chosen['branch']=='S' else 'FALSE'))
    assert len(endpoints)==55 and len(paired)==5
    freeze_dir=RESULT/'blind_freeze'; freeze_dir.mkdir()
    write_csv(freeze_dir/'ALL_ENDPOINTS_BLIND.csv',endpoints,TERMINAL_COLUMNS)
    write_csv(freeze_dir/'PAIRED_SELECTION_BLIND.csv',paired)
    failed=[r for r in endpoints if r['status']!='COMPLETE']
    write_csv(freeze_dir/'FAILURE_RECORDS.csv',failed,TERMINAL_COLUMNS)
    files=sorted(p for p in RESULT.rglob('*') if p.is_file())
    files+=sorted(p for p in BASE.rglob('*') if p.is_file() and not p.name.endswith('.log') and '__pycache__' not in p.parts)
    write_csv(freeze_dir/'BLIND_FILE_HASHES.csv',[dict(path=str(p),sha256=sha(p),bytes=p.stat().st_size) for p in files])
    markers={n:sha(freeze_dir/n) for n in ['ALL_ENDPOINTS_BLIND.csv','PAIRED_SELECTION_BLIND.csv','FAILURE_RECORDS.csv','BLIND_FILE_HASHES.csv']}
    write_text(freeze_dir/'BLIND_ENDPOINTS_FROZEN.txt',json.dumps(dict(
      status='STAGE7L_BLIND_ENDPOINTS_FROZEN',frozen_at=stamp(),truth_opened=False,
      endpoints=55,selection_uses_truth=False,hashes=markers),indent=2))
    os.sync()

def run():
    assert (BASE/'PREFLIGHT_PASS.txt').exists(),'Preflight not completed'
    assert not RESULT.exists(),'Refusing experiment restart'
    subprocess.run(r_command('generate'),env=ENV,check=True)
    deadline=time.monotonic()+BATCH_SECONDS; started=stamp()
    original=read_csv(BASE/'FIT_MANIFEST.csv')
    with futures.ThreadPoolExecutor(max_workers=MAX_WORKERS) as pool:
        list(pool.map(lambda row:bounded_task(row,'U',deadline),original))
    assert not CRITICAL.is_set(),'Critical evidence/interface conflict; dependent work stopped'
    subprocess.run(r_command('select_u'),env=ENV,check=True)
    winners=read_csv(RESULT/'U_WINNERS_FROZEN.csv'); s_rows=[]
    for u in winners:
        s=dict(u,fit_id=u['S_fit_id'],source_U_fit_id=u['fit_id'])
        if not u['fit_id']:
            failure(s,'S','NO_ELIGIBLE_U','No frozen eligible U winner'); continue
        if time.monotonic()>=deadline:
            failure(s,'S','BATCH_CAP_NOT_STARTED','Wall cap reached before S construction'); continue
        construction_log=RESULT/'logs'/(s['fit_id']+'_construction.log')
        with open(construction_log,'x') as log:
            try:
                completed=subprocess.run(r_command('construct_s',u['data_id']),env=ENV,
                   stdout=log,stderr=subprocess.STDOUT,timeout=min(600,deadline-time.monotonic()))
                code=completed.returncode
            except subprocess.TimeoutExpired: code=124
        if code:
            failure(s,'S','CONSTRUCTION_FAILED','Fixed helper failed; see construction evidence'); continue
        s_rows.append(s)
    with futures.ThreadPoolExecutor(max_workers=MAX_WORKERS) as pool:
        list(pool.map(lambda row:bounded_task(row,'S',deadline),s_rows))
    assert not CRITICAL.is_set(),'Critical evidence/interface conflict; dependent work stopped'
    write_text(RESULT/'FITTING_COST.txt',json.dumps(dict(started_at=started,
      finished_at=stamp(),wall_seconds=BATCH_SECONDS-(deadline-time.monotonic()),
      U_planned=50,S_planned=5,S_constructed=len(s_rows),max_workers=8),indent=2))
    freeze()
    subprocess.run(['Rscript',str(BASE/'evaluate_stage7l.R'),str(BASE),str(RESULT),LIB,SUPPORT],env=ENV,check=True)
    write_text(RESULT/'EXPERIMENT_COMPLETE.txt','STAGE7L_EXPERIMENT_COMPLETE\n'+stamp())
    final_files=sorted(p for p in RESULT.rglob('*') if p.is_file())
    write_csv(RESULT/'COMPLETE_FILE_HASHES.csv',[dict(path=str(p),sha256=sha(p),bytes=p.stat().st_size) for p in final_files])
    small=[p for p in RESULT.rglob('*') if p.is_file() and p.suffix in ('.csv','.txt','.md','.json')]
    with tarfile.open(RESULT/'stage7l_small_readout.tar.gz','x:gz') as bundle:
        for p in sorted(small): bundle.add(p,arcname=str(p.relative_to(RESULT)))
    os.sync()

if __name__=='__main__':
    {'register':register,'preflight':preflight,'run':run}[sys.argv[1]]()
