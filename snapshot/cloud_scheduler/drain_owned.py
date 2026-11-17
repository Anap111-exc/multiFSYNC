"""Wait only for authenticated registered worker groups after supervisor failure.

No signal, retry, deadline, process-name inference, or blind PID wait is used.
An unprovable live group is reported as unresolved and must block automatic stop.
"""
from pathlib import Path
import hashlib,json,os,time

def read(path):return json.loads(Path(path).read_text(encoding='utf-8-sig'))
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def persist(path,value):
    path=Path(path);tmp=path.with_name(path.name+'.part')
    with tmp.open('w',encoding='utf-8') as stream:
        json.dump(value,stream,indent=2);stream.write('\n');stream.flush();os.fsync(stream.fileno())
    os.replace(str(tmp),str(path))

def inventory():
    rows=[]
    for path in Path('/proc').iterdir():
        if not path.name.isdigit():continue
        try:
            raw=(path/'stat').read_text();fields=raw[raw.rfind(')')+2:].split()
            rows.append(dict(pid=int(path.name),ppid=int(fields[1]),pgrp=int(fields[2]),
                             start_ticks=int(fields[19]),state=fields[0]))
        except (FileNotFoundError,ProcessLookupError):continue
    return rows

def identity(row):return {k:int(row[k]) for k in ('pid','start_ticks','pgrp')}
def matches(row,expected):return bool(row and all(row[k]==expected[k] for k in ('pid','start_ticks','pgrp')))
def live(row):return not row['state'].startswith('Z')

def discover(output_root,registered_ids):
    out=Path(output_root);scheduler=out/'scheduler';groups=[];unresolved=[]
    planfile=out/'fit_plan_v1/PLAN.json';releasefile=out/'fit_plan_v1/RELEASE.json'
    if not scheduler.exists():return groups,unresolved
    if not planfile.is_file() or not releasefile.is_file():
        return [],[dict(reason='scheduler_exists_without_bound_plan_release')]
    psha,rsha=sha(planfile),sha(releasefile)
    plan=read(planfile);tasks={t['id']:t for t in plan['tasks']}
    if set(tasks)!=set(registered_ids):return [],[dict(reason='drain_plan_registration_mismatch')]
    ledger={}
    try:
        ledger=read(scheduler/'LEDGER.json')
        if ledger.get('plan_sha256')!=psha or ledger.get('release_sha256')!=rsha:
            raise RuntimeError('ledger binding mismatch')
    except Exception as error:
        unresolved.append(dict(reason='ledger_unreadable_or_unbound',error_type=type(error).__name__))
        ledger={}
    for task_id in registered_ids:
        folder=scheduler/task_id;record=None;item=ledger.get('tasks',{}).get(task_id,{})
        try:
            wpath=folder/'WORKER_IDENTITY.json'
            if wpath.is_file():
                w=read(wpath)
                if w.get('plan_sha256')!=psha or w.get('release_sha256')!=rsha:raise RuntimeError('worker binding mismatch')
                record=identity(w)
            elif item.get('identity'):record=identity(item['identity'])
            if record is None:
                if item.get('attempts',0)>0:unresolved.append(dict(task_id=task_id,reason='attempted_worker_has_no_bound_identity'))
                continue
            if record['pid']!=record['pgrp']:raise RuntimeError('worker was not a session leader')
            known=[record]
            fpath=folder/'FIT_IDENTITY.json'
            if fpath.is_file():
                fit=read(fpath);proc=fit.get('proc')
                if fit.get('command')!=tasks[task_id]['command']:raise RuntimeError('fit command binding mismatch')
                if proc:
                    member=identity(proc)
                    if member['pid']!=fit['pid'] or member['pgrp']!=record['pgrp']:raise RuntimeError('fit process group mismatch')
                    known.append(member)
            groups.append(dict(task_id=task_id,leader=record,known_members=known,state='UNOBSERVED'))
        except Exception as error:
            unresolved.append(dict(task_id=task_id,reason='invalid_owned_identity',error_type=type(error).__name__))
    return groups,unresolved

def drain(output_root,registered_ids,reason,inventory_fn=inventory,sleep_fn=time.sleep,poll_seconds=5):
    """Persist progress until known groups exit; never wait on a reused/unproven PID."""
    out=Path(output_root);log=out/'pipeline';log.mkdir(parents=True,exist_ok=True)
    groups,unresolved=discover(out,registered_ids)
    started=time.time();observations=0
    persist(log/'DRAIN_START.json',dict(reason=reason,started_unix=started,owned_groups=groups,
        unresolved=unresolved,timeout_seconds=None,signals_sent=0,retries=0))
    while True:
        rows=inventory_fn() if groups else [];observations+=1
        by_pid={r['pid']:r for r in rows};waiting=[]
        for group in groups:
            if group['state'] in ('EXITED','UNRESOLVED'):continue
            expected=group['leader'];members=[r for r in rows if r['pgrp']==expected['pgrp'] and live(r)]
            if not members:
                group['state']='EXITED';continue
            leader=by_pid.get(expected['pid'])
            if leader and live(leader) and not matches(leader,expected):
                group['state']='UNRESOLVED'
                unresolved.append(dict(task_id=group['task_id'],reason='leader_PID_reused',expected=expected,observed=leader))
                continue
            proof=any(matches(member,known) for member in members for known in group['known_members'])
            if not proof:
                group['state']='UNRESOLVED'
                unresolved.append(dict(task_id=group['task_id'],reason='live_group_without_bound_surviving_identity',members=members))
                continue
            # A verified survivor proves the original process group still exists;
            # remember its current members so leader exit does not orphan the drain.
            known={(m['pid'],m['start_ticks'],m['pgrp']) for m in group['known_members']}
            for member in members:
                key=(member['pid'],member['start_ticks'],member['pgrp'])
                if key not in known:group['known_members'].append(identity(member));known.add(key)
            group['state']='WAITING_OWNED_GROUP';waiting.append(group['task_id'])
        state=dict(status='WAITING' if waiting else ('UNRESOLVED' if unresolved else 'DRAINED'),
            reason=reason,started_unix=started,observed_unix=time.time(),observations=observations,
            waiting_task_ids=waiting,groups=groups,unresolved=unresolved,signals_sent=0,retries=0,
            timeout_seconds=None,owned_process_groups=[g['leader'] for g in groups])
        persist(log/'DRAIN_STATE.json',state)
        if not waiting:
            persist(log/'DRAIN_COMPLETE.json',state)
            return state
        sleep_fn(poll_seconds)
