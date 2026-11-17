"""Pure, deterministic resource decisions; independent of process execution."""
GiB=1024**3

def memory_budget(peak_bytes):
    return max(GiB,2*max(0,int(peak_bytes)))

def can_submit(*,available_bytes,owned_rss_bytes,running_budgets,candidate_peak_bytes,
               pilot_ready,running_count,max_concurrency=8,reserve_bytes=3*GiB):
    limit=max_concurrency if pilot_ready else 2
    if running_count>=limit:return False
    capacity=max(0,int(available_bytes))+max(0,int(owned_rss_bytes))
    required=reserve_bytes+sum(running_budgets)+memory_budget(candidate_peak_bytes)
    return capacity>=required

def pilot_observed(task_states,pilot_ids,minimum_seconds=300):
    return all(i in task_states and task_states[i].get('elapsed_seconds',0)>=minimum_seconds
               and task_states[i].get('sample_count',0)>=2 for i in pilot_ids)

def pilot_witnesses(task_states,tasks,minimum_seconds=300):
    witnesses=[]
    for task in tasks[:2]:
        row=task_states.get(task['id'],{})
        if (row.get('elapsed_seconds',0)>=minimum_seconds
                and row.get('sample_count',0)>=2
                and not row.get('monitoring_errors')):
            witnesses.append(task['id'])
    return witnesses


def pilot_ready(task_states,tasks,minimum_seconds=300):
    return len(pilot_witnesses(task_states,tasks,minimum_seconds))==2

def candidate_peak(task,states,tasks):
    by_id={t['id']:t for t in tasks}
    peaks=[s.get('peak_rss_bytes',0) for i,s in states.items() if i in by_id and by_id[i]['arm']==task['arm']]
    return max(peaks,default=0)
