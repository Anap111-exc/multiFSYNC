#!/usr/bin/env python3
"""Read-only process inventory; operator must review and explicitly attest it."""
import json, os, pathlib, subprocess, sys, time

def main():
    if not sys.platform.startswith('linux'): raise RuntimeError('Linux process audit only')
    raw=subprocess.check_output(['ps','-eo','pid=,ppid=,user=,comm=,stat=,etimes=,pcpu=,pmem='],
                                universal_newlines=True,timeout=10)
    rows=[]
    for line in raw.splitlines():
        a=line.split()
        if len(a)!=8: raise RuntimeError('Unexpected ps columns')
        pid=int(a[0]); started=None;pgrp=None;exe=None
        try:
            fields=pathlib.Path('/proc/{}/stat'.format(pid)).read_text().rsplit(')',1)[1].split()
            started=int(fields[19]);pgrp=int(fields[2])
            try:exe=os.readlink('/proc/{}/exe'.format(pid))
            except OSError:pass
        except (OSError,ValueError,IndexError): pass
        rows.append(dict(pid=pid,ppid=int(a[1]),user=a[2],command=a[3],state=a[4],
                         elapsed_seconds=int(a[5]),cpu_percent=float(a[6]),memory_percent=float(a[7]),start_ticks=started,
                         pgrp=pgrp,comm=a[3],exe=exe))
    print(json.dumps(dict(experiment_id='cavi_true222_fit444_m4_20261117_v1',
                         instance_id='i-bp1egmsairexm9w1rnks',region='cn-hangzhou',
                         observed_unix=time.time(),utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),
                         collector_pid=os.getpid(),processes=rows,full_ps_without_arguments=raw,
                         all_processes_operator_reviewed=False,no_other_tasks_require_instance=False,
                         all_experiment_tasks_complete_or_safe_exit=False,
                         approved_system_runtime_processes=[],mutations_performed=False),indent=2))
    return 0

if __name__=='__main__': sys.exit(main())
