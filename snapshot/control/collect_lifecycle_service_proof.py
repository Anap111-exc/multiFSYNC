"""Read-only Linux preflight: fingerprint only timedatex; no ECS API or secrets."""
import json
from pathlib import Path
import sys

EX=Path('/root/cavi-true222-fit444-m4-exec-20261117-v1')
sys.path.insert(0,str(EX))
from lifecycle import system_service
try:
    profile=system_service.enroll()
    result=dict(status='PASS',approved_restarting_system_services=[profile],network_calls=0)
except Exception as error:
    result=dict(status='UNAVAILABLE_KEEP_STRICT_PID_CHECK',approved_restarting_system_services=[],
                error_type=type(error).__name__,network_calls=0)
print(json.dumps(result,ensure_ascii=False))
