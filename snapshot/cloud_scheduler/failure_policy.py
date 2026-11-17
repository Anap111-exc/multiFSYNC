"""Pre-registered narrow numerical-failure classification; never a retry."""
import re

POLICY=dict(id='known_guard_status_csv_v1',continue_remaining_registered=True,retry=False,replacement_seed=False)
PATTERNS=(
    r'[^\r\n]+: precision is not numerically symmetric\.',
    r'[^\r\n]+: Cholesky factorisation failed after adaptive jitter up to [0-9.eE+\-]+ times the diagonal scale\.',
    r'Expected RSS is non-finite at s=\d+, i=\d+, j=\d+\.',
)

def classify_failure(*,exit_code,fit_id,status_rows,bindings_ok,monitoring_errors,worker_error,resource_errors=()):
    if exit_code!=2 or bindings_ok is not True or monitoring_errors or worker_error or resource_errors:
        return 'FAILED'
    if len(status_rows)!=1:return 'FAILED'
    row=status_rows[0]
    if row.get('fit_id')!=fit_id or row.get('status')!='failed':return 'FAILED'
    error=row.get('error','')
    if not isinstance(error,str) or not any(re.fullmatch(pattern,error) for pattern in PATTERNS):return 'FAILED'
    return 'NUMERICAL_FAILED_NO_RETRY'
