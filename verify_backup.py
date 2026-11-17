"""Read-only byte verification. No R, shell, cloud API or network execution."""
import hashlib,json
from pathlib import Path
ROOT=Path(__file__).resolve().parent
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
m=json.loads((ROOT/'BACKUP_MANIFEST.json').read_text(encoding='utf-8'))
for row in m['files']:
    p=ROOT/row['path']
    assert p.is_file() and not p.is_symlink() and '..' not in Path(row['path']).parts,row['path']
    assert p.stat().st_size==row['bytes'] and sha(p)==row['sha256'],row['path']
base=ROOT/'snapshot'
for line in (base/'CODE_SHA256SUMS.txt').read_text(encoding='utf-8').splitlines():
    assert line[64:66]=='  ' and sha(base/line[66:])==line[:64],line[66:]
package=json.loads((base/'PACKAGE_MANIFEST.json').read_text(encoding='utf-8'))
assert sha(base/'CODE_SHA256SUMS.txt')==package['code_manifest_sha256']
for row in package['files']:assert sha(base/row['relative_path'])==row['sha256']
for folder in ['ops_v3','probe_call_v2']:
    a=base/'amendments'/folder
    for row in json.loads((a/'AMENDMENT.json').read_text(encoding='utf-8'))['files']:
        assert sha(a/row['relative_path'])==row['sha256']
cfg=json.loads((base/'amendments/ops_v3/EXECUTION_CONFIG.json').read_text(encoding='utf-8'))
assert cfg['new_fit_count']==120 and cfg['arms']==['original']
assert cfg['fitted_factor_capacities']==[4,4,4] and cfg['fitted_fpca_capacity']==4
assert cfg['annealing_sweeps']==99 and cfg['ordinary_t1_sweeps']==400
assert cfg['factor_shrinkage_by_arm']['original']['mode']=='none'
assert cfg['original_config_sha256']==sha(base/'config/EXECUTION_CONFIG.json')
print(json.dumps({'status':'PASS','backup_files_verified':len(m['files']),'frozen_source_manifest_sha256':sha(base/'CODE_SHA256SUMS.txt'),'numerical_fits':0,'cloud_calls':0}))
