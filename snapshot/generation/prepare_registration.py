from pathlib import Path
import csv,hashlib,json

ROOT=Path(__file__).resolve().parents[1]
EXPERIMENT=ROOT.name
OLD=Path('E:/multiFSYNC-results/cavi_mgp_overspec_l3_m4_20261105_v1')
NC=Path('E:/multiFSYNC-results/mgp_factor_count_pilot_m3_20261101_v3/package_src/generation')
GEN=ROOT/'generation';CFG=ROOT/'config'
CFG.mkdir(parents=True,exist_ok=True)
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def put(p,data):
    assert not p.exists(),f'Refuse overwrite {p}'
    p.write_bytes(data)
def csv_new(p,rows):
    assert not p.exists(),f'Refuse overwrite {p}'
    with p.open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
oldseeds=list(csv.DictReader((OLD/'config/SEED_MANIFEST.csv').open(encoding='utf-8-sig')))
oldplan=list(csv.DictReader((OLD/'config/PLAN.csv').open(encoding='utf-8-sig')))
used={int(r['seed']) for r in oldseeds}|{202611161,202611162}
rows=[];fitrows=[];seedrows=[]
roles=['data_seed','truth_basis_seed','truth_loading_seed','mean_truth_seed','test_time_seed','test_noise_seed']
for rep in range(1,11):
    ident=f'NC01_TRUE222_20261117__rep{rep:02d}'
    row=dict(data_id=ident,scenario_id='NC01_TRUE222',scenario_index=1,rep_id=rep,S=2,n_each=60,p=500,
             L_f=2,L_s_by_study='[2,2]',M=2,K=9,d_generated=0,d_fit=0,cross_scenario_base_realization_paired=False)
    for role in roles:
        counter=0
        while True:
            label=f'{EXPERIMENT}|{rep:02d}|{role}|{counter}'
            seed=int.from_bytes(hashlib.sha256(label.encode('ascii')).digest()[:8],'big')%2147483646+1
            if seed not in used:break
            counter+=1
        used.add(seed);row[role]=seed
        seedrows.append(dict(data_id=ident,seed_role=role,index=0,seed=seed,collision_counter=counter,status='REGISTERED_NOT_GENERATED'))
    previous=f'NC01_FACTOR_COUNT_20261101__rep{rep:02d}'
    selected=sorted((r for r in oldplan if r['dataset_id']==previous and r['arm']=='original'),key=lambda r:int(r['start']))
    assert len(selected)==12 and [int(r['start']) for r in selected]==list(range(1,13))
    for r in selected:
        fitrows.append(dict(dataset_id=ident,start=int(r['start']),seed=int(r['seed']),reference_dataset_id=previous,
                            reference_experiment=OLD.name,reference_arm='original'))
        seedrows.append(dict(data_id=ident,seed_role='fit_seed',index=int(r['start']),seed=int(r['seed']),collision_counter=0,status='REUSED_REGISTERED_START_SEED'))
    rows.append(row)
csv_new(CFG/'DATASET_REGISTRATION.csv',rows)
csv_new(GEN/'FIT_SEED_REFERENCE.csv',fitrows)
csv_new(GEN/'SEED_REGISTRATION.csv',seedrows)
source_files={}
for src,name in [(NC/'nc_common.R','nc_common.R'),
                 (NC/'stage7m_reconstruct.R','stage7m_reconstruct.R'),
                 (OLD/'scripts/project_fit_inputs.R','project_fit_inputs.R')]:
    put(GEN/name,src.read_bytes());source_files[name]={'source':str(src),'source_sha256':digest(src),'copied_sha256':digest(GEN/name)}
original=(NC/'nc_generation.R').read_text(encoding='utf-8')
core=original.split('\nnc_generate_scenario <- function',1)[0]
core=core.replace("assert(length(ls) == S && all(ls == 1L), 'Unregistered study-specific dimension')",
                  "assert(S == 2L && lf == 2L && length(ls) == S && all(ls == 2L), 'Unregistered true factor dimension')")
put(GEN/'nc_generation_core.R',core.encode('utf-8'))
runtime=(NC/'nc_runtime.R').read_text(encoding='utf-8')
canonical=runtime[runtime.index('nc_canonical_function_sha256 <- function'):runtime.index('\nnc_load_runtime <- function')]
fixed=runtime[runtime.index('nc_fixed_metadata <- function'):runtime.index('\nnc_metadata_checks <- function')]
put(GEN/'nc_fixed_metadata.R',(canonical+'\n'+fixed).encode('utf-8'))
source_files['nc_generation_core.R']={'source':str(NC/'nc_generation.R'),'source_sha256':digest(NC/'nc_generation.R'),
 'change':'Keep four generation/QC functions; only dimension guard changes to S2,Lf2,Ls2. Calibration and random generator calls unchanged.'}
source_files['nc_fixed_metadata.R']={'source':str(NC/'nc_runtime.R'),'source_sha256':digest(NC/'nc_runtime.R'),
 'change':'Extract canonical SHA and fixed-metadata functions without body modification.'}
contract=dict(schema='CAVI_TRUE222_GENERATION_V1',experiment_id=EXPERIMENT,status='REGISTERED_BEFORE_GENERATION',
 independent_datasets=10,registered_fits=120,starts_each=12,generated_on='Linux_only',
 baseline_library='/root/stage7i-exec-20260929-v1/stage7i_lib_ecs_20260929_v1',
 generator_canonical_sha256='24998da4ca8d80dfb32e2687183b976d599e4d86dfacb921f5165d57e42c9ec7',
 reconstruct_sha256='51d83c1864dbc346bb4fee0e25feee7e099698080e2cfc38199970f7fd92ea27',
 registration_sha256=digest(CFG/'DATASET_REGISTRATION.csv'),fit_seed_reference_sha256=digest(GEN/'FIT_SEED_REFERENCE.csv'),
 seed_registration_sha256=digest(GEN/'SEED_REGISTRATION.csv'),source_bindings=source_files,
 design=dict(S=2,n_s=[60,60],p=500,d=0,true_factor_counts=[2,2,2],true_fpca_rank=2,fitted_factor_capacities=[4,4,4],fitted_fpca_capacity=4,
 nonzero_loadings_each_column=50,loading_squared_norm_each_column=45,within_study_disjoint_support=True,
 truth_basis_dimension=8,truth_basis='iid_common_splinet_stiefel',fpca_weights=[.8,.2],noise_sd=.3,training_visits=[6,9],new_times_each_subject=3),
 seeds=dict(population_seeds='SHA256 of experiment_id|rep|role|collision_counter; first8bytes big-endian modulo2147483646 plus1; reject collisions with previous experiment and current registration.',
 fit_seeds='Exactly reuse the previous original-CAVI 12-start seed vector by rep_id, without reference to outcomes.',
 previous_seed_manifest_sha256=digest(OLD/'config/SEED_MANIFEST.csv'),previous_plan_sha256=digest(OLD/'config/PLAN.csv')),
 information_boundary=dict(fit='Observation-only projected files; no true L/M, coefficients, targets or responses.',
 evaluation='Separate sealed truth and test truth; read after winner selection is frozen.',
 generation_qc='Legality and source checks only; never regenerate inconvenient data or change seeds.'),
 automatic_retry=False,source_model_changes=False)
put(GEN/'GENERATION_CONTRACT.json',(json.dumps(contract,indent=2)+'\n').encode('utf-8'))
print(json.dumps({'registered_datasets':len(rows),'registered_start_seeds':len(fitrows),'first_data_id':rows[0]['data_id'],
 'last_data_id':rows[-1]['data_id'],'registration_sha256':contract['registration_sha256']},indent=2))
