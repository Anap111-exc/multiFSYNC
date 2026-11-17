#!/usr/bin/env python3
"""Summarize frozen original-CAVI science; no fitting, selection or truth access."""
import csv, hashlib, json, math, statistics, sys
from collections import defaultdict
from pathlib import Path

def read(path):
    with path.open(encoding='utf-8-sig',newline='') as stream:
        return list(csv.DictReader(stream))

def number(value):
    if str(value).lower() in ('true','false'):
        return float(str(value).lower()=='true')
    try:
        value=float(value)
        return value if math.isfinite(value) else None
    except (TypeError,ValueError):
        return None

def mean(values):
    values=[x for x in values if x is not None]
    return statistics.mean(values) if values else None

def sd(values):
    values=[x for x in values if x is not None]
    return statistics.stdev(values) if len(values)>1 else None

def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda:stream.read(1024*1024),b''):h.update(block)
    return h.hexdigest()

def save_csv(path, rows):
    fields=list(dict.fromkeys(key for row in rows for key in row))
    with path.open('x',encoding='utf-8',newline='') as stream:
        writer=csv.DictWriter(stream,fields);writer.writeheader();writer.writerows(rows)

def summarize(output):
    output=Path(output).resolve();plan=read(output/'PLAN.csv')
    ids=list(dict.fromkeys(r['dataset_id'] for r in plan))
    assert len(ids)==10 and len(plan)==120, 'Expected ten datasets, one arm and twelve starts'
    assert {row['arm'] for row in plan}=={'original'}
    for data_id in ids:
        arm_rows={arm:sorted([row for row in plan if row['dataset_id']==data_id and row['arm']==arm],
                            key=lambda row:int(row['start'])) for arm in ('original',)}
        assert all([int(row['start']) for row in rows]==list(range(1,13)) for rows in arm_rows.values())
    target=output/'aggregate'
    assert not target.exists(), 'Refuse aggregate overwrite or automatic retry'
    target.mkdir()
    inputs=[];pool=[];metrics=[];statuses=[];fpca_roles=[]
    def load(path):
        inputs.append(dict(path=str(path),sha256=sha(path),bytes=path.stat().st_size))
        return read(path)
    def add(data_id,arm,metric,value,source,scope='all10'):
        metrics.append(dict(data_id=data_id,arm=arm,metric=metric,value=number(value),scope=scope,source=source))
    def numeric_table(rows,data_id,source,keys,exclude=()):
        for row in rows:
            label='.'.join(str(row[k]) for k in keys)
            for key,value in row.items():
                if key in set(keys)|set(exclude)|{'data_id','dataset_id','arm','fit_id','start','ELBO'}:continue
                if number(value) is not None:
                    add(data_id,row['arm'],source+'.'+label+'.'+key,value,source)
    for data_id in ids:
        root=output/'science'/data_id
        marker=root/'SCIENCE_COMPLETE.json'
        assert marker.is_file(), 'Missing scientific completion: '+data_id
        status=json.loads(marker.read_text(encoding='utf-8'))
        assert status['status']=='SCIENCE_COMPLETE' and status['data_id']==data_id
        assert status['fit_updates']==0 and status['truth_used_for_selection'] is False
        assert status['selection_frozen_before_truth'] is True and status['registered']==12
        inputs.append(dict(path=str(marker),sha256=sha(marker),bytes=marker.stat().st_size));statuses.append(status)
        p=load(root/'POOL_SCIENCE.csv');assert len(p)==1 and {x['arm'] for x in p}=={'original'}
        pool.extend(p)
        rank_rows=load(root/'WINNER_FPCA_ROLE_RANKS.csv')
        for arm in ('original',):
            arm_pool=next(row for row in p if row['arm']==arm)
            arm_ranks=[row for row in rank_rows if row['arm']==arm]
            winner_exists=number(arm_pool['winner_exists'])==1
            assert len(arm_ranks)==(6 if winner_exists else 0)
            for block,study,role in (('shared','0','shared'),('specific','1','specific_1'),('specific','2','specific_2')):
                for true_factor in (1,2):
                    matched=[row for row in arm_ranks if row['block']==block and row['study']==study and int(row['true_factor'])==true_factor]
                    assert len(matched)==(1 if winner_exists else 0)
                    rank=number(matched[0]['matched_fpca_M_09']) if matched else None
                    assert rank is None or (rank.is_integer() and 0<=rank<=4)
                    category='missing_winner' if not winner_exists else ('missing_role' if rank is None else str(int(rank)))
                    fpca_roles.append(dict(data_id=data_id,arm=arm,role=role,true_factor=true_factor,truth_rank=2,
                        matched_rank=rank,category=category,winner_exists=winner_exists,exact_rank=rank==2,
                        matching='one_to_one_own_block_factor_PPI_ge_0.5_cosine_ge_0.9_max_cardinality_then_cosine'))
        for filename,keys in [('WINNER_RAW_COMPONENT_METRICS.csv',['scope','object']),
                              ('WINNER_Q_INTERVAL_METRICS.csv',['scope','object']),
                              ('WINNER_HISTORICAL_SUMMARY.csv',[]),
                              ('WINNER_SCIENCE.csv',[])]:
            numeric_table(load(root/filename),data_id,filename[:-4],keys)
        for arm in ('original',):
            folder=root/arm
            if not folder.exists():continue  # An unavailable winner remains in its registered pool denominator.
            for filename,keys in [('ACTIVE_INTERVALS.csv',['scope','object','mask']),
                                  ('CENTERED_SIGNAL.csv',['scope','object']),('NOISE.csv',['mask','scale'])]:
                rows=load(folder/filename)
                for row in rows:row['arm']=arm
                numeric_table(rows,data_id,filename[:-4],keys)
            comp=load(folder/'COMPONENT_METRICS.csv')
            for component in ('1','2','all'):
                rr=[r for r in comp if component=='all' or r['true_component']==component]
                selected=['stage7m_function_ise','stage7m_function_cosine','score_nrmse','correlation',
                          'stage7m_eigenvalue_error','stage7m_estimated_eigenvalue','stage7m_estimated_pve']
                for metric in selected:
                    groups=defaultdict(list)
                    for row in rr:
                        key=(row['true_role'],'0' if row['true_role']=='shared' else row['study'],row['true_factor'])
                        groups[key].append(number(row[metric]))
                    add(data_id,arm,'canonical_unique_factor.component_'+component+'.'+metric,
                        mean([mean(values) for values in groups.values()]),'COMPONENT_METRICS')
                groups=defaultdict(list)
                for row in rr:
                    key=(row['true_role'],'0' if row['true_role']=='shared' else row['study'],row['true_factor'])
                    value=number(row['stage7m_eigenvalue_error']);groups[key].append(abs(value) if value is not None else None)
                add(data_id,arm,'canonical_unique_factor.component_'+component+'.absolute_eigenvalue_error',
                    mean([mean(values) for values in groups.values()]),'COMPONENT_METRICS')
            loads=load(folder/'LOADING_METRICS.csv')
            for rep in sorted({r['representation'] for r in loads}):
                rr=[r for r in loads if r['representation']==rep]
                for metric in ('loading_cosine','loading_norm_ratio','loading_relative_l2','roc_auc','pr_average_precision',
                               'precision_ppi05','recall_ppi05','f1_ppi05'):
                    groups=defaultdict(list)
                    for row in rr:
                        key=(row['true_role'],'0' if row['true_role']=='shared' else row['study'],row['true_factor'])
                        groups[key].append(number(row[metric]))
                    add(data_id,arm,'loading.'+rep+'.'+metric,mean([mean(v) for v in groups.values()]),'LOADING_METRICS')
            blocks=load(folder/'BLOCK_SUBSPACE.csv')
            for scope in sorted({r['scope'] for r in blocks}):
                rr=[r for r in blocks if r['scope']==scope and r['metric_class']=='correct_block']
                values={('shared' if r['truth_block']=='shared' else 'specific_'+r['truth_study']):
                        number(r['truth_energy_trace_coverage']) for r in rr}
                for block,label in [('shared','MTA'),('specific_1','MTB1'),('specific_2','MTB2')]:
                    if block in values:add(data_id,arm,'loading.'+scope+'.'+label,values[block],'BLOCK_SUBSPACE')
                add(data_id,arm,'loading.'+scope+'.MTB',mean([v for k,v in values.items() if k.startswith('specific_')]),'BLOCK_SUBSPACE')
    assert len(statuses)==10 and len(pool)==10
    seen=set()
    for row in metrics:
        key=(row['data_id'],row['arm'],row['metric']);assert key not in seen, 'Duplicate metric '+str(key);seen.add(key)
    grouped=defaultdict(list)
    for row in metrics:grouped[row['metric']].append(row)
    scientific=[]
    for metric,rows in sorted(grouped.items()):
        values=[row['value'] for row in rows if row['value'] is not None]
        spread=sd(values)
        scientific.append(dict(metric=metric,arm='original',registered_datasets=10,finite_n=len(values),
            mean=mean(values),sd=spread,median=statistics.median(values) if values else None,
            mcse_mean=spread/math.sqrt(len(values)) if spread is not None else None))
    structural=[];counts=[]
    for arm in ('original',):
        rows=[r for r in pool if r['arm']==arm];assert len(rows)==10
        structural.append(dict(arm=arm,datasets=10,registered_starts=120,valid_starts=sum(int(r['valid']) for r in rows),
          failed_or_invalid=sum(int(r['failed_or_invalid']) for r in rows),correct_candidates=sum(int(r['correct_candidates']) for r in rows),
          pools_with_correct_structure=sum(number(r['pool_has_correct_structure']) for r in rows),
          winners_available=sum(number(r['winner_exists']) for r in rows),
          winners_correct=sum(number(r['winner_correct']) for r in rows),
          winner_hit_rate_registered=sum(number(r['winner_correct']) for r in rows)/10))
        for name in ('shared','specific_1','specific_2','total'):
            truth_count=6 if name=='total' else 2
            values=[]
            for row in rows:
                blocks=[number(row['winner_'+block]) for block in ('shared','specific_1','specific_2')]
                values.append(sum(blocks) if name=='total' and all(value is not None for value in blocks)
                              else number(row.get('winner_'+name)))
            finite=[value for value in values if value is not None]
            spread=sd(finite)
            counts.append(dict(arm=arm,role=name,truth=truth_count,registered_datasets=10,available=len(finite),missing=10-len(finite),
                mean_count=mean(finite),sd_count=spread,mcse_mean=spread/math.sqrt(len(finite)) if spread is not None else None,
                mean_bias=mean([value-truth_count for value in finite]),MAE=mean([abs(value-truth_count) for value in finite]),
                RMSE=math.sqrt(mean([(value-truth_count)**2 for value in finite])) if finite else None,
                underselected=sum(value<truth_count for value in finite),exact=sum(value==truth_count for value in finite),
                overselected=sum(value>truth_count for value in finite)))
    rank_frequencies=[];rank_dataset=[]
    for arm in ('original',):
        arm_roles=[row for row in fpca_roles if row['arm']==arm]
        assert len(arm_roles)==60
        for role in ('all_true_roles','shared','specific_1','specific_2'):
            records=[row for row in arm_roles if role=='all_true_roles' or row['role']==role]
            for category in ('0','1','2','3','4','missing_role','missing_winner'):
                count=sum(row['category']==category for row in records)
                rank_frequencies.append(dict(arm=arm,role=role,category=category,count=count,
                    registered_true_roles=len(records),fraction_registered=count/len(records)))
        for data_id in ids:
            records=[row for row in arm_roles if row['data_id']==data_id]
            rank_dataset.append(dict(data_id=data_id,arm=arm,registered_roles=6,
                matched_roles=sum(row['matched_rank'] is not None for row in records),
                exact_rank_roles=sum(row['exact_rank'] for row in records),
                all_six_true_roles_rank2=all(row['exact_rank'] for row in records)))
    save_csv(target/'CANDIDATE_POOLS.csv',pool);save_csv(target/'STRUCTURE_RECOVERY_SUMMARY.csv',structural)
    save_csv(target/'FACTOR_COUNT_SUMMARY.csv',counts);save_csv(target/'WINNER_METRICS_LONG.csv',metrics)
    save_csv(target/'SCIENTIFIC_METRIC_SUMMARY.csv',scientific)
    save_csv(target/'FPCA_ROLE_RANKS.csv',fpca_roles)
    save_csv(target/'FPCA_RANK_FREQUENCIES.csv',rank_frequencies)
    save_csv(target/'FPCA_RANK_DATASET_SUMMARY.csv',rank_dataset)
    save_csv(target/'INPUT_SCIENCE_HASHES.csv',inputs)
    lines=['# Original CAVI: true (2,2,2), fitted (4,4,4), FPCA 2 to 4',
      '', 'Ten independent generated datasets; two studies, 60 subjects per study, 500 variables; twelve original CAVI starts each. Each pool selects maximum complete ordinary ELBO before reading truth.',
      '', '| Method | Correct candidates / 120 | Pools with a correct candidate / 10 | Correct winner / 10 | Invalid starts |',
      '|---|---:|---:|---:|---:|']
    for r in structural:lines.append('| {arm} | {correct_candidates} | {pools_with_correct_structure:g} | {winners_correct:g} | {failed_or_invalid} |'.format(**r))
    lines+=['', 'Factor PPI >= 0.5 is the formal count rule. Count bias/MAE/RMSE and under/exact/over counts include explicit unavailable-winner totals.',
      'Each dataset has six distinct true global roles; missing roles and missing winners remain in FPCA denominators. Shared roles are not double counted across studies.',
      'Scientific means and data-level SD/MCSE use independent datasets, not starts. Full raw contributions include every fitted factor; no truth-selected pruning. Reserved-time prediction refers to existing subjects, not independent new subjects.',
      'MTA/MTB use the frozen energy-weighted loading-subspace trace definitions. Function/score alignment is for evaluation only. Intervals are pointwise variational moment intervals.',
      'Source hashes, finite denominators, and separate pool reachability versus winner accuracy are retained. No MGP arm or artificial paired contrast.']
    (target/'RESULT_REPORT.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    marker=dict(status='SCIENTIFIC_EVALUATION_COMPLETE',datasets=10,registered_full_fits=120,
      new_original_fits_registered=120,new_MGP_fits_registered=0,selection_frozen_before_truth=True,
      fit_updates_by_science=0,scientific_metric_count=len(scientific),truth_used_for_selection=False,
      incomplete_original_pools=sum(int(r['valid'])<12 for r in pool),
      total_evaluation_seconds=sum(s['seconds'] for s in statuses))
    (target/'SUMMARY.json').write_text(json.dumps(marker,indent=2),encoding='utf-8')
    (target/'SCIENTIFIC_EVALUATION_COMPLETE.json').write_text(json.dumps(marker,indent=2),encoding='utf-8')
    artifacts=[dict(path=str(p),sha256=sha(p),bytes=p.stat().st_size) for p in sorted(target.iterdir()) if p.is_file()]
    (target/'OUTPUT_MANIFEST.json').write_text(json.dumps(artifacts,indent=2),encoding='utf-8')
    print(json.dumps(marker))

if __name__=='__main__':
    assert len(sys.argv)==2, 'Usage: summarize.py OUTPUT_ROOT'
    summarize(sys.argv[1])
