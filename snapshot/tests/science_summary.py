"""Synthetic output-contract checks, including unavailable winners; zero fits."""
import csv, importlib.util, json, pathlib, sys, tempfile
sys.dont_write_bytecode=True

root=pathlib.Path(sys.argv[1]).resolve()
testout=pathlib.Path(sys.argv[2]).resolve() if len(sys.argv)>2 else root/'tests'
if len(sys.argv)>2:testout.mkdir(parents=True,exist_ok=False)
spec=importlib.util.spec_from_file_location('summary',root/'science/summarize.py')
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
def write(path,rows,fields=None):
    path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fields or (list(rows[0]) if rows else []));w.writeheader();w.writerows(rows)
def read(path):
    with path.open(encoding='utf-8',newline='') as f:return list(csv.DictReader(f))
plan=read(root/'config/PLAN.csv');ids=list(dict.fromkeys(r['dataset_id'] for r in plan))
checks=[]
with tempfile.TemporaryDirectory(prefix='science_summary_',dir=testout) as tmp:
    for case,missing in [('all_available',set()),('one_missing',{ids[-1]}),('all_missing',set(ids))]:
        out=pathlib.Path(tmp)/case;out.mkdir();write(out/'PLAN.csv',plan)
        for did in ids:
            dest=out/'science'/did;dest.mkdir(parents=True)
            unavailable=did in missing
            status=dict(status='SCIENCE_COMPLETE',data_id=did,fit_updates=0,truth_used_for_selection=False,selection_frozen_before_truth=True,registered=12,seconds=0)
            (dest/'SCIENCE_COMPLETE.json').write_text(json.dumps(status))
            write(dest/'POOL_SCIENCE.csv',[dict(data_id=did,arm='original',registered=12,valid=0 if unavailable else 12,
                failed_or_invalid=12 if unavailable else 0,correct_candidates=0 if unavailable else 3,
                pool_has_correct_structure=not unavailable,winner_exists=not unavailable,winner_correct=not unavailable,
                winner_shared='' if unavailable else 2,winner_specific_1='' if unavailable else 2,winner_specific_2='' if unavailable else 2)])
            ranks=[]
            if not unavailable:
                for block,study in [('shared',0),('specific',1),('specific',2)]:
                    for factor in (1,2):
                        ranks.append(dict(data_id=did,arm='original',block=block,study=study,true_factor=factor,
                                          matched_fpca_M_09=2 if factor==1 else ('' if study==2 else 3)))
            write(dest/'WINNER_FPCA_ROLE_RANKS.csv',ranks,['data_id','arm','block','study','true_factor','matched_fpca_M_09'])
            for file,keys in [('WINNER_RAW_COMPONENT_METRICS.csv',dict(scope='grid201',object='signal')),
                              ('WINNER_Q_INTERVAL_METRICS.csv',dict(scope='grid201',object='signal')),
                              ('WINNER_HISTORICAL_SUMMARY.csv',{}),('WINNER_SCIENCE.csv',{})]:
                values=[] if unavailable else [dict(data_id=did,arm='original',**keys,value=0.25)]
                write(dest/file,values,['data_id','arm',*keys,'value'])
            if not unavailable:
                folder=dest/'original';folder.mkdir()
                for filename in ['ACTIVE_INTERVALS.csv','CENTERED_SIGNAL.csv','NOISE.csv','COMPONENT_METRICS.csv','LOADING_METRICS.csv','BLOCK_SUBSPACE.csv']:
                    write(folder/filename,[])
        module.summarize(out)
        summary=json.loads((out/'aggregate/SUMMARY.json').read_text())
        assert summary['registered_full_fits']==120 and summary['new_MGP_fits_registered']==0
        structures=read(out/'aggregate/STRUCTURE_RECOVERY_SUMMARY.csv');assert len(structures)==1
        assert int(structures[0]['registered_starts'])==120
        assert float(structures[0]['winners_available'])==10-len(missing)
        assert int(structures[0]['failed_or_invalid'])==12*len(missing)
        counts=read(out/'aggregate/FACTOR_COUNT_SUMMARY.csv')
        assert all(int(r['missing'])==len(missing) for r in counts)
        assert {r['role']:int(r['truth']) for r in counts}==dict(shared=2,specific_1=2,specific_2=2,total=6)
        ranks=read(out/'aggregate/FPCA_ROLE_RANKS.csv');assert len(ranks)==60
        assert sum(r['category']=='missing_winner' for r in ranks)==6*len(missing)
        assert sum(r['category']=='missing_role' for r in ranks)==10-len(missing)
        freqs=[r for r in read(out/'aggregate/FPCA_RANK_FREQUENCIES.csv') if r['role']=='all_true_roles']
        assert sum(int(r['count']) for r in freqs)==60
        assert all(int(r['registered_true_roles'])==60 for r in freqs)
        scientific=read(out/'aggregate/SCIENTIFIC_METRIC_SUMMARY.csv')
        if case=='all_missing':assert scientific==[]
        else:assert all(int(r['finite_n'])==10-len(missing) for r in scientific if r['mean'])
        try:module.summarize(out)
        except AssertionError as e:assert 'overwrite' in str(e)
        else:raise AssertionError('overwrite was accepted')
        checks.append(dict(case=case,status='PASS',registered_datasets=10,registered_starts=120,true_role_denominator=60))
(testout/'science_SUMMARY_PASS.json').write_text(json.dumps(dict(status='PASS',cases=checks,new_fits=0),indent=2),encoding='utf-8')
print('SCIENCE_SUMMARY_PASS',len(checks),'cases')
