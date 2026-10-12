from pathlib import Path
import pandas as pd
import json
ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'paper/tables';OUT.mkdir(exist_ok=True)
def table(name,header,rows,align=None):
    if align is None:align='l'+'r'*(len(header)-1)
    s='\\begin{tabular}{'+align+'}\n\\toprule\n'+' & '.join(header)+'\\\\\n\\midrule\n'
    s+='\n'.join(' & '.join(map(str,r))+'\\\\' for r in rows)+'\n\\bottomrule\n\\end{tabular}\n'
    (OUT/name).write_text(s)
flow=pd.read_csv(ROOT/'public_audit/filter_flow.csv').iloc[:7]
names=['Frozen release','Human species','TRB chain','Class-I MHC','Two-field HLA source-label syntax','Canonical CDR3-beta','Canonical peptide, length 8--15']
table('filter_flow.tex',['Filter retained sequentially','Rows'],[(label,f'{int(n):,}') for label,n in zip(names,flow.n)])
meta=json.loads((ROOT/'public_audit/split_sensitivity/split_provenance.json').read_text())['splits']
table('distance_splits.tex',['Seed','Train','Validation','Test','Retained test','Retained (\\%)'],[(x['split_seed'],x['partition_row_counts']['train'],x['partition_row_counts']['validation'],x['partition_row_counts']['test'],x['test_retained_n'],f"{100*x['test_retention_fraction']:.2f}") for x in meta])
rows=[]
for protocol,path in [('Split first','split_first_results'),('Construct first','results')]:
    d=pd.read_csv(ROOT/f'crossfit_audit/{path}/summary.csv')
    for sampler,label in [('source_release_order','Release'),('source_peptide_sorted','Sorted'),('source_row_shuffled','Shuffled'),('exact_closed_matching','Assignment')]:
        x=d[d.feature.eq('peptide_hla')&d.sampler.eq(sampler)].iloc[0]
        rows.append([protocol,label,f'{x.auroc_mean:.4f}',f'{x.ap_mean:.4f}',f'{x.brier_mean:.4f}'])
table('trained_controls.tex',['Protocol','Sampler','AUROC','AP','Brier'],rows,'llrrr')
d=pd.read_csv(ROOT/'public_audit/matching_retention_by_hla.csv').sort_values('hla')
rows=[]
for x in d.itertuples():rows.append([x.hla.replace('*',r'\textasteriskcentered{}'),x.n_positive,x.n_retained,x.n_distinct_peptides,f'{100*x.n_retained/x.n_positive:.1f}'])
s='\\begin{longtable}{lrrrr}\n\\caption{Complete retention by source HLA annotation stratum. Flagged labels A*08:01 and A*24:01 are discussed in Supplementary Methods.}\\label{tab:fullhla}\\\\\n\\toprule\nSource label & Original & Retained & Peptides & Retained (\\%)\\\\\n\\midrule\\endfirsthead\n\\toprule Source label & Original & Retained & Peptides & Retained (\\%)\\\\\\midrule\\endhead\n'
s+='\n'.join(' & '.join(map(str,r))+'\\\\' for r in rows)+'\n\\bottomrule\n\\end{longtable}\n'
(OUT/'hla_retention.tex').write_text(s)
pub=pd.read_csv(ROOT/'public_audit/aggregate_results.csv')
rows=[]
for cohort in ['score_ge_0','score_ge_1','score_eq_3']:
    for order,label in [('release_order','Release'),('peptide_sorted','Sorted'),('row_shuffled','Shuffled')]:
        g=pub[pub.cohort.eq(cohort)&pub.order.eq(order)&pub.batch_size.eq(64)].set_index('metric')
        vals=[]
        for m,scale in [('known_positive_collisions_fraction',100),('peptide_tv',1),('peptide_identity_oracle_auc',1)]:
            vals.append(f"${scale*g.loc[m,'mean']:.4f}\\pm{scale*g.loc[m,'sd']:.4f}$")
        rows.append([{'score_ge_0':'$\\geq0$','score_ge_1':'$\\geq1$','score_eq_3':'$=3$'}[cohort],label]+vals)
table('confidence_full.tex',['Score','Order','Contradictions (\\%)','Peptide TV','Oracle AUROC'],rows,'llrrr')
print('Generated',len(list(OUT.glob('*.tex'))),'tables')
