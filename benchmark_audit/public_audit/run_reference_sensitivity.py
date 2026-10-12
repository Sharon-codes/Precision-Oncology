#!/usr/bin/env python3
"""Confidence-reference sensitivity for known-positive contradictions."""
import json
import pandas as pd
from run_public_audit import HERE,prepare,original,corrected,metrics,save_manifest

def main():
    c,known,kp,_=prepare(HERE.parent/'research/data/vdjdb.txt');main=c['score_ge_1']
    strong=set(main[['beta','peptide','hla']].itertuples(index=False,name=None));rows=[]
    for order in ['release_order','peptide_sorted','row_shuffled']:
        for seed in range(30):
            ordered,neg,_,_=original(main,HERE.parent/'precision_repo',order,64,seed,known)
            rows.append({'order':order,'seed':seed,'n_negative':len(neg),'all_eligible_known_collisions':sum(t in known for t in neg),'score_ge_1_known_collisions':sum(t in strong for t in neg)})
    pd.DataFrame(rows).to_csv(HERE/'confidence_collision_sensitivity.csv',index=False)
    r,n,s=corrected(c['score_eq_3'],known);pd.DataFrame(s).to_csv(HERE/'strict_matching_retention_by_hla.csv',index=False);save_manifest(HERE,'strict_corrected_seed0',r,n,c['score_ge_0'])
    r1,n1,s1=corrected(main,strong)
    (HERE/'additional_sensitivity.json').write_text(json.dumps({'score3_corrected':{'retained_n':len(r),'original_n':len(c['score_eq_3']),'metrics':metrics(r,n,known,kp)},'main_confidence_reference_only':{'retained_n':len(r1),'known_collisions_against_all_scores':sum(t in known for t in n1)},'dominant_peptide_fractions':{name:float(d.peptide.value_counts().iloc[0]/len(d)) for name,d in c.items()}},indent=2)+'\n')
    print(pd.DataFrame(rows).groupby('order')[['all_eligible_known_collisions','score_ge_1_known_collisions']].mean().to_string())
if __name__=='__main__':main()
