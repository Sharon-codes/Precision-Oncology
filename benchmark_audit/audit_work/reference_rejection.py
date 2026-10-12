#!/usr/bin/env python3
"""Explicit rejection reference on the same fixed retained positive cohort.

For each positive row independently, draw uniformly WITH REPLACEMENT from
retained donor rows of the same HLA whose beta/peptide/HLA triple is absent
from the declared all-score known-positive registry. This is a new, explicitly
defined reference construction, not an official biological negative sampler.
No prediction model is trained or tested.
"""
import argparse, hashlib, importlib.util, json, sys
from collections import Counter,defaultdict
from pathlib import Path
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parent

def load_public(path):
    spec=importlib.util.spec_from_file_location('public_reference',path);module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);return module

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--source',type=Path,default=ROOT.parent/'research/data/vdjdb.txt');ap.add_argument('--public-source',type=Path,default=ROOT.parent/'public_audit/run_public_audit.py');ap.add_argument('--manifest',type=Path,default=ROOT.parent/'public_audit/manifests/common_corrected_seed0.csv');ap.add_argument('--seeds',type=int,default=30);ap.add_argument('--out',type=Path,default=ROOT);a=ap.parse_args();a.out.mkdir(parents=True,exist_ok=True)
    public=load_public(a.public_source);_,known,known_pair,_=public.prepare(a.source)
    raw=pd.read_csv(a.source,sep='\t',keep_default_na=False);manifest=pd.read_csv(a.manifest)
    ix=manifest.receiver_source_row_0based.to_numpy(dtype=int)
    d=raw.iloc[ix][['cdr3','antigen.epitope','mhc.a']].rename(columns={'cdr3':'beta','antigen.epitope':'peptide','mhc.a':'hla'}).reset_index(drop=True)
    d['source_row_0based']=ix
    assert len(d)==2691 and len(set(ix))==len(ix)
    assert all(t in known for t in d[['beta','peptide','hla']].itertuples(index=False,name=None))
    eligible={};donor_by_hla={h:sub.index.to_numpy() for h,sub in d.groupby('hla')}
    for b,h in d[['beta','hla']].itertuples(index=False,name=None):
        if (b,h) not in eligible:
            eligible[b,h]=np.array([j for j in donor_by_hla[h] if (b,d.loc[j,'peptide'],h) not in known],dtype=int)
            assert len(eligible[b,h])>0
    # Exact expected class marginals (sampling is over donor ROWS, so repeated
    # peptide identities retain their multiplicity in these probabilities).
    expected=defaultdict(float)
    for b,h in d[['beta','hla']].itertuples(index=False,name=None):
        candidates=eligible[b,h];cnt=Counter(d.loc[candidates,'peptide'])
        for p,c in cnt.items():expected[p,h]+=c/len(candidates)
    expected_pep=defaultdict(float)
    for (p,h),value in expected.items():expected_pep[p]+=value
    observed=Counter(d.peptide);observed_ph=Counter(zip(d.peptide,d.hla));n=len(d)
    exp_tv=.5*sum(abs(observed[p]-expected_pep[p]) for p in set(observed)|set(expected_pep))/n
    exp_ph_tv=.5*sum(abs(observed_ph[p]-expected[p]) for p in set(observed_ph)|set(expected))/n
    rows=[];md=a.out/'reference_rejection_manifests';md.mkdir(exist_ok=True)
    for seed in range(a.seeds):
        rng=np.random.default_rng(seed)
        donor_indices=np.array([int(rng.choice(eligible[b,h])) for b,h in d[['beta','hla']].itertuples(index=False,name=None)])
        neg=[(d.loc[i,'beta'],d.loc[j,'peptide'],d.loc[j,'hla']) for i,j in enumerate(donor_indices)]
        result=public.metrics(d,neg,known,known_pair)
        assert result['known_positive_collisions']==0
        assert result['hla_tv']==result['beta_tv']==0
        assert all(d.loc[i,'hla']==t[2] for i,t in enumerate(neg))
        result.update(seed=seed,unique_negative_triplets=len(set(neg)),duplicate_negative_rows=n-len(set(neg)),unique_donor_rows=len(set(donor_indices)))
        rows.append(result)
        pd.DataFrame({'receiver_source_row_0based':ix,'negative_antigen_source_row_0based':d.loc[donor_indices,'source_row_0based'].to_numpy()}).to_csv(md/f'reference_rejection_seed{seed}.csv',index=False)
    frame=pd.DataFrame(rows);frame.to_csv(a.out/'reference_rejection_seed_results.csv',index=False)
    summary={}
    for col in frame.columns:
        if col=='seed':continue
        x=frame[col].to_numpy(dtype=float);summary[col]=dict(mean=float(x.mean()),sd=float(x.std(ddof=1)) if len(x)>1 else 0,min=float(x.min()),max=float(x.max()))
    matched_ix=manifest.negative_antigen_source_row_0based.to_numpy(dtype=int)
    matched=[(d.loc[i,'beta'],raw.iloc[j]['antigen.epitope'],raw.iloc[j]['mhc.a']) for i,j in enumerate(matched_ix)]
    match_metrics=public.metrics(d,matched,known,known_pair)
    match_metrics.update(unique_negative_triplets=len(set(matched)),duplicate_negative_rows=n-len(set(matched)))
    result=dict(scope='New explicitly defined reference rejection construction; no official-model or biological-negative claims.',
      protocol='Independently for every retained receiver, uniformly sample one eligible retained donor ROW with replacement within the same HLA, excluding every registry-positive beta/peptide/HLA triple.',
      source_sha256=hashlib.sha256(a.source.read_bytes()).hexdigest(),public_source_sha256=hashlib.sha256(a.public_source.read_bytes()).hexdigest(),retained_manifest_sha256=hashlib.sha256(a.manifest.read_bytes()).hexdigest(),
      n_positive=n,n_hla=int(d.hla.nunique()),n_seeds=a.seeds,seeds=list(range(a.seeds)),rng='numpy.default_rng(seed)',
      registry='Same all-score eligible human-TRB/MHCI/exact-two-field-HLA reference used by public audit.',
      assertions_passed=['every_receiver_has_eligible_donor','no_known_positive_collision','exact_HLA_marginal_equality','exact_beta_marginal_equality','within_HLA_pairing'],
      exact_expected_peptide_tv=float(exp_tv),exact_expected_peptide_hla_tv=float(exp_ph_tv),
      summary=summary,closed_matching_same_cohort=match_metrics,
      limits=['Expected marginal TV is positive because conditioning on each receiver alters donor inclusion probabilities; with-replacement sampling adds finite-sample fluctuation.',
      'Identity-oracle AUROC is a same-dataset resubstitution marginal diagnostic, not held-out prediction performance.',
      'Zero registry collisions does not establish biological nonbinding.',
      'This comparison uses the already-retained common cohort and does not resolve the matching selection cost.'])
    (a.out/'reference_rejection_results.json').write_text(json.dumps(result,indent=2)+'\n')
    keep=['peptide_tv','peptide_identity_oracle_auc','peptide_hla_tv','peptide_hla_identity_oracle_auc','known_positive_collisions','hla_tv','beta_tv','unique_negative_triplets','duplicate_negative_rows','unique_donor_rows']
    print(json.dumps({k:summary[k] for k in keep},indent=2));print('EXPECTED_TV',exp_tv,exp_ph_tv)
if __name__=='__main__':main()
