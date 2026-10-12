from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE=Path(__file__).resolve().parent
order=['source_release_order','source_peptide_sorted','source_row_shuffled','exact_closed_matching']
labels=['Source order','Peptide sorted','Row shuffled','Closed matching']
colors=['#54758B','#C76332','#40877A','#635589']
fig,axes=plt.subplots(1,2,figsize=(10,3.8),sharey=True)
for ax,folder,title in zip(axes,['results','split_first_results'],
 ['A  Negatives before receptor folds','B  Negatives within receptor folds']):
 d=pd.read_csv(HERE/folder/'oof_metrics.csv').query("feature=='peptide_hla'")
 for i,key in enumerate(order):
  vals=d.loc[d.sampler.eq(key),'auroc'].to_numpy()
  jitter=np.random.default_rng(99+i).uniform(-.11,.11,len(vals))
  ax.scatter(i+jitter,vals,s=10,alpha=.35,color=colors[i],edgecolors='none')
  ax.plot([i-.18,i+.18],[vals.mean(),vals.mean()],color=colors[i],lw=2.5)
  ax.text(i,vals.max()+.021,f'{vals.mean():.3f}',ha='center',fontsize=9,color=colors[i])
 ax.axhline(.5,color='#777777',ls='--',lw=1,zorder=0)
 ax.set_xticks(range(4),labels,rotation=25,ha='right')
 ax.set_title(title,loc='left',fontsize=10,fontweight='bold')
 ax.spines[['top','right']].set_visible(False)
 ax.set_ylim(.26,.8)
 ax.grid(axis='y',alpha=.12)
axes[0].set_ylabel('Out-of-fold AUROC (peptide + HLA control)')
fig.text(.5,.015,'Within-antigen-collection construction diagnostic; no receptor sequence input. Points are fold-map / decoy-seed runs.',
 ha='center',fontsize=8,color='#444444')
fig.tight_layout(rect=(0,.04,1,1))
fig.savefig(HERE/'metadata_crossfit_protocol.pdf',bbox_inches='tight')
fig.savefig(HERE/'metadata_crossfit_protocol.png',dpi=180,bbox_inches='tight')
