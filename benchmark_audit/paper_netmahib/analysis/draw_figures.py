"""Draw all manuscript graphics directly from executed result files."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Circle

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'paper/figures';OUT.mkdir(exist_ok=True)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.titlesize':10,
 'axes.labelsize':9,'xtick.labelsize':8,'ytick.labelsize':8,'pdf.fonttype':42,
 'ps.fonttype':42,'axes.spines.top':False,'axes.spines.right':False,
 'savefig.facecolor':'white','figure.dpi':150})
COL={'release_order':'#426B91','peptide_sorted':'#CC6644','row_shuffled':'#B39742',
     'rejection':'#8064A2','matching':'#268477'}
LABEL={'release_order':'Release order','peptide_sorted':'Peptide sorted','row_shuffled':'Row shuffled',
       'rejection':'Within-HLA rejection','matching':'Closed matching'}
def save(fig,name):
    for ext in ['pdf','eps','png']:fig.savefig(OUT/f'{name}.{ext}',bbox_inches='tight',dpi=240)
    plt.close(fig)
def panel(ax,letter,title):
    ax.set_title(f'{letter}  {title}',loc='left',fontweight='bold',pad=12)
def clean(ax):
    ax.grid(axis='y',color='#E5E8EC',lw=.6);ax.set_axisbelow(True)

# Mechanism and exact assignment schematic. The graph is a declared synthetic example.
res=json.loads((ROOT/'audit_work/source_audit_results.json').read_text())['synthetic_batch_order_simulation']
fig,axs=plt.subplots(1,2,figsize=(7.05,3.1),gridspec_kw={'width_ratios':[1,1.15]})
ax=axs[0];y=[res[x]['majority_score_auroc']['mean'] for x in ['peptide_sorted','row_shuffled']]+[.5]
sd=[res[x]['majority_score_auroc']['sd'] for x in ['peptide_sorted','row_shuffled']]+[0]
ax.bar(range(3),np.array(y)-.48,bottom=.48,color=[COL['peptide_sorted'],COL['row_shuffled'],COL['matching']],width=.58)
ax.errorbar(range(3),y,yerr=sd,fmt='none',ecolor='black',capsize=3,lw=.8)
for i,v in enumerate(y):ax.text(i,v+.011,f'{v:.3f}',ha='center',fontsize=9)
ax.set_xticks(range(3),['Sorted','Shuffled','Exact cyclic\ncontrol']);ax.set_ylim(.48,.84)
ax.axhline(.5,color='#555555',ls='--',lw=.8);ax.set_ylabel('Fixed peptide-only score AUROC');clean(ax)
panel(ax,'a','A shortcut without a binding model')
ax.text(.03,.94,'Synthetic: 4,096 rows; 9 peptides\n100 sampler seeds; batch size 64',transform=ax.transAxes,fontsize=8,va='top')
ax=axs[1];ax.set_xlim(0,1);ax.set_ylim(0,1);ax.axis('off');panel(ax,'b','Closed assignment with abstention')
pos={'1':(.25,.77),'2':(.64,.77),'3':(.44,.42),'4':(.85,.35)}
for a,b in [('1','2'),('2','3'),('3','1')]:
    ax.add_patch(FancyArrowPatch(pos[a],pos[b],arrowstyle='-|>',mutation_scale=13,shrinkA=15,shrinkB=15,lw=1.8,color=COL['matching']))
for k,(x,y0) in pos.items():
    ax.add_patch(Circle((x,y0),.056,facecolor='#ECF4F2' if k!='4' else '#F6E8E2',edgecolor=COL['matching'] if k!='4' else COL['peptide_sorted'],lw=1.3))
    ax.text(x,y0,k,ha='center',va='center',fontsize=10)
ax.text(.44,.96,'One HLA stratum; nodes are positive rows',ha='center',fontsize=8)
ax.text(.27,.18,'Rows 1–3: retain cycle\nSame receptor and donor row set',ha='center',fontsize=8,color=COL['matching'])
ax.text(.85,.18,'Row 4:\nabstain',ha='center',fontsize=8,color=COL['peptide_sorted'])
ax.text(.48,.015,'Allowed edge: cost 0    Diagonal: 1    Forbidden: n + 1',ha='center',fontsize=8)
fig.tight_layout(w_pad=2);save(fig,'fig1_mechanism')

data=pd.read_csv(ROOT/'public_audit/seed_results.csv')
common=data[data.cohort.eq('score_ge_1_retained')]
# Same retained cohort; independent rejection comparator appended when available.
groups=[common[common.order.eq(x)] for x in ['release_order','peptide_sorted','row_shuffled']]
names=['release_order','peptide_sorted','row_shuffled']
reject_path=ROOT/'audit_work/reference_rejection_seed_results.csv'
if reject_path.exists():groups.append(pd.read_csv(reject_path));names.append('rejection')
groups.append(common[common.sampler.eq('exact_closed_matching')]);names.append('matching')
fig,axs=plt.subplots(1,3,figsize=(7.05,3.15))
for ax,metric,scale,title in zip(axs,['known_positive_collisions_fraction','peptide_tv','peptide_identity_oracle_auc'],[100,1,1],['Registry contradictions (%)','Peptide total variation','Identity-oracle AUROC']):
    vals=np.array([g[metric].mean()*scale for g in groups]);sds=np.array([g[metric].std(ddof=1)*scale if len(g)>1 else 0 for g in groups])
    yy=np.arange(len(names));ax.barh(yy,vals,color=[COL[n] for n in names],height=.6)
    ax.errorbar(vals,yy,xerr=sds,fmt='none',ecolor='#222222',capsize=2,lw=.7)
    ax.set_yticks(yy,[LABEL[n] for n in names] if ax is axs[0] else ['']*len(names));ax.invert_yaxis()
    ax.set_xlabel(title);ax.grid(axis='x',color='#E5E8EC',lw=.6);ax.set_axisbelow(True)
    if metric.endswith('auc'):ax.set_xlim(.48,.83);ax.axvline(.5,color='#555555',ls='--',lw=.8)
    else:ax.set_xlim(0,max(vals)*1.24 if max(vals)>0 else 1)
    for j,v in enumerate(vals):ax.text(v+sds[j]+.015*(ax.get_xlim()[1]-ax.get_xlim()[0]),j,f'{v:.3f}' if scale==1 else f'{v:.2f}',va='center',fontsize=7)
for ax,l in zip(axs,'abc'):ax.set_title(l,loc='left',fontweight='bold')
fig.suptitle('Identical retained positives: 2,691 rows',fontsize=10,y=1.03)
fig.tight_layout(w_pad=.7);save(fig,'fig2_common_cohort')

# Evidence-threshold sensitivity on the full (not matched) cohorts.
fig,axs=plt.subplots(1,2,figsize=(7.05,3.0))
cohorts=['score_ge_0','score_ge_1','score_eq_3'];orders=['release_order','peptide_sorted','row_shuffled']
for ax,metric,title,fmt in zip(axs,['known_positive_collisions_fraction','peptide_identity_oracle_auc'],['Registry contradictions (%)','Peptide identity-oracle AUROC'],['.2f','.3f']):
    matrix=np.array([[data.loc[data.cohort.eq(c)&data.order.eq(o)&data.batch_size.eq(64),metric].mean() for c in cohorts] for o in orders])
    if metric.endswith('fraction'):matrix*=100
    im=ax.imshow(matrix,cmap='YlOrRd' if metric.endswith('fraction') else 'Blues',aspect='auto',vmin=0 if metric.endswith('fraction') else .5)
    ax.set_xticks(range(3),['Score ≥ 0\nn = 29,306','Score ≥ 1\nn = 3,182','Score = 3\nn = 320'])
    ax.set_yticks(range(3),[LABEL[o] for o in orders] if ax is axs[0] else ['']*3)
    for i in range(3):
        for j in range(3):ax.text(j,i,format(matrix[i,j],fmt),ha='center',va='center',color='white' if matrix[i,j]> (matrix.max()+matrix.min())*.5 else '#17212B',fontsize=9)
    ax.tick_params(length=0);panel(ax,'a' if ax is axs[0] else 'b',title)
fig.tight_layout(w_pad=1.7);save(fig,'fig3_confidence')

# Coverage loss is shown, not hidden by the common-cohort comparison.
ret=pd.read_csv(ROOT/'public_audit/matching_retention_by_hla.csv').sort_values('n_positive',ascending=False)
fig,axs=plt.subplots(1,2,figsize=(7.05,3.8),gridspec_kw={'width_ratios':[1.4,1]})
ax=axs[0];sub=ret.head(12).iloc[::-1];yy=np.arange(len(sub));ax.barh(yy,sub.n_retained,color=COL['matching'],label='Retained');ax.barh(yy,sub.n_positive-sub.n_retained,left=sub.n_retained,color='#DBCAC2',label='Excluded')
ax.set_yticks(yy,[x.replace('HLA-','') for x in sub.hla]);ax.set_xlabel('Positive records');panel(ax,'a','Coverage in the 12 largest HLA strata');ax.legend(frameon=False,fontsize=8,loc='lower right')
meta=json.loads((ROOT/'public_audit/split_sensitivity/split_provenance.json').read_text())['splits']
ax=axs[1];coverage=[100*x['test_retention_fraction'] for x in meta];ax.bar(range(5),coverage,color=COL['matching'],width=.65)
for i,v in enumerate(coverage):ax.text(i,v+2,f'{v:.1f}',ha='center',fontsize=8)
ax.set_ylim(0,100);ax.set_xticks(range(5),range(5));ax.set_xlabel('Peptide-component split seed');ax.set_ylabel('Retained test positives (%)');panel(ax,'b','Coverage after distance splitting');clean(ax)
ax.text(.02,.94,'All cross-partition\nedit distances ≥ 3',transform=ax.transAxes,va='top',fontsize=8)
fig.tight_layout(w_pad=2);save(fig,'fig4_coverage')

first=pd.read_csv(ROOT/'crossfit_audit/split_first_results/summary.csv');late=pd.read_csv(ROOT/'crossfit_audit/results/summary.csv')
fig,axs=plt.subplots(1,2,figsize=(7.05,3.1))
samplers=['source_release_order','source_peptide_sorted','source_row_shuffled','exact_closed_matching'];nm=['release_order','peptide_sorted','row_shuffled','matching']
for ax,frame,title in zip(axs,[first,late],['Split receptors, then construct decoys','Construct globally, then split receptors']):
    q=frame[frame.feature.eq('peptide_hla')].set_index('sampler').loc[samplers]
    ax.bar(range(4),q.auroc_mean,color=[COL[x] for x in nm],width=.6)
    ax.errorbar(range(4),q.auroc_mean,yerr=q.auroc_sd_between_sampler_seeds.fillna(0),fmt='none',ecolor='#222222',capsize=3,lw=.8)
    for i,v in enumerate(q.auroc_mean):ax.text(i,v+.014,f'{v:.3f}',ha='center',fontsize=9)
    ax.set_xticks(range(4),['Release','Sorted','Shuffled','Matching']);ax.set_ylim(.25,.9);ax.axhline(.5,color='#555555',ls='--',lw=.8);clean(ax)
    panel(ax,'a' if ax is axs[0] else 'b',title)
axs[0].set_ylabel('Peptide + HLA held-out-receptor AUROC')
axs[1].text(.03,.95,'Protocol-placement sensitivity:\ncomplementary folds are coupled',transform=axs[1].transAxes,va='top',fontsize=8)
fig.tight_layout(w_pad=1.3);save(fig,'fig5_partition_order')

# Supplement: batch sizes, no interpolation between tested conditions.
fig,axs=plt.subplots(1,3,figsize=(7.05,2.8))
for ax,metric,scale,label in zip(axs,['known_positive_collisions_fraction','peptide_tv','fallback_count_fraction'],[100,1,100],['Registry contradictions (%)','Peptide total variation','Fallback use (%)']):
    for o in orders:
        q=data[data.cohort.eq('score_ge_1')&data.order.eq(o)].groupby('batch_size')[metric].agg(['mean','std']).loc[[16,64,256]]
        ax.errorbar(range(3),q['mean']*scale,yerr=q['std']*scale,marker='o',ms=4,lw=1.2,capsize=2,color=COL[o],label=LABEL[o])
    ax.set_xticks(range(3),[16,64,256]);ax.set_xlabel('Batch size');ax.set_ylabel(label);clean(ax)
axs[1].legend(frameon=False,fontsize=7,loc='upper right');fig.tight_layout(w_pad=1);save(fig,'figS1_batch_sizes')

# Author-facing geometric counterexample: not empirical evidence about neural attention.
geo=pd.read_csv(ROOT/'research/results/geometry_control.csv')
print('Drew manuscript figures:',len(list(OUT.glob('fig*.pdf'))))
