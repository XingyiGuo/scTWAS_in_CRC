import numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.lines import Line2D
import matplotlib.cm as cm
from spec import build
plt.rcParams.update({'font.family':'Liberation Sans','pdf.fonttype':42,'axes.linewidth':0.5,
 'xtick.labelsize':6,'ytick.labelsize':6,'axes.labelsize':7})
INK='#1a1a1a'; INK2='#52514e'
SEQ=LinearSegmentedColormap.from_list('p',['#e9eef5','#a8c6e8','#4a8ed6','#1b5490','#0c2f52'])
ORDER={'Normal':['STM','TAC','ABS','GOB','EE','TUF','CT'],
       'Adenoma':['STM','TAC','ABS','GOB','ASC','SSC']}

def panel(ax,tis,show_y=True):
    cts,ng,deg,rate,R,S=build(tis)
    cts=[c for c in ORDER[tis] if c in cts]
    R=R[R.cell_type.isin(cts)]
    hits=R[(R.p<0.05)&(R.n_linked>=3)]
    tfs=sorted(set(hits.TF))
    best=hits.sort_values('p').drop_duplicates('TF').set_index('TF')
    tfs=sorted(tfs,key=lambda t:(cts.index(best.loc[t,'cell_type']),best.loc[t,'p']))
    RM=R.set_index(['TF','cell_type'])
    vmax=6.0
    for i,t in enumerate(tfs):
        for j,c in enumerate(cts):
            r=RM.loc[(t,c)]
            if r.n_linked==0: continue
            sz=8+330*r.rate
            lp=min(-np.log10(max(r.p,1e-12)),vmax)
            sig=(r.q<0.05)and(r.n_linked>=3)
            ax.scatter([j],[i],s=sz,c=[SEQ(lp/vmax)],edgecolors=('#101010' if sig else '#b8b7b2'),
                       linewidths=(0.8 if sig else 0.35),zorder=3)
    # dominant-cell-type separators
    doms=[cts.index(best.loc[t,'cell_type']) for t in tfs]
    for i in range(1,len(tfs)):
        if doms[i]!=doms[i-1]:
            ax.axhline(i-0.5,color='#c9c8c3',lw=0.5,zorder=1)
    ax.set_xticks(range(len(cts))); ax.set_xticklabels(cts,fontsize=6.4,rotation=0)
    ax.xaxis.set_ticks_position('top'); ax.xaxis.set_label_position('top')
    ax.set_yticks(range(len(tfs)))
    ax.set_yticklabels(tfs if show_y else [],fontsize=5.8)
    ax.set_ylim(len(tfs)-0.2,-0.8); ax.set_xlim(-0.7,len(cts)-0.3)
    ax.tick_params(length=0)
    ax.grid(axis='x',lw=0.4,color='#f2f1ed',zorder=0)
    ax.set_axisbelow(True)
    for s in ax.spines.values(): s.set_color('#d8d7d2')
    return cts,tfs

fig=plt.figure(figsize=(7.09,6.2))
gs=fig.add_gridspec(1,2,width_ratios=[1.0,0.92],wspace=0.30,left=0.085,right=0.985,top=0.915,bottom=0.115)
axa=fig.add_subplot(gs[0,0]); panel(axa,'Normal')
axa.set_title('Normal colon',fontsize=7.6,fontweight='bold',color=INK,pad=16)
fig.text(0.008,0.955,'a',fontsize=10,fontweight='bold',va='bottom')
axb=fig.add_subplot(gs[0,1]); panel(axb,'Adenoma')
axb.set_title('Adenoma and serrated lesions',fontsize=7.6,fontweight='bold',color=INK,pad=16)
fig.text(0.505,0.955,'b',fontsize=10,fontweight='bold',va='bottom')


cax=fig.add_axes([0.115,0.045,0.20,0.010])
cb=fig.colorbar(cm.ScalarMappable(cmap=SEQ,norm=Normalize(0,6)),cax=cax,orientation='horizontal')
cb.set_label('Enrichment of TF–risk gene links in cell type, $-$log$_{10}$ $P$',fontsize=5.6,labelpad=1.5,color=INK2)
cb.ax.tick_params(labelsize=5,length=2,width=0.4); cb.outline.set_linewidth(0.4)
cb.set_ticks([0,2,4,6]); cb.ax.set_xticklabels(['0','2','4','≥6'])

h=[Line2D([],[],marker='o',ls='',mfc='#cfd6de',mec='#b8b7b2',mew=0.35,ms=np.sqrt(8+330*v),
          label=f'{int(v*100)}%') for v in (0.05,0.15,0.30)]
l1=fig.legend(handles=h,loc='lower center',bbox_to_anchor=(0.60,0.008),ncol=3,frameon=False,
          fontsize=5.6,handletextpad=0.25,columnspacing=1.4,
          title='Risk genes of that cell type linked to TF',title_fontsize=5.8)
l1.get_title().set_color(INK2)
h2=[Line2D([],[],marker='o',ls='',mfc='#7fb0e0',mec='#101010',mew=0.8,ms=4,label='FDR < 5%')]
fig.legend(handles=h2,loc='lower center',bbox_to_anchor=(0.885,0.028),frameon=False,fontsize=5.8,handletextpad=0.25)

fig.savefig('Fig2_celltype_specific_TFs.pdf'); fig.savefig('Fig2_celltype_specific_TFs.png',dpi=230)
print('ok')
