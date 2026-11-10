import numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.path import Path
from matplotlib.patches import PathPatch
from matplotlib.lines import Line2D
plt.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Liberation Sans','Arial','DejaVu Sans'],'pdf.fonttype':42,'axes.linewidth':0.5})
INK='#1a1a1a'; INK2='#52514e'
TFD='#1b5490'; TFL='#2a78d6'; GEX='#eb6834'
from keys import spec
from base import C,E
SP={t:spec(t)[3].set_index(['TF','cell_type']) for t in ['Normal','Adenoma']}
MODS=[('Normal','ABS','Absorptive'),('Normal','GOB','Goblet'),
      ('Adenoma','ASC','Adenoma-specific'),('Adenoma','SSC','Serrated-specific')]

def module(tis,ct,maxtf=6):
    r=SP[tis].reset_index()
    s=r[(r.cell_type==ct)&(r.p<0.05)&(r.n_linked>=3)].sort_values('p').head(maxtf)
    tfs=list(s.TF)
    d=E[(E.tissue==tis)&(E.cell_type==ct)&(E.TF.isin(tfs))]
    cnt=d.groupby('gene').TF.nunique()
    genes=set(cnt[cnt>=2].index)
    if len(genes)<5: genes=set(cnt.index)
    for tf in tfs:
        sub=d[d.TF==tf]
        if not set(sub.gene)&genes: genes|=set(sub.sort_values('r',ascending=False).gene.head(2))
    genes=sorted(genes)
    return tfs,genes,d[d.gene.isin(genes)],s.set_index('TF')

def draw(ax,tis,ct,full):
    tfs,genes,pairs,S=module(tis,ct)
    gct=C[C.tissue==tis].groupby('gene').cell_type.nunique()
    ty=np.linspace(0,1,len(tfs)+2)[1:-1] if len(tfs)>1 else np.array([0.5])
    gy=np.linspace(0,1,len(genes)+2)[1:-1] if len(genes)>1 else np.array([0.5])
    TP=dict(zip(tfs,ty)); GP=dict(zip(genes,gy))
    for _,p in pairs.sort_values('r').iterrows():
        y0=TP[p.TF]; y1=GP[p.gene]; strong=p.r>0.7
        ax.add_patch(PathPatch(Path([(0,y0),(0.5,y0),(0.5,y1),(1,y1)],
            [Path.MOVETO,Path.CURVE4,Path.CURVE4,Path.CURVE4]),fill=False,
            lw=0.55 if strong else 0.30,color='#1b4f72' if strong else '#a7aeb5',
            alpha=0.72 if strong else 0.30,zorder=3 if strong else 2))
    for t in tfs:
        y=TP[t]; row=S.loc[t]; fdr=row.q<0.05
        ax.scatter([0],[y],s=34,c=[TFD if fdr else TFL],marker='s',edgecolors='white',linewidths=0.6,zorder=4)
        ax.text(-0.045,y,t,ha='right',va='center',fontsize=5.4,color=INK,fontweight='bold',zorder=5)
    for g in genes:
        y=GP[g]; ex=gct[g]==1
        ax.scatter([1],[y],s=16,c=[GEX if ex else 'white'],edgecolors='#4f5a65',linewidths=0.45,zorder=4)
        ax.text(1.042,y,g,ha='left',va='center',fontsize=4.9,color=INK,zorder=5)
    ax.set_xlim(-0.62,1.62); ax.set_ylim(-0.055,1.055); ax.axis('off')
    ax.set_title(f'{full} ({ct}) · {tis.lower()}\n{len(pairs)} links · {len(tfs)} TFs · {len(genes)} risk genes',
                 fontsize=6.6,color=INK,fontweight='bold',pad=3)

fig=plt.figure(figsize=(12,4.1))
gs=fig.add_gridspec(1,4,wspace=0.06,hspace=0.20,left=0.02,right=0.985,top=0.930,bottom=0.070)
for k,(tis,ct,full) in enumerate(MODS):
    ax=fig.add_subplot(gs[0,k]); draw(ax,tis,ct,full)
    ax.text(-0.02,1.09,'fghi'[k],transform=ax.transAxes,fontsize=9.5,fontweight='bold',va='bottom')
h=[Line2D([],[],marker='s',ls='',mfc=TFD,mec='white',ms=5,label='Cell-type-enriched TF, FDR < 5%'),
   Line2D([],[],marker='s',ls='',mfc=TFL,mec='white',ms=5,label='Cell-type-enriched TF, $P$ < 0.05'),
   Line2D([],[],marker='o',ls='',mfc=GEX,mec='#4f5a65',mew=0.45,ms=4,label='Risk gene linked only in this cell type'),
   Line2D([],[],marker='o',ls='',mfc='white',mec='#4f5a65',mew=0.45,ms=4,label='Risk gene shared with other cell types'),
   Line2D([],[],color='#1b4f72',lw=1.0,alpha=0.75,label='$r$ > 0.7'),
   Line2D([],[],color='#a7aeb5',lw=0.6,alpha=0.5,label=r'$0.5 < r \leq 0.7$')]
fig.legend(handles=h,loc='lower center',ncol=3,frameon=False,fontsize=5.9,
           handletextpad=0.35,columnspacing=1.5,bbox_to_anchor=(0.5,0.004))
fig.savefig('Figure4fi_celltype_TF_modules.pdf'); fig.savefig('Figure4fi_celltype_TF_modules.png',dpi=230)
print('ok')
