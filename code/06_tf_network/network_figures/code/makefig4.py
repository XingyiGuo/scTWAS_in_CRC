import numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.path import Path
from matplotlib.patches import PathPatch
from matplotlib.colors import LinearSegmentedColormap, Normalize
import matplotlib.cm as cm
from spec import build
plt.rcParams.update({'font.family':'Liberation Sans','pdf.fonttype':42,'axes.linewidth':0.5})
INK='#1a1a1a'; INK2='#52514e'
TFC='#2a78d6'; TFC2='#1b5490'
GEX='#eb6834'   # exclusive gene
GSH='#ffffff'   # shared gene

df=pd.read_csv('tf_gene_correlations_annot.csv'); e=df[df.r>0.5]
MODS=[('Normal','TAC','Transit-amplifying'),('Normal','ABS','Absorptive'),('Normal','GOB','Goblet'),
      ('Normal','CT','Crypt-top'),('Adenoma','ASC','Adenoma-specific'),('Adenoma','SSC','Serrated-specific')]
R={t:build(t)[4] for t in ['Normal','Adenoma']}

def module(tis,ct,maxtf=6):
    r=R[tis]; s=r[(r.cell_type==ct)&(r.p<0.05)&(r.n_linked>=3)].sort_values('p').head(maxtf)
    tfs=list(s.TF)
    d=e[(e.tissue==tis)&(e.cell_type==ct)&(e.TF.isin(tfs))]
    cnt=d.groupby('gene').TF.nunique()
    genes=set(cnt[cnt>=2].index)
    if len(genes)<5: genes=set(cnt.index)
    for tf in tfs:                       # guarantee every module TF has at least one shown target
        sub=d[d.TF==tf]
        if not set(sub.gene)&genes:
            genes|=set(sub.sort_values('r',ascending=False).gene.head(2))
    genes=sorted(genes)
    gct=df[df.tissue==tis].groupby('gene').cell_type.nunique()
    excl={g:(gct[g]==1) for g in genes}
    pairs=d[d.gene.isin(genes)][['TF','gene','r']]
    return tfs,genes,excl,pairs,s.set_index('TF')

def draw(ax,tis,ct,full):
    tfs,genes,excl,pairs,S=module(tis,ct)
    ny=max(len(tfs),len(genes))
    ty=np.linspace(0,1,len(tfs)+2)[1:-1] if len(tfs)>1 else np.array([0.5])
    ty=0.5+(ty-0.5)*min(1.0,len(tfs)/ny*1.35)
    gy=np.linspace(0,1,len(genes)+2)[1:-1]
    TP={t:(0.0,y) for t,y in zip(tfs,ty)}
    GP={g:(1.0,y) for g,y in zip(genes,gy)}
    for _,p in pairs.sort_values('r').iterrows():
        x0,y0=TP[p.TF]; x1,y1=GP[p.gene]
        c=[(x0,y0),(0.5,y0),(0.5,y1),(x1,y1)]
        ax.add_patch(PathPatch(Path(c,[Path.MOVETO,Path.CURVE4,Path.CURVE4,Path.CURVE4]),fill=False,
            lw=0.55 if p.r>0.7 else 0.35,color='#27343f' if p.r>0.7 else '#9aa3ab',
            alpha=0.55 if p.r>0.7 else 0.30,zorder=2))
    for t in tfs:
        x,y=TP[t]
        n=int(S.loc[t,'n_linked'])
        ax.scatter([x],[y],s=40,c=[TFC2 if S.loc[t,'q']<0.05 else TFC],edgecolors='white',linewidths=0.6,zorder=4,marker='s')
        ax.text(x-0.045,y,t,ha='right',va='center',fontsize=5.6,color=INK,fontweight='bold',zorder=5)
    for g in genes:
        x,y=GP[g]
        ax.scatter([x],[y],s=16,c=[GEX if excl[g] else GSH],edgecolors='#4f5a65',linewidths=0.45,zorder=4)
        ax.text(x+0.042,y,g,ha='left',va='center',fontsize=4.9,color=INK,zorder=5)
    ax.set_xlim(-0.60,1.62); ax.set_ylim(-0.06,1.06); ax.axis('off')
    ax.set_title(f'{full} ({ct}) · {tis.lower()}',fontsize=6.9,color=INK,fontweight='bold',pad=3)

fig=plt.figure(figsize=(7.09,6.4))
gs=fig.add_gridspec(2,3,wspace=0.06,hspace=0.16,left=0.02,right=0.985,top=0.945,bottom=0.075)
for k,(tis,ct,full) in enumerate(MODS):
    ax=fig.add_subplot(gs[k//3,k%3]); draw(ax,tis,ct,full)
    ax.text(-0.02,1.055,'abcdef'[k],transform=ax.transAxes,fontsize=9.5,fontweight='bold',va='bottom')

from matplotlib.lines import Line2D
h=[Line2D([],[],marker='s',ls='',mfc=TFC2,mec='white',ms=5,label='Cell-type-enriched TF (FDR < 5%)'),
   Line2D([],[],marker='s',ls='',mfc=TFC,mec='white',ms=5,label='Cell-type-enriched TF (P < 0.05)'),
   Line2D([],[],marker='o',ls='',mfc=GEX,mec='#4f5a65',mew=0.45,ms=4,label='Risk gene unique to this cell type'),
   Line2D([],[],marker='o',ls='',mfc='white',mec='#4f5a65',mew=0.45,ms=4,label='Risk gene shared with other cell types'),
   Line2D([],[],color='#27343f',lw=0.8,alpha=0.6,label='$r$ > 0.7'),
   Line2D([],[],color='#9aa3ab',lw=0.6,alpha=0.6,label='0.5 < $r$ ≤ 0.7')]
fig.legend(handles=h,loc='lower center',ncol=3,frameon=False,fontsize=5.9,
           handletextpad=0.35,columnspacing=1.6,bbox_to_anchor=(0.5,0.005))
fig.savefig('Fig3_celltype_modules.pdf'); fig.savefig('Fig3_celltype_modules.png',dpi=230)
print('ok')
