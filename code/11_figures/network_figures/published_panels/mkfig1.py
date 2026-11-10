import numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from net3 import *
from matplotlib.colors import Normalize
import matplotlib.cm as cm
from matplotlib.lines import Line2D
from base import E,OCC_T

KEY=list(pd.read_csv('key_tfs.csv').TF)
ORD_N=['STM','TAC','ABS','GOB','EE','TUF','CT']
ORD_A=['STM','TAC','ABS','GOB','EE','TUF','CT','ASC','SSC']
OCCT={t:OCC_T[OCC_T.tissue==t].set_index('TF').occ_total for t in ['Normal','Adenoma']}
VMAX=max(OCCT['Normal'].max(),OCCT['Adenoma'].max())

def prep(tis):
    d=E[(E.tissue==tis)&(E.TF.isin(KEY))]
    pairs=d.groupby(['TF','gene']).agg(rmax=('r','max'),nct=('cell_type','nunique'),
                                       occ=('occupies_lead_variant','max')).reset_index()
    gct=d.groupby('gene').cell_type.apply(lambda s:sorted(set(s)))
    prim={g:(gct[g][0] if len(gct[g])==1 else 'Multiple') for g in pairs.gene.unique()}
    return pairs,prim

fig=plt.figure(figsize=(7.09,4.15))
gs=fig.add_gridspec(1,2,wspace=0.015,left=0.005,right=0.995,top=0.985,bottom=0.12)
for i,(tis,order) in enumerate([('Normal',ORD_N),('Adenoma',ORD_A)]):
    ax=fig.add_subplot(gs[0,i]); pairs,prim=prep(tis)
    nm='Normal colon epithelium' if tis=='Normal' else 'Adenoma and serrated lesions'
    draw(ax,pairs,prim,order,OCCT[tis],VMAX,nm,
         f'{len(pairs)} TF–risk gene links · {pairs.TF.nunique()} key TFs · {pairs.gene.nunique()} risk genes',
         gene_fs=4.2,tf_fs=4.9)
    ax.text(0.005,0.995,'bc'[i],transform=ax.transAxes,fontsize=10,fontweight='bold',va='top')

cax=fig.add_axes([0.345,0.072,0.135,0.012])
cb=fig.colorbar(cm.ScalarMappable(cmap=TFRAMP,norm=Normalize(0,VMAX)),cax=cax,orientation='horizontal')
cb.set_label('TF fill: risk-gene lead variants occupied',fontsize=5.3,labelpad=1.4,color=INK2)
cb.ax.tick_params(labelsize=5,length=2,width=0.4); cb.outline.set_linewidth(0.4)
h=[Line2D([],[],color='#1b4f72',lw=1.0,alpha=0.75,label='TF also occupies the gene’s lead variant'),
   Line2D([],[],color='#a7aeb5',lw=0.6,alpha=0.5,label='Co-expression only')]
fig.legend(handles=h,loc='center left',bbox_to_anchor=(0.545,0.074),frameon=False,fontsize=5.3,
           handlelength=1.6,handletextpad=0.4)

fig.savefig('Figure4bc_key_TF_network.pdf')
fig.savefig('Figure4bc_key_TF_network.png',dpi=300)
