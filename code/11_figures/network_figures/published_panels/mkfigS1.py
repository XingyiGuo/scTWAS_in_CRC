import numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
import matplotlib.cm as cm
from base import C,E,OCC
plt.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Liberation Sans','Arial','DejaVu Sans'],'pdf.fonttype':42,'axes.linewidth':0.5})
INK='#1a1a1a'; INK2='#52514e'
GRN=LinearSegmentedColormap.from_list('o',['#f4f9f7','#b9e3d4','#5ec6a4','#1baf7a','#0d6d4b'])
BLU=LinearSegmentedColormap.from_list('c',['#f7f8f6','#cddff0','#7fb0e0','#2a78d6','#123f72'])
ORD={'Normal':['STM','TAC','ABS','GOB','EE','TUF','CT'],
     'Adenoma':['STM','TAC','ABS','GOB','EE','TUF','CT','ASC','SSC']}
tfs=sorted(C.TF.unique())
def mats(tis):
    cts=ORD[tis]
    co=E[E.tissue==tis].groupby(['TF','cell_type']).gene.nunique().unstack().reindex(index=tfs,columns=cts).fillna(0)
    oc=OCC[OCC.tissue==tis].pivot(index='TF',columns='cell_type',values='n_occupied').reindex(index=tfs,columns=cts).fillna(0)
    return co,oc
cN,oN=mats('Normal'); cA,oA=mats('Adenoma')
order=(cN.sum(axis=1).rank()+cA.sum(axis=1).rank()+oN.sum(axis=1).rank()*0.6).sort_values(ascending=False).index
PAN=[(cN,'Risk genes co-expressed\nwith TF — normal',BLU,'a'),
     (oN,'Risk genes whose lead variant\nthe TF occupies — normal',GRN,'b'),
     (cA,'Risk genes co-expressed\nwith TF — adenoma',BLU,'c'),
     (oA,'Risk genes whose lead variant\nthe TF occupies — adenoma',GRN,'d')]
fig=plt.figure(figsize=(7.09,7.6))
gs=fig.add_gridspec(1,4,width_ratios=[7,7,9,9],wspace=0.10,left=0.105,right=0.982,top=0.885,bottom=0.075)
for k,(M,ttl,cmap,lab) in enumerate(PAN):
    M=M.loc[order]; ax=fig.add_subplot(gs[0,k]); vmax=max(M.values.max(),1)
    im=ax.imshow(M.values,cmap=cmap,vmin=0,vmax=vmax,aspect='auto')
    ax.set_xticks(range(M.shape[1])); ax.set_xticklabels(M.columns,fontsize=5.8,rotation=90)
    if k==0: ax.set_yticks(range(len(M))); ax.set_yticklabels(M.index,fontsize=5.4)
    else: ax.set_yticks([])
    ax.set_title(ttl,fontsize=6.3,color=INK,pad=4)
    ax.text(-0.02 if k else -0.42,1.055,lab,transform=ax.transAxes,fontsize=9.5,fontweight='bold',va='bottom')
    ax.set_xticks(np.arange(-.5,M.shape[1],1),minor=True); ax.set_yticks(np.arange(-.5,len(M),1),minor=True)
    ax.grid(which='minor',color='white',lw=0.5); ax.tick_params(which='minor',length=0); ax.tick_params(length=0)
    for s in ax.spines.values(): s.set_color('#c9c8c3')
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            v=int(M.values[i,j])
            if v>0: ax.text(j,i,str(v),ha='center',va='center',fontsize=3.8,color=('white' if v/vmax>0.55 else INK))
    cb=fig.colorbar(im,ax=ax,orientation='horizontal',fraction=0.03,pad=0.05,aspect=20)
    cb.ax.tick_params(labelsize=4.6,length=2,width=0.4); cb.outline.set_linewidth(0.4)
    cb.set_label('Number of risk genes',fontsize=4.9,labelpad=1,color=INK2)
fig.savefig('ExtendedDataFigure7_TF_cell_type_matrix.pdf'); fig.savefig('ExtendedDataFigure7_TF_cell_type_matrix.png',dpi=220)
print('ok')
