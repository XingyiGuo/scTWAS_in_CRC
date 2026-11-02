import numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
import matplotlib.cm as cm
plt.rcParams.update({'font.family':'Liberation Sans','pdf.fonttype':42,'axes.linewidth':0.5})
INK='#1a1a1a'; INK2='#52514e'
SEQ=LinearSegmentedColormap.from_list('s',['#f7f8f6','#cddff0','#7fb0e0','#2a78d6','#123f72'])
OCC=LinearSegmentedColormap.from_list('o',['#f4f9f7','#b9e3d4','#5ec6a4','#1baf7a','#0d6d4b'])

df=pd.read_csv('tf_gene_correlations.csv'); occ=pd.read_csv('tf_leadvariant_occupancy.csv')
e=df[df.r>0.5]
N=e[e.tissue=='Normal'].groupby(['TF','cell_type']).gene.nunique().unstack(fill_value=0)
A=e[e.tissue=='Adenoma'].groupby(['TF','cell_type']).gene.nunique().unstack(fill_value=0)
O=occ.pivot(index='TF',columns='cell_type',values='n_occupied').fillna(0)
tfs=sorted(set(df.TF))
ordn=['STM','TAC','ABS','GOB','EE','TUF','CT']; orda=ordn+['ASC','SSC']
N=N.reindex(index=tfs,columns=ordn).fillna(0); A=A.reindex(index=tfs,columns=orda).fillna(0)
O=O.reindex(index=tfs,columns=orda).fillna(0)
order=(N.sum(axis=1).rank()+A.sum(axis=1).rank()+O.sum(axis=1).rank()).sort_values(ascending=False).index
N,A,O=N.loc[order],A.loc[order],O.loc[order]

fig=plt.figure(figsize=(7.09,7.4))
gs=fig.add_gridspec(1,3,width_ratios=[7,9,9],wspace=0.10,left=0.115,right=0.975,top=0.90,bottom=0.075)
vmaxC=max(N.values.max(),A.values.max())
for k,(M,ttl,cmp_,vmax,lab) in enumerate([
      (N,'Co-expressed risk genes\nnormal colon',SEQ,vmaxC,'a'),
      (A,'Co-expressed risk genes\nadenoma / serrated',SEQ,vmaxC,'b'),
      (O,'Risk-gene lead variants\noccupied by TF',OCC,O.values.max(),'c')]):
    ax=fig.add_subplot(gs[0,k])
    im=ax.imshow(M.values,cmap=cmp_,vmin=0,vmax=vmax,aspect='auto')
    ax.set_xticks(range(M.shape[1])); ax.set_xticklabels(M.columns,fontsize=6,rotation=90)
    if k==0:
        ax.set_yticks(range(len(M))); ax.set_yticklabels(M.index,fontsize=5.6)
    else: ax.set_yticks([])
    ax.set_title(ttl,fontsize=7,color=INK,pad=4)
    ax.text(-0.02 if k else -0.30,1.035,lab,transform=ax.transAxes,fontsize=10,fontweight='bold',va='bottom')
    for s in ax.spines.values(): s.set_color('#c9c8c3')
    ax.set_xticks(np.arange(-.5,M.shape[1],1),minor=True); ax.set_yticks(np.arange(-.5,len(M),1),minor=True)
    ax.grid(which='minor',color='white',lw=0.5); ax.tick_params(which='minor',length=0)
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            v=int(M.values[i,j])
            if v>0: ax.text(j,i,str(v),ha='center',va='center',fontsize=4.0,
                            color=('white' if v/vmax>0.55 else INK))
    cax=fig.inset_axes if False else None
    cb=fig.colorbar(im,ax=ax,orientation='horizontal',fraction=0.03,pad=0.055,aspect=26)
    cb.ax.tick_params(labelsize=5,length=2,width=0.4); cb.outline.set_linewidth(0.4)
    cb.set_label('Number of risk genes',fontsize=5.4,labelpad=1,color=INK2)
fig.savefig('FigS_TFmatrix.pdf'); fig.savefig('FigS_TFmatrix.png',dpi=220)
print('ok')
