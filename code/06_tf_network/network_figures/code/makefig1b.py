import numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from netfig2 import *
from matplotlib.colors import Normalize
import matplotlib.cm as cm

df=pd.read_csv('tf_gene_correlations.csv'); occ=pd.read_csv('tf_leadvariant_occupancy.csv')
K=key_tfs(df,occ,20); KEY=sorted(set(K['Normal'])|set(K['Adenoma']))
pd.Series(KEY).to_csv('key_tfs.csv',index=False,header=['TF'])
ORD_N=['STM','TAC','ABS','GOB','EE','TUF','CT']
ORD_A=['STM','TAC','ABS','GOB','EE','TUF','CT','ASC','SSC']
o_all=occ.groupby('TF').n_occupied.sum(); VMAX=o_all.max()

def prep(tis):
    d=df[(df.tissue==tis)&(df.r>0.5)&(df.TF.isin(KEY))]
    pairs=d.groupby(['TF','gene']).agg(rmax=('r','max'),nct=('cell_type','nunique')).reset_index()
    gct=df[df.tissue==tis].groupby('gene').cell_type.apply(lambda s:sorted(set(s)))
    prim={g:(gct[g][0] if len(gct[g])==1 else 'Multiple') for g in pairs.gene.unique()}
    return pairs,prim

fig=plt.figure(figsize=(7.09,7.15))
gs=fig.add_gridspec(1,2,wspace=0.015,left=0.005,right=0.995,top=0.995,bottom=0.475)
gs2=fig.add_gridspec(1,2,wspace=0.36,left=0.085,right=0.985,top=0.415,bottom=0.065)
info={}
for i,(tis,order) in enumerate([('Normal',ORD_N),('Adenoma',ORD_A)]):
    ax=fig.add_subplot(gs[0,i]); pairs,prim=prep(tis); info[tis]=pairs
    nm='Normal colon epithelium' if tis=='Normal' else 'Adenoma and serrated lesions'
    draw(ax,pairs,prim,order,o_all,VMAX,nm,
         f'{len(pairs)} TF–risk gene links · {pairs.TF.nunique()} key TFs · {pairs.gene.nunique()} risk genes')
    ax.text(0.005,0.995,'ab'[i],transform=ax.transAxes,fontsize=10,fontweight='bold',va='top')

cax=fig.add_axes([0.435,0.487,0.145,0.0075])
cb=fig.colorbar(cm.ScalarMappable(cmap=TFRAMP,norm=Normalize(0,VMAX)),cax=cax,orientation='horizontal')
cb.set_label('TF fill: risk-gene lead variants occupied',fontsize=5.3,labelpad=1.5,color=INK2)
cb.ax.tick_params(labelsize=5,length=2,width=0.5); cb.outline.set_linewidth(0.4)

# ---- c
axc=fig.add_subplot(gs2[0,0])
for tis,col,mk in [('Normal',BLUE,'o'),('Adenoma',ORANGE,'^')]:
    dt=df[(df.tissue==tis)&(df.r>0.5)].groupby('TF').gene.nunique()
    axc.scatter([o_all.get(t,0) for t in dt.index],dt.values,s=15,facecolor=col,
                edgecolor='white',linewidth=0.35,alpha=0.9,marker=mk,label=tis,zorder=3)
seen=set()
for t in ['POLR2A','RAD21','SP1','MAX','CBX3','GATA6','SMC3','YY1','CEBPB','FOXA2','TFAP4','ELF1']:
    dn=df[(df.tissue=='Normal')&(df.r>0.5)].groupby('TF').gene.nunique()
    if t in dn.index:
        axc.annotate(t,(o_all.get(t,0),dn[t]),textcoords='offset points',xytext=(3.5,2.5),
                     fontsize=4.9,color=INK2)
axc.set_xlabel('Risk-gene lead variants occupied by TF')
axc.set_ylabel('Risk genes co-expressed with TF ($r$ > 0.5)')
axc.spines[['top','right']].set_visible(False)
axc.grid(axis='y',lw=0.4,color='#eceae5',zorder=0); axc.set_axisbelow(True)
axc.legend(frameon=False,loc='upper right',handletextpad=0.3)
axc.text(-0.185,1.05,'c',transform=axc.transAxes,fontsize=10,fontweight='bold',va='bottom')

# ---- d
axd=fig.add_subplot(gs2[0,1])
dn=df[(df.tissue=='Normal')&(df.r>0.5)].groupby('TF').gene.nunique()
da=df[(df.tissue=='Adenoma')&(df.r>0.5)].groupby('TF').gene.nunique()
m=pd.DataFrame({'Normal':dn,'Adenoma':da,'occ':o_all}).fillna(0).loc[KEY]
m['sc']=m.Normal.rank()+m.Adenoma.rank()+m.occ.rank()
m=m.sort_values('sc')
y=np.arange(len(m))
axd.hlines(y,0,m[['Normal','Adenoma','occ']].max(axis=1),color='#eceae5',lw=0.7,zorder=1)
axd.scatter(m.occ,y,s=19,color=AQUA,marker='s',label='Lead variants occupied',zorder=3,edgecolor='white',linewidth=0.35)
axd.scatter(m.Normal,y,s=17,color=BLUE,label='Co-expressed, normal',zorder=4,edgecolor='white',linewidth=0.35)
axd.scatter(m.Adenoma,y,s=17,color=ORANGE,marker='^',label='Co-expressed, adenoma',zorder=5,edgecolor='white',linewidth=0.35)
axd.set_yticks(y); axd.set_yticklabels(m.index,fontsize=5.4)
axd.set_ylim(-0.8,len(m)-0.2)
axd.set_xlabel('Number of cell-type risk genes')
axd.spines[['top','right','left']].set_visible(False)
axd.tick_params(axis='y',length=0)
axd.grid(axis='x',lw=0.4,color='#f3f2ee',zorder=0); axd.set_axisbelow(True)
axd.legend(frameon=False,loc='lower right',handletextpad=0.3,borderpad=0.2)
axd.text(-0.25,1.05,'d',transform=axd.transAxes,fontsize=10,fontweight='bold',va='bottom')

fig.savefig('Fig_TFnetwork.pdf'); fig.savefig('Fig_TFnetwork.png',dpi=230)
print('key TFs',len(KEY),KEY)
