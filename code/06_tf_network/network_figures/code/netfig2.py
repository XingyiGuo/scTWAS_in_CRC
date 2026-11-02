import numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.path import Path
from matplotlib.patches import PathPatch, Wedge
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.patheffects as pe

plt.rcParams.update({'font.family':'Liberation Sans','pdf.fonttype':42,'ps.fonttype':42,
 'axes.linewidth':0.6,'xtick.major.width':0.6,'ytick.major.width':0.6,
 'xtick.labelsize':6,'ytick.labelsize':6,'axes.labelsize':7,'legend.fontsize':6})
BLUE='#2a78d6'; ORANGE='#eb6834'; AQUA='#1baf7a'
INK='#1a1a1a'; INK2='#52514e'
TFRAMP=LinearSegmentedColormap.from_list('occ',['#eef3fa','#a9c8ea','#4a8ed6','#12406f'])

def key_tfs(df,occ,ntop=20):
    o=occ.groupby('TF').n_occupied.sum(); out={}
    for tis in ['Normal','Adenoma']:
        e=df[(df.tissue==tis)&(df.r>0.5)]
        deg=e.groupby('TF').gene.nunique()
        m=pd.DataFrame({'deg':deg}).join(o.rename('occ')).fillna(0)
        m['sc']=(m.deg.rank(pct=True)+m.occ.rank(pct=True))/2
        out[tis]=list(m.sort_values('sc',ascending=False).head(ntop).index)
    return out

def sectors_layout(genes,prim,order,gap_deg=5.0):
    groups={}
    for g in genes: groups.setdefault(prim[g],[]).append(g)
    keys=[k for k in order if k in groups]+(['Multiple'] if 'Multiple' in groups else [])
    n=len(genes); gap=np.deg2rad(gap_deg)
    total=2*np.pi-gap*len(keys)
    pos={}; secs=[]; a=np.pi/2
    for k in keys:
        gs=sorted(groups[k]); w=total*len(gs)/n
        for i,g in enumerate(gs): pos[g]=a-w*(i+0.5)/len(gs)
        secs.append((k,a,a-w)); a=a-w-gap
    return pos,secs

def bez(ax,p0,p1,curv,**kw):
    c=(np.array(p0)+np.array(p1))/2*(1-curv)
    ax.add_patch(PathPatch(Path([p0,tuple(c),p1],[Path.MOVETO,Path.CURVE3,Path.CURVE3]),fill=False,**kw))

def rot_ha(a):
    r=np.rad2deg(a); ha='left'
    if r>90 or r<-90: r+=180; ha='right'
    return r,ha

def draw(ax,pairs,prim,order,occ_sum,occ_vmax,title,sub):
    genes=sorted(pairs.gene.unique()); tfs=sorted(pairs.TF.unique())
    deg_g=pairs.groupby('gene').size(); deg_t=pairs.groupby('TF').size()
    ang,secs=sectors_layout(genes,prim,order)
    gp={g:(np.cos(ang[g]),np.sin(ang[g])) for g in genes}
    tang={}
    for tf,grp in pairs.groupby('TF'):
        w=grp.rmax.values; a=np.array([ang[g] for g in grp.gene])
        tang[tf]=np.arctan2((w*np.sin(a)).sum(),(w*np.cos(a)).sum())
    dmax=deg_t.max()
    P=np.array([[ (0.60-0.36*(deg_t[t]/dmax)**0.55)*np.cos(tang[t]),
                  (0.60-0.36*(deg_t[t]/dmax)**0.55)*np.sin(tang[t])] for t in tfs])
    for _ in range(400):
        D=P[:,None,:]-P[None,:,:]; d=np.linalg.norm(D,axis=2)+1e-9
        np.fill_diagonal(d,1e9)
        f=np.where(d<0.245,(0.245-d)/d,0)[:,:,None]*D
        P=P+0.22*f.sum(1)
        rr=np.linalg.norm(P,axis=1); ov=rr>0.70
        P[ov]=P[ov]/rr[ov,None]*0.70
    tp={t:tuple(P[i]) for i,t in enumerate(tfs)}
    # sector bands + radial separators
    for i,(k,a0,a1) in enumerate(secs):
        ax.add_patch(Wedge((0,0),1.045,np.rad2deg(a1),np.rad2deg(a0),width=0.045,
            facecolor=('#d3d2cc' if i%2==0 else '#eceae5'),edgecolor='none',zorder=1))
        am=(a0+a1)/2; rot,ha=rot_ha(am)
        ax.text(1.38*np.cos(am),1.38*np.sin(am),k,rotation=rot,rotation_mode='anchor',
                ha=ha,va='center',fontsize=6.2,color=INK,zorder=8,fontweight='bold')
        ax.plot([1.05*np.cos(am),1.36*np.cos(am)],[1.05*np.sin(am),1.36*np.sin(am)],
                lw=0.35,color='#c7c6c0',zorder=1)
    for _,e in pairs.sort_values('rmax').iterrows():
        strong=e.rmax>0.7
        bez(ax,tp[e.TF],gp[e.gene],0.52,lw=0.32 if strong else 0.18,
            color='#2b3947' if strong else '#a4acb4',alpha=0.45 if strong else 0.20,zorder=2)
    ax.scatter([gp[g][0] for g in genes],[gp[g][1] for g in genes],
               s=[5+2.6*deg_g[g] for g in genes],c='white',edgecolors='#4f5a65',linewidths=0.4,zorder=4)
    for g in genes:
        rot,ha=rot_ha(ang[g])
        ax.text(1.065*np.cos(ang[g]),1.065*np.sin(ang[g]),g,rotation=rot,rotation_mode='anchor',
                ha=ha,va='center',fontsize=4.1,color=INK,zorder=7)
    for t in tfs:
        x,y=tp[t]
        ax.scatter([x],[y],s=34+72*(deg_t[t]/dmax)**0.7,c=[TFRAMP(occ_sum.get(t,0)/occ_vmax)],
                   edgecolors='#0f2f4d',linewidths=0.5,zorder=5)
        ax.text(x,y-0.075,t,ha='center',va='top',fontsize=4.7,color=INK,zorder=8,fontweight='bold',
                path_effects=[pe.withStroke(linewidth=1.4,foreground='white')])
    ax.set_xlim(-1.60,1.60); ax.set_ylim(-1.60,1.74); ax.set_aspect('equal'); ax.axis('off')
    ax.text(0,1.73,title,ha='center',va='top',fontsize=7.6,color=INK,fontweight='bold')
    ax.text(0,1.645,sub,ha='center',va='top',fontsize=6.0,color=INK2)
