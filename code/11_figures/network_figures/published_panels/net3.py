import numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.path import Path
from matplotlib.patches import PathPatch, Wedge
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.patheffects as pe
plt.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Liberation Sans','Arial','DejaVu Sans'],'pdf.fonttype':42,'ps.fonttype':42,
 'axes.linewidth':0.6,'xtick.labelsize':6,'ytick.labelsize':6,'axes.labelsize':7,'legend.fontsize':6})
BLUE='#2a78d6'; ORANGE='#eb6834'; AQUA='#1baf7a'
INK='#1a1a1a'; INK2='#52514e'
TFRAMP=LinearSegmentedColormap.from_list('occ',['#eef3fa','#a9c8ea','#4a8ed6','#12406f'])

def sectors_layout(genes,prim,order,gap_deg=6.0):
    groups={}
    for g in genes: groups.setdefault(prim[g],[]).append(g)
    keys=[k for k in order if k in groups]+(['Multiple'] if 'Multiple' in groups else [])
    n=len(genes); gap=np.deg2rad(gap_deg); total=2*np.pi-gap*len(keys)
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

def draw(ax,pairs,prim,order,occ_sum,occ_vmax,title,sub,gene_fs=4.6,tf_fs=5.0):
    genes=sorted(pairs.gene.unique()); tfs=sorted(pairs.TF.unique())
    deg_g=pairs.groupby('gene').TF.nunique(); deg_t=pairs.groupby('TF').gene.nunique()
    ang,secs=sectors_layout(genes,prim,order)
    gp={g:(np.cos(ang[g]),np.sin(ang[g])) for g in genes}
    tang={}
    for tf,grp in pairs.groupby('TF'):
        w=grp.rmax.values; a=np.array([ang[g] for g in grp.gene])
        tang[tf]=np.arctan2((w*np.sin(a)).sum(),(w*np.cos(a)).sum())
    dmax=deg_t.max()
    P=np.array([[ (0.55-0.34*(deg_t[t]/dmax)**0.55)*np.cos(tang[t]),
                  (0.55-0.34*(deg_t[t]/dmax)**0.55)*np.sin(tang[t])] for t in tfs])
    # box-based de-overlap: node + its label below it
    HW=np.array([0.0275*len(t)+0.034 for t in tfs])
    HH=np.full(len(tfs),0.085)
    for it in range(900):
        moved=0.0
        for i in range(len(tfs)):
            for j in range(i+1,len(tfs)):
                dx=P[j,0]-P[i,0]; dy=P[j,1]-P[i,1]
                ox=HW[i]+HW[j]-abs(dx); oy=HH[i]+HH[j]-abs(dy)
                if ox>0 and oy>0:
                    if ox/(HW[i]+HW[j]) < oy/(HH[i]+HH[j]):
                        sh=(ox/2+0.002)*(1 if dx>=0 else -1)
                        P[j,0]+=sh; P[i,0]-=sh; moved+=abs(sh)
                    else:
                        sh=(oy/2+0.002)*(1 if dy>=0 else -1)
                        P[j,1]+=sh; P[i,1]-=sh; moved+=abs(sh)
        rr=np.linalg.norm(P,axis=1); ov=rr>0.635
        P[ov]=P[ov]/rr[ov,None]*0.635
        if moved<1e-4: break
    tp={t:tuple(P[i]) for i,t in enumerate(tfs)}
    for i,(k,a0,a1) in enumerate(secs):
        ax.add_patch(Wedge((0,0),1.045,np.rad2deg(a1),np.rad2deg(a0),width=0.045,
            facecolor=('#d3d2cc' if i%2==0 else '#eceae5'),edgecolor='none',zorder=1))
        am=(a0+a1)/2; rot,ha=rot_ha(am)
        ax.text(1.40*np.cos(am),1.40*np.sin(am),k,rotation=rot,rotation_mode='anchor',
                ha=ha,va='center',fontsize=6.4,color=INK,zorder=8,fontweight='bold')
        ax.plot([1.05*np.cos(am),1.37*np.cos(am)],[1.05*np.sin(am),1.37*np.sin(am)],
                lw=0.35,color='#c7c6c0',zorder=1)
    for _,e in pairs.sort_values('occ').iterrows():
        occ=bool(e.occ)
        w=0.26+0.55*max(min((e.rmax-0.5)/0.4,1.0),0.0)
        bez(ax,tp[e.TF],gp[e.gene],0.52,lw=w,
            color='#1b4f72' if occ else '#a7aeb5',alpha=0.70 if occ else 0.28,zorder=3 if occ else 2)
    ax.scatter([gp[g][0] for g in genes],[gp[g][1] for g in genes],
               s=[10+7.0*deg_g[g] for g in genes],c='white',edgecolors='#4f5a65',linewidths=0.45,zorder=4)
    for g in genes:
        rot,ha=rot_ha(ang[g])
        ax.text(1.07*np.cos(ang[g]),1.07*np.sin(ang[g]),g,rotation=rot,rotation_mode='anchor',
                ha=ha,va='center',fontsize=gene_fs,color=INK,zorder=7)
    for t in tfs:
        x,y=tp[t]
        ax.scatter([x],[y],s=30+80*(deg_t[t]/dmax)**0.7,c=[TFRAMP(min(occ_sum.get(t,0)/occ_vmax,1))],
                   edgecolors='#0f2f4d',linewidths=0.5,zorder=5)
        ax.text(x,y-0.080,t,ha='center',va='top',fontsize=tf_fs,color=INK,zorder=8,fontweight='bold',
                path_effects=[pe.withStroke(linewidth=1.5,foreground='white')])
    ax.set_xlim(-1.66,1.66); ax.set_ylim(-1.66,1.96); ax.set_aspect('equal'); ax.axis('off')
    ax.text(0,1.95,title,ha='center',va='top',fontsize=7.8,color=INK,fontweight='bold')
    ax.text(0,1.855,sub,ha='center',va='top',fontsize=6.1,color=INK2)
