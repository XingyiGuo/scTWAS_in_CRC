import numpy as np, itertools, json, sys
from PIL import Image
Image.MAX_IMAGE_PIXELS=None
BLUE=(76,120,168); RED=(228,87,86)

def runs(mask1d):
    idx=np.where(mask1d)[0]
    out=[]
    if len(idx)==0: return out
    s=idx[0]; p=idx[0]
    for i in idx[1:]:
        if i>p+1: out.append((s,p)); s=i
        p=i
    out.append((s,p)); return out

def panel_bars(a, y0,y1, x0,x1, legend_h=0.10):
    sub=a[y0:y1, x0:x1]
    red=np.all(sub==RED,axis=2); blue=np.all(sub==BLUE,axis=2)
    # remove legend swatches: small squares in top region
    H=sub.shape[0]
    cut=int(H*legend_h)
    red[:cut]=False; blue[:cut]=False
    any_=red|blue
    if any_.sum()==0: return None
    baseline=np.where(any_.any(axis=1))[0].max()
    bars=[]
    for name,m in (('red',red),('blue',blue)):
        colhas=m.any(axis=0)
        for s,e in runs(colhas):
            if e-s < 3: continue
            seg=m[:, s:e+1]
            top=np.where(seg.any(axis=1))[0].min()
            bot=np.where(seg.any(axis=1))[0].max()
            if bot < baseline-8: continue   # not a bar anchored at axis
            bars.append(dict(kind=name, xs=s+x0, xe=e+x0, xc=(s+e)/2+x0, h=baseline-top+1))
    bars.sort(key=lambda b:b['xc'])
    return dict(baseline=baseline+y0, bars=bars)

def fit_scale(hs):
    best=None
    hmin=min(hs)
    for s in np.arange(hmin*0.45, hmin*1.6, 0.02):
        v=np.array(hs)/s
        err=np.sum((v-np.round(v))**2)
        if np.any(np.round(v)<1): continue
        if best is None or err<best[1]: best=(s,err)
    return best
