import re,os,glob,json
import numpy as np, pandas as pd
from extract2 import parse

def matrix(pdf):
    rects,texts=parse(pdf)
    from collections import Counter
    cnt=Counter((round(r[2],2),round(r[3],2)) for r in rects)
    (cw,ch),_=cnt.most_common(1)[0]
    cells=[r for r in rects if round(r[2],2)==cw and round(r[3],2)==ch]
    legend=sorted([r for r in rects if (round(r[2],2),round(r[3],2))!=(cw,ch)],key=lambda r:r[1])
    lcol=np.array([r[4] for r in legend]); n=len(legend)
    lx=np.median([r[0] for r in legend]); lw=np.median([r[2] for r in legend])
    ticks=[]
    for t in texts:
        s=t['text'].replace('−','-')
        if re.fullmatch(r'-?\d+(\.\d+)?',s) and lx+lw-2 < t['x'] < lx+lw+40 and t['b']==0:
            ticks.append((t['y'],float(s)))
    ticks=sorted(set(ticks))
    vmin,vmax=min(v for _,v in ticks),max(v for _,v in ticks)
    lval=vmin+(np.arange(n)+0.5)*(vmax-vmin)/n
    xs=sorted(set(round(r[0],2) for r in cells)); ys=sorted(set(round(r[1],2) for r in cells))
    xi={v:i for i,v in enumerate(xs)}; yi={v:i for i,v in enumerate(ys)}
    M=np.full((len(ys),len(xs)),np.nan); worst=0
    for x,y,w,h,c in cells:
        d=np.sum((lcol-np.array(c))**2,axis=1); j=int(np.argmin(d))
        worst=max(worst,d[j])
        if d[j]>0.01: continue
        M[yi[round(y,2)],xi[round(x,2)]]=lval[j]
    rowtx=sorted([t for t in texts if t['b']==0 and abs(t['a']-6.0)<0.01 and t['x']>max(xs)+cw-2],key=lambda t:t['y'])
    coltx=sorted([t for t in texts if abs(t['b']-t['a'])<0.01 and t['a']<6],key=lambda t:t['x'])
    rows=[t['text'] for t in rowtx]; cols=[t['text'] for t in coltx]
    title=[t['text'] for t in texts if t['a']>=12]
    ok=(len(rows)==len(ys)) and (len(cols)==len(xs))
    return dict(M=M,rows=rows,cols=cols,nrow=len(ys),ncol=len(xs),title=title,
                ticks=ticks,ok=ok,worst=float(worst),nlegend=n)

out={}
for p in sorted(glob.glob('/mnt/user-data/uploads/TFgenes/*.pdf')):
    d=matrix(p); b=os.path.basename(p)
    print(f"{b:60s} {d['nrow']}x{d['ncol']} labels ok={d['ok']} (rows {len(d['rows'])}, cols {len(d['cols'])}) worstcol={d['worst']:.4f} ticks={[v for _,v in d['ticks']]}")
    out[b]=d
np.save('mats.npy',out,allow_pickle=True)
