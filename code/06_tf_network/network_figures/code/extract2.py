import re,subprocess,os,glob,json
import numpy as np

def content(pdf):
    out=subprocess.run(['qpdf','--qdf','--object-streams=disable',pdf,'-'],capture_output=True).stdout
    parts=sorted(re.findall(rb'stream\r?\n(.*?)\r?\nendstream',out,re.S),key=len,reverse=True)
    return parts[0].decode('latin-1')

def unesc(s):
    return s.replace('\\(','(').replace('\\)',')').replace('\\\\','\\')

def parse(pdf):
    s=content(pdf)
    rects=[];fill=None
    for m in re.finditer(r'([\d.]+) ([\d.]+) ([\d.]+) scn|([-\d.]+) ([-\d.]+) ([-\d.]+) ([-\d.]+) re\s+f',s):
        if m.group(1) is not None:
            fill=tuple(float(m.group(i)) for i in (1,2,3))
        else:
            rects.append(tuple(float(m.group(i)) for i in (4,5,6,7))+(fill,))
    texts=[]
    for b in re.findall(r'BT\n(.*?)\nET\n',s,re.S):
        m=re.search(r'/F\d+ 1 Tf\s+([-\d.]+) ([-\d.]+) ([-\d.]+) ([-\d.]+) ([-\d.]+) ([-\d.]+) Tm\s+(.*?)(?:Tj|TJ)',b,re.S)
        if not m: continue
        a,bb,c,d,e,f=[float(m.group(i)) for i in range(1,7)]
        body=m.group(7)
        if body.strip().startswith('['):
            txt=''.join(unesc(x) for x in re.findall(r'\((.*?)(?<!\\)\)',body))
        else:
            mm=re.search(r'\((.*)\)\s*$',body.strip(),re.S)
            txt=unesc(mm.group(1)) if mm else ''
        texts.append(dict(x=e,y=f,a=a,b=bb,text=txt))
    return rects,texts

def matrix(pdf):
    rects,texts=parse(pdf)
    from collections import Counter
    cnt=Counter((round(r[2],2),round(r[3],2)) for r in rects)
    (cw,ch),_=cnt.most_common(1)[0]
    cells=[r for r in rects if round(r[2],2)==cw and round(r[3],2)==ch]
    legend=sorted([r for r in rects if (round(r[2],2),round(r[3],2))!=(cw,ch)],key=lambda r:r[1])
    lx=np.median([r[0] for r in legend]); lw=np.median([r[2] for r in legend])
    lyc=np.array([r[1]+r[3]/2 for r in legend]); lcol=np.array([r[4] for r in legend])
    # tick labels: numeric text to the right of legend
    ticks=[]
    for t in texts:
        s=t['text'].replace('−','-')
        if re.fullmatch(r'-?\d+(\.\d+)?',s) and t['x']>lx+lw-2 and t['x']<lx+lw+40 and t['b']==0:
            v=float(s)
            if -1.01<=v<=1.01: ticks.append((t['y'],v))
    ticks=sorted(set(ticks))
    ty=np.array([t[0] for t in ticks]); tv=np.array([t[1] for t in ticks])
    A=np.polyfit(ty,tv,1)
    lval=np.polyval(A,lyc-3.0)  # text baseline ~3pt below center
    # refine: use center directly, offset estimated from fit residual is small
    lval=np.polyval(A,lyc)-np.polyval(A,0)*0 
    lval=A[0]*lyc+A[1]
    # baseline offset correction: tick text baseline sits ~ -2 below the y it labels
    xs=sorted(set(round(r[0],2) for r in cells)); ys=sorted(set(round(r[1],2) for r in cells))
    xi={v:i for i,v in enumerate(xs)}; yi={v:i for i,v in enumerate(ys)}
    M=np.full((len(ys),len(xs)),np.nan)
    for x,y,w,h,c in cells:
        d=np.sum((lcol-np.array(c))**2,axis=1)
        j=int(np.argmin(d))
        if d[j]>0.02:   # grey/NA
            continue
        # linear interp between neighbours by color distance
        M[yi[round(y,2)],xi[round(x,2)]]=lval[j]
    # labels
    rowlab={}; collab={}
    for t in texts:
        if t['b']==0 and t['a']==6.0 and t['x']>max(xs)+cw-1:
            rowlab[round(t['y'],2)]=t['text']
        elif abs(t['b']-t['a'])<0.01 and t['a']<6:
            collab[round(t['x'],2)]=t['text']
    rows=[]
    for y in ys:
        best=min(rowlab,key=lambda ry:abs(ry-(y+ch/2-2.2))) if rowlab else None
        rows.append(rowlab.get(best) if best is not None and abs(best-(y+ch/2-2.2))<4 else None)
    cols=[]
    ck=sorted(collab)
    for x in xs:
        best=min(ck,key=lambda cx:abs(cx-(x+cw/2))) if ck else None
        cols.append(collab.get(best) if best is not None and abs(best-(x+cw/2))<8 else None)
    title=[t['text'] for t in texts if t['a']>=12]
    return dict(M=M,rows=rows,cols=cols,title=title,ticks=ticks,nlegend=len(legend),
                lval=(float(lval.min()),float(lval.max())))
