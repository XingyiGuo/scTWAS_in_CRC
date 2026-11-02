import numpy as np, bars2, json
from PIL import Image
Image.MAX_IMAGE_PIXELS=None
NL='/mnt/user-data/uploads/TFgenes/NL_ALL7CT_TF_leadVariantOccupation_vs_corrGT0.5_orderedByLeadVariant_large.png'
AD='/mnt/user-data/uploads/TFgenes/ADSER_ALL9CT_TF_leadVariantOccupation_vs_corrGT0.5_orderedByLeadVariant_large.png'
NLCT=['ABS','GOB','CT','EE','TUF','TAC','STM']
ADCT=['ABS','GOB','CT','EE','TUF','TAC','STM','ASC','SSC']

def panel_slots(a,rgbsum,r,c):
    y0,y1,x0,x1=r*2400,(r+1)*2400,c*2400,(c+1)*2400
    d=bars2.panel_bars(a,y0,y1,x0,x1)
    if d is None: return None
    base=d['baseline']
    # axis line = first full-width dark row below base
    axis=None
    for y in range(base+60,base+160):
        if (rgbsum[y,x0:x1]<600).sum()>1500: axis=y; break
    ty=axis+6
    row=rgbsum[ty,x0:x1]<600
    runs=bars2.runs(row)
    ticks=[(s+e)/2+x0 for s,e in runs]
    hs=[b['h'] for b in d['bars']]
    sc,err=bars2.fit_scale(hs)
    w=np.median([b['xe']-b['xs'] for b in d['bars']])
    red=[0]*len(ticks); blue=[0]*len(ticks)
    for b in d['bars']:
        cen=b['xc']+(w/2 if b['kind']=='red' else -w/2)
        i=int(np.argmin([abs(cen-t) for t in ticks]))
        v=int(round(b['h']/sc))
        if b['kind']=='red': red[i]=v
        else: blue[i]=v
    return dict(ticks=ticks,red=red,blue=blue,axis=axis,base=base,scale=sc,err=err,box=(y0,y1,x0,x1))

def all_panels(path,cts):
    im=Image.open(path).convert('RGB'); a=np.array(im); s=a.sum(axis=2)
    out={};k=0
    for r in range(3):
        for c in range(3):
            if k>=len(cts): break
            d=panel_slots(a,s,r,c)
            if d: out[cts[k]]=d
            k+=1
    return out
if __name__=='__main__':
    res={}
    for path,cts,tis in ((NL,NLCT,'Normal'),(AD,ADCT,'Adenoma')):
        for ct,d in all_panels(path,cts).items():
            res[f'{tis}|{ct}']={k:v for k,v in d.items()}
            print(tis,ct,'nticks',len(d['ticks']),'sum red',sum(d['red']),'sum blue',sum(d['blue']),'axis',d['axis'])
    json.dump(res,open("slot_values.json","w"),default=lambda o:int(o) if hasattr(o,"item") else o)
