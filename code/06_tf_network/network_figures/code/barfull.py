import numpy as np, pandas as pd, bars2, json
from PIL import Image
Image.MAX_IMAGE_PIXELS=None
NL='/mnt/user-data/uploads/TFgenes/NL_ALL7CT_TF_leadVariantOccupation_vs_corrGT0.5_orderedByLeadVariant_large.png'
AD='/mnt/user-data/uploads/TFgenes/ADSER_ALL9CT_TF_leadVariantOccupation_vs_corrGT0.5_orderedByLeadVariant_large.png'
NLCT=['ABS','GOB','CT','EE','TUF','TAC','STM']
ADCT=['ABS','GOB','CT','EE','TUF','TAC','STM','ASC','SSC']
def panels(path,cts):
    a=np.array(Image.open(path).convert('RGB'))
    out={}
    k=0
    for r in range(3):
        for c in range(3):
            if k>=len(cts): break
            d=bars2.panel_bars(a,r*2400,(r+1)*2400,c*2400,(c+1)*2400)
            if d is None: k+=1; continue
            hs=[b['h'] for b in d['bars']]
            s,err=bars2.fit_scale(hs)
            for b in d['bars']: b['val']=int(round(b['h']/s))
            out[cts[k]]=dict(bars=d['bars'],scale=s,err=err,box=(r*2400,(r+1)*2400,c*2400,(c+1)*2400))
            k+=1
    return a,out
if __name__=='__main__':
    res={}
    for path,cts,tis in ((NL,NLCT,'Normal'),(AD,ADCT,'Adenoma')):
        a,out=panels(path,cts)
        for ct,d in out.items():
            red=sorted([b['val'] for b in d['bars'] if b['kind']=='red'],reverse=True)
            blue=sorted([b['val'] for b in d['bars'] if b['kind']=='blue'],reverse=True)
            res[(tis,ct)]=dict(red=red,blue=blue,err=d['err'])
    df=pd.read_csv('tf_gene_correlations.csv')
    for (tis,ct),d in res.items():
        c=df[(df.tissue==tis)&(df.cell_type==ct)].groupby('TF').r.apply(lambda s:(s>0.5).sum())
        c=sorted(c[c>0].values.tolist(),reverse=True)
        match = (c==d['red'])
        print(f"{tis:8s} {ct:4s} redMATCH={match}  sum_red_png={sum(d['red'])} sum_red_pdf={sum(c)}  n_blue={len(d['blue'])} sum_blue={sum(d['blue'])} fit_err={d['err']:.4f}")
        if not match: print('   png',d['red']); print('   pdf',c)
    json.dump({f'{k[0]}|{k[1]}':v for k,v in res.items()},open('bar_values.json','w'),indent=1)
