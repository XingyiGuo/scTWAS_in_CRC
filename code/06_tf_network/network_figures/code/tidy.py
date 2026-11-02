import numpy as np, pandas as pd, glob, os, re
from build_mats import matrix
rows=[]
meta=[]
for p in sorted(glob.glob('/mnt/user-data/uploads/TFgenes/*.pdf')):
    b=os.path.basename(p); d=matrix(p)
    ct=b.split('_')[0]
    tissue='Adenoma' if ('ADSER' in b or '_AD_' in b or '_SER_' in b) else 'Normal'
    meta.append(dict(file=b,cell_type=ct,tissue=tissue,n_tf=d['nrow'],n_gene=d['ncol'],title=' | '.join(d['title'])))
    M=d['M']
    for i,tf in enumerate(d['rows']):
        for j,g in enumerate(d['cols']):
            rows.append((tissue,ct,tf,g,M[i,j]))
df=pd.DataFrame(rows,columns=['tissue','cell_type','TF','gene','r'])
df.to_csv('tf_gene_correlations.csv',index=False)
pd.DataFrame(meta).to_csv('file_meta.csv',index=False)
print(df.shape, df.r.isna().sum())
print(pd.DataFrame(meta)[['cell_type','tissue','n_tf','n_gene']].to_string())
print(df.groupby(['tissue','cell_type']).agg(genes=('gene','nunique'),tfs=('TF','nunique'),
     n_r_gt05=('r',lambda s:(s>0.5).sum()), n_absr_gt05=('r',lambda s:(s.abs()>0.5).sum()),
     n_r_lt_m05=('r',lambda s:(s<-0.5).sum())).to_string())
