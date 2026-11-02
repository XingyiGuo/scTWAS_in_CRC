import pandas as pd, numpy as np, json
from scipy import stats
df=pd.read_csv('tf_gene_correlations.csv'); occ=pd.read_csv('tf_leadvariant_occupancy.csv')
meta=pd.read_csv('file_meta.csv')
meta['n']=meta.title.str.extract(r'n = (\d+)').astype(int)
N=meta.set_index(['tissue','cell_type']).n.to_dict()
out={}
rows=[]
for (tis,ct),d in df.groupby(['tissue','cell_type']):
    n=N[(tis,ct)]; k=len(d)
    obs=int((d.r>0.5).sum()); neg=int((d.r<-0.5).sum())
    t=0.5*np.sqrt(n-2)/np.sqrt(1-0.25); p1=stats.t.sf(t,n-2)
    exp=k*p1
    # binomial test
    pb=stats.binomtest(obs,k,p1,alternative='greater').pvalue
    rows.append(dict(tissue=tis,cell_type=ct,n_samples=n,n_genes=d.gene.nunique(),n_tests=k,
        obs_pos=obs,obs_neg=neg,exp_null=round(exp,2),fold=round(obs/exp,1) if exp>0 else np.nan,
        p_binom=pb, pct=round(100*obs/k,1)))
E=pd.DataFrame(rows)
E.to_csv('enrichment_stats.csv',index=False)
print(E.to_string())
print()
# FDR per pair using per-cell-type n
def pr(r,n):
    r=np.clip(np.abs(r),0,0.9999); t=r*np.sqrt(n-2)/np.sqrt(1-r**2); return 2*stats.t.sf(t,n-2)
df['n']=[N[(t,c)] for t,c in zip(df.tissue,df.cell_type)]
df['p']=pr(df.r.values,df.n.values)
def bh(p):
    p=np.asarray(p); o=np.argsort(p); q=np.empty_like(p)
    m=len(p); prev=1
    for i in range(m-1,-1,-1):
        prev=min(prev,p[o[i]]*m/(i+1)); q[o[i]]=prev
    return q
df['q']=bh(df.p.values)
df['sig']=(df.q<0.05)&(df.r>0.5)
print('pairs r>0.5 & FDR<5%:', df.sig.sum(), 'of', (df.r>0.5).sum())
print(df[df.r>0.5].groupby('tissue').apply(lambda d: pd.Series({'r>0.5':len(d),'FDR<5%':int(d.sig.sum())}),include_groups=False).to_string())
df.to_csv('tf_gene_correlations_annot.csv',index=False)
