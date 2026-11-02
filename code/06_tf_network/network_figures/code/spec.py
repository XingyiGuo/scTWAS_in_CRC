import pandas as pd, numpy as np
from scipy import stats
df=pd.read_csv('tf_gene_correlations_annot.csv'); occ=pd.read_csv('tf_leadvariant_occupancy.csv')
e=df[df.r>0.5]
CTS={'Normal':['STM','TAC','ABS','GOB','EE','TUF','CT'],
     'Adenoma':['STM','TAC','ABS','GOB','ASC','SSC','EE','TUF','CT']}
MINGENE=10
def tau(x):
    x=np.asarray(x,float)
    return np.nan if x.max()<=0 else np.sum(1-x/x.max())/(len(x)-1)
def bh(p):
    p=np.asarray(p,float); o=np.argsort(p); q=np.empty_like(p); m=len(p); prev=1.0
    for i in range(m-1,-1,-1):
        prev=min(prev,p[o[i]]*m/(i+1)); q[o[i]]=prev
    return q
def build(tis):
    cts=[c for c in CTS[tis] if df[(df.tissue==tis)&(df.cell_type==c)].gene.nunique()>=MINGENE]
    d=df[(df.tissue==tis)&(df.cell_type.isin(cts))]
    ng=d.groupby('cell_type').gene.nunique()[cts]
    deg=e[(e.tissue==tis)&(e.cell_type.isin(cts))].groupby(['TF','cell_type']).gene.nunique().unstack()
    deg=deg.reindex(index=sorted(set(df.TF)),columns=cts).fillna(0).astype(int)
    rate=deg.div(ng,axis=1)
    rows=[]
    for tf in deg.index:
        for c in cts:
            a=deg.loc[tf,c]; b=ng[c]-a
            cc=deg.loc[tf,[x for x in cts if x!=c]].sum(); dd=ng[[x for x in cts if x!=c]].sum()-cc
            p=stats.fisher_exact([[a,b],[cc,dd]],alternative='greater')[1]
            rows.append(dict(tissue=tis,TF=tf,cell_type=c,n_linked=int(a),n_genes=int(ng[c]),
                             rate=a/ng[c], rate_other=cc/(cc+dd) if (cc+dd)>0 else 0, p=p))
    R=pd.DataFrame(rows); R['q']=bh(R.p.values)
    T=rate.apply(tau,axis=1); dom=rate.idxmax(axis=1)
    S=pd.DataFrame({'tau':T,'dominant':dom,'max_rate':rate.max(axis=1),'total_links':deg.sum(axis=1)})
    return cts,ng,deg,rate,R,S
if __name__=='__main__':
    for tis in ['Normal','Adenoma']:
        cts,ng,deg,rate,R,S=build(tis)
        sig=R[(R.q<0.05)&(R.n_linked>=3)]
        print('===',tis,'cell types used:',cts)
        print('significant TF x cell-type enrichments (q<0.05, >=3 genes):',len(sig))
        for c in cts:
            s=sig[sig.cell_type==c].sort_values('q')
            print(f'  {c:4s} ({ng[c]} genes):', ', '.join(f'{r.TF} {r.n_linked}/{r.n_genes} q={r.q:.1e}' for _,r in s.head(8).iterrows()))
        print()
