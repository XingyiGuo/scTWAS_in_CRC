import pandas as pd, numpy as np
from scipy import stats
from base import C,E,OCC,OCC_T,O,OP

def bh(p):
    p=np.asarray(p,float); o=np.argsort(p); q=np.empty_like(p); m=len(p); prev=1.0
    for i in range(m-1,-1,-1):
        prev=min(prev,p[o[i]]*m/(i+1)); q[o[i]]=prev
    return q

TFS=sorted(C.TF.unique())
def key_tfs(ntop=20):
    out={}
    for tis in ['Normal','Adenoma']:
        deg=E[E.tissue==tis].groupby('TF').gene.nunique()
        o=OCC_T[OCC_T.tissue==tis].set_index('TF').occ_total
        m=pd.DataFrame({'deg':deg}).join(o).fillna(0)
        m=m.reindex(TFS).fillna(0)
        m['sc']=(m.deg.rank(pct=True)+m.occ_total.rank(pct=True))/2
        out[tis]=list(m.sort_values('sc',ascending=False).head(ntop).index)
    return out

MINGENE=10
CTS={'Normal':['STM','TAC','ABS','GOB','EE','TUF','CT'],
     'Adenoma':['STM','TAC','ABS','GOB','ASC','SSC','EE','TUF','CT']}
def spec(tis):
    cts=[c for c in CTS[tis] if C[(C.tissue==tis)&(C.cell_type==c)].gene.nunique()>=MINGENE]
    d=C[(C.tissue==tis)&(C.cell_type.isin(cts))]
    ng=d.groupby('cell_type').gene.nunique()[cts]
    deg=E[(E.tissue==tis)&(E.cell_type.isin(cts))].groupby(['TF','cell_type']).gene.nunique().unstack()
    deg=deg.reindex(index=TFS,columns=cts).fillna(0).astype(int)
    rows=[]
    for tf in TFS:
        for c in cts:
            a=deg.loc[tf,c]; b=ng[c]-a
            oc=deg.loc[tf,[x for x in cts if x!=c]].sum(); od=ng[[x for x in cts if x!=c]].sum()-oc
            p=stats.fisher_exact([[a,b],[oc,od]],alternative='greater')[1]
            rows.append(dict(tissue=tis,TF=tf,cell_type=c,n_linked=int(a),n_genes=int(ng[c]),
                             rate=a/ng[c],rate_other=oc/(oc+od) if (oc+od)>0 else 0,p=p))
    R=pd.DataFrame(rows); R['q']=bh(R.p.values)
    return cts,ng,deg,R

if __name__=='__main__':
    K=key_tfs(20); KEY=sorted(set(K['Normal'])|set(K['Adenoma']))
    print('key TFs normal:',K['Normal']); print('key TFs adenoma:',K['Adenoma'])
    print('union',len(KEY),'shared',len(set(K['Normal'])&set(K['Adenoma'])),sorted(set(K['Normal'])&set(K['Adenoma'])))
    print(' normal-only',sorted(set(K['Normal'])-set(K['Adenoma'])),' adenoma-only',sorted(set(K['Adenoma'])-set(K['Normal'])))
    pd.Series(KEY,name='TF').to_csv('key_tfs.csv',index=False)
    print()
    for tis in ['Normal','Adenoma']:
        cts,ng,deg,R=spec(tis); R.to_csv(f'spec_{tis}.csv',index=False)
        s=R[(R.p<0.05)&(R.n_linked>=3)]
        print('===',tis,'cell types:',cts,' nominal:',len(s),' FDR<5%:',int(((R.q<0.05)&(R.n_linked>=3)).sum()))
        for c in cts:
            x=s[s.cell_type==c].sort_values('p')
            print(f'  {c:4s} ({ng[c]}):', ', '.join(f'{r.TF} {r.n_linked}/{r.n_genes} P={r.p:.1e}{"*" if r.q<0.05 else ""}' for _,r in x.head(9).iterrows()))
        print()
    # Spearman occupancy vs co-expression degree, per cell type
    rows=[]
    for tis in ['Normal','Adenoma']:
        for c in CTS[tis]:
            d=C[(C.tissue==tis)&(C.cell_type==c)]
            if d.empty: continue
            deg=E[(E.tissue==tis)&(E.cell_type==c)].groupby('TF').gene.nunique()
            oc=OCC[(OCC.tissue==tis)&(OCC.cell_type==c)].set_index('TF').n_occupied
            j=pd.DataFrame(index=TFS).join(deg.rename('deg')).join(oc).fillna(0)
            if j.deg.std()==0: continue
            rho,p=stats.spearmanr(j.n_occupied,j.deg)
            rows.append((tis,c,round(rho,3),p))
    S=pd.DataFrame(rows,columns=['tissue','cell_type','rho','p'])
    print('Spearman(occupancy, co-expression degree) per cell type:'); print(S.to_string(index=False))
    for tis in ['Normal','Adenoma']:
        print(tis,'Fisher-combined P =',round(stats.combine_pvalues(S[S.tissue==tis].p)[1],4))
