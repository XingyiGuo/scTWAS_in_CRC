import pandas as pd, numpy as np
from scipy import stats
C=pd.read_csv('corr_long.csv')          # exact Pearson r, every TF x gene x cell type
OP=pd.read_csv('occ_pairs.csv')         # TF x gene x cell type occupancy events (from the new tables)
O=pd.read_csv('occ_gene.csv')           # per gene: lead variant, GWAS P, n TFs

# ---- per-TF occupancy summary: how many risk genes' lead variants each TF occupies ----
OCC = OP.groupby(['tissue','cell_type','TF']).gene.nunique().rename('n_occupied').reset_index()
OCC_T = OCC.groupby(['tissue','TF']).n_occupied.sum().rename('occ_total').reset_index()

# ---- co-expression layer over ALL risk genes ----
C['coexpr'] = C.r > 0.5
def pval(r,n):
    r=np.clip(np.abs(r),0,0.99999); t=r*np.sqrt(n-2)/np.sqrt(1-r**2); return 2*stats.t.sf(t,n-2)
C['p']=pval(C.r.values,C.n_donors.values)
def bh(p):
    p=np.asarray(p,float); ok=~np.isnan(p); q=np.full(len(p),np.nan)
    pp=p[ok]; o=np.argsort(pp); qq=np.empty_like(pp); m=len(pp); prev=1.0
    for i in range(m-1,-1,-1):
        prev=min(prev,pp[o[i]]*m/(i+1)); qq[o[i]]=prev
    q[ok]=qq; return q
C['q']=bh(C.p.values)
E = C[C.coexpr].copy()                  # co-expression edges

# annotate each edge with whether that TF also occupies that gene's lead variant
occset=set(zip(OP.tissue,OP.cell_type,OP.TF,OP.gene))
E['occupies_lead_variant']=[(t,c,f,g) in occset for t,c,f,g in zip(E.tissue,E.cell_type,E.TF,E.gene)]

if __name__=='__main__':
    print('co-expression edges (r>0.5):')
    print(E.groupby('tissue').agg(links=('r','size'),pairs=('gene','size'),
          TFs=('TF','nunique'),genes=('gene','nunique')).to_string())
    print(' unique TF-gene pairs:',{t:E[E.tissue==t].groupby(['TF','gene']).ngroups for t in ['Normal','Adenoma']})
    print(' all pass FDR<5%:',bool((E.q<0.05).all()),' max q=',round(E.q.max(),5))
    print()
    print('per-cell-type:')
    print(C.groupby(['tissue','cell_type']).apply(lambda d: pd.Series({
        'genes':d.gene.nunique(),'tests':len(d),'links':int(d.coexpr.sum()),
        'pct':round(100*d.coexpr.mean(),1),'TFs':d[d.coexpr].TF.nunique(),
        'genes_linked':d[d.coexpr].gene.nunique(),'neg':int((d.r<-0.5).sum())}),include_groups=False).to_string())
    print()
    print('per-TF occupancy totals (top 12):')
    for tis in ['Normal','Adenoma']:
        print(' ',tis, OCC_T[OCC_T.tissue==tis].sort_values('occ_total',ascending=False).head(12)
              .set_index('TF').occ_total.to_dict())
    print()
    print('occupancy grand totals:', OCC_T.groupby('tissue').occ_total.sum().to_dict())
