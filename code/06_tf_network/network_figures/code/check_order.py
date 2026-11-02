import json,sys,pandas as pd,numpy as np
S=json.load(open('slot_values.json'))
df=pd.read_csv('tf_gene_correlations.csv')
def check(tis,ct,order):
    d=S[f'{tis}|{ct}']
    assert len(order)==51, len(order)
    red=d['red']; blue=d['blue']
    c=df[(df.tissue==tis)&(df.cell_type==ct)].groupby('TF').r.apply(lambda s:int((s>0.5).sum())).to_dict()
    bad=[(i,tf,red[i],c.get(tf,0)) for i,tf in enumerate(order) if red[i]!=c.get(tf,0)]
    mono=all(blue[i]>=blue[i+1] for i in range(50))
    tfs_pdf=set(df[(df.tissue==tis)&(df.cell_type==ct)].TF)
    print(f'{tis} {ct}: red-mismatch={len(bad)} blue_monotone={mono} unknown_names={sorted(set(order)-tfs_pdf)} dup={len(order)!=len(set(order))}')
    if bad: print('   ',bad[:12])
    return bad
