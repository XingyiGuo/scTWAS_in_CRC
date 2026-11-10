"""Rebuild specimen x cell-type pseudobulk count matrices from the four VUMC/HTAN h5ad files.

Reproduces the aggregation behind res/new_DE/pseudobulk_units.csv:
  * PCA83 cluster "20" dropped
  * ABN cell types transferred from DIS/VAL on (HTAN Specimen ID, barcode before first dash);
    VAL barcodes live in obs['Cells'], not the obs index
  * MSI-H / MSS specimens in ABN labelled CARCINOMA
Asserts the rebuilt units match the shipped table before the expensive read.
Writes pseudobulk_raw.npz (per-cohort mat / units / genes).
"""
import os
import time
import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp

ROOT = r'L:\scRNA-seq\VUMC_HTAN\data'
UNITS_CSV = r'L:\scRNA-seq\VUMC_HTAN\res\new_DE\pseudobulk_units.csv'
FILES = {'DIS': 'VUMC_HTAN_DIS_EPI_V2.h5ad', 'VAL': 'VUMC_HTAN_VAL_EPI_V2.h5ad',
         'PCA83': 'PCA83samples.h5ad', 'ABN': 'VUMC_ABNORMALS_EPI_V2.h5ad'}
CHUNK = 1000


def read_obs(path):
    with h5py.File(path, 'r') as h:
        g = h['obs']
        cats = g['__categories'] if '__categories' in g else None
        d = {}
        for k in g.keys():
            if k == '__categories':
                continue
            v = g[k][:]
            if cats is not None and k in cats:
                c = np.array([x.decode() if isinstance(x, bytes) else x for x in cats[k][:]])
                v = pd.Categorical.from_codes(v, c)
            elif v.dtype.kind == 'S':
                v = np.array([x.decode() for x in v])
            d[k] = v
        idx = d.pop('_index', None)
        df = pd.DataFrame(d)
        if idx is not None:
            df.index = [x.decode() if isinstance(x, bytes) else x for x in idx]
        return df


def read_genes(path):
    with h5py.File(path, 'r') as h:
        k = 'Genes' if 'Genes' in h['var'] else '_index'
        return np.array([x.decode() if isinstance(x, bytes) else x for x in h['var'][k][:]])


def bare(arr):
    return np.array([str(x).split('-')[0] for x in arr])


def build_cell_map(OB):
    parts = []
    for coh, scol, ccol in [('DIS', 'HTAN Specimen ID', 'Cell_Type'),
                            ('VAL', 'HTAN Specimen ID', 'Cell_Type'),
                            ('PCA83', 'SPECID.htan', 'Cell Type')]:
        o = OB[coh]
        parts.append(pd.DataFrame({'cohort': coh, 'specimen': o[scol].astype(str).values,
                                   'cell_type': o[ccol].astype(str).values,
                                   'row': np.arange(len(o))}))
    bc_of = {'DIS': OB['DIS'].index.values, 'VAL': OB['VAL']['Cells'].astype(str).values}
    ref = pd.concat([pd.DataFrame({'specimen': OB[c]['HTAN Specimen ID'].astype(str).values,
                                   'bc': bare(bc_of[c]),
                                   'ct': OB[c]['Cell_Type'].astype(str).values})
                     for c in ['DIS', 'VAL']])
    ref = ref.drop_duplicates(['specimen', 'bc']).set_index(['specimen', 'bc'])['ct']
    oa = OB['ABN']
    key = pd.MultiIndex.from_arrays([oa['HTAN Specimen ID'].astype(str).values,
                                     bare(oa.index.values)])
    ct = pd.Series(ref.reindex(key).values, index=np.arange(len(oa)), dtype=object)
    is_carc = pd.Series(oa['Sample_Classification'].astype(str).values).isin(['MSI-H', 'MSS']).values
    ct[is_carc & ct.isna()] = 'CARCINOMA'
    parts.append(pd.DataFrame({'cohort': 'ABN', 'specimen': oa['HTAN Specimen ID'].astype(str).values,
                               'cell_type': ct.values, 'row': np.arange(len(oa))}))
    cell = pd.concat(parts, ignore_index=True)
    cell = cell[cell.cell_type.notna() & (cell.cell_type != 'nan') & (cell.cell_type != '20')]
    cell['unit'] = cell.cohort + '|' + cell.specimen + '|' + cell.cell_type
    return cell


def aggregate(cell, GN, log=print):
    PB = {}
    t0 = time.time()
    for coh, fn in FILES.items():
        sub = cell[cell.cohort == coh]
        ulist = sorted(sub.unit.unique())
        umap = {u: i for i, u in enumerate(ulist)}
        rows = sub['row'].values
        cols = sub['unit'].map(umap).values
        acc = np.zeros((len(ulist), len(GN[coh])), dtype=np.int64)
        with h5py.File(os.path.join(ROOT, fn), 'r') as h:
            X = h['X']
            n = X.shape[0]
            for a in range(0, n, CHUNK):
                b = min(a + CHUNK, n)
                m = (rows >= a) & (rows < b)
                if not m.any():
                    continue
                blk = X[a:b, :]
                ind = sp.csr_matrix((np.ones(m.sum(), dtype=np.float32),
                                     (cols[m], rows[m] - a)), shape=(len(ulist), b - a))
                acc += np.rint(ind.dot(blk)).astype(np.int64)
                del blk, ind
        PB[coh] = dict(units=np.array(ulist), genes=GN[coh], mat=acc)
        log(f'{coh} {acc.shape} total={acc.sum():.4g} {time.time() - t0:.0f}s')
    return PB


if __name__ == '__main__':
    OB = {k: read_obs(os.path.join(ROOT, v)) for k, v in FILES.items()}
    GN = {k: read_genes(os.path.join(ROOT, v)) for k, v in FILES.items()}
    cell = build_cell_map(OB)
    U = pd.read_csv(UNITS_CSV)
    mine = cell.groupby('unit').size()
    ref = U.set_index('unit')['n_cells']
    assert len(mine) == len(ref) == 1844, (len(mine), len(ref))
    assert mine.reindex(ref.index).equals(ref), 'unit cell counts do not match the shipped table'
    del OB
    PB = aggregate(cell, GN)
    np.savez_compressed('pseudobulk_raw.npz',
                        **{f'{k}_mat': v['mat'] for k, v in PB.items()},
                        **{f'{k}_units': v['units'] for k, v in PB.items()},
                        **{f'{k}_genes': v['genes'] for k, v in PB.items()})
    tot = {k: int(v['mat'].sum()) for k, v in PB.items()}
    print('saved pseudobulk_raw.npz', tot)
