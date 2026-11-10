"""New contrast sets requested on top of res/new_DE:

  C1-C7  normal colon cell type  vs  MSI-H carcinoma   and  vs MSS carcinoma
  E1-E7  normal colon cell type  vs  ASC in adenoma    and  vs SSC in serrated polyp
  L1-L2  ASC in adenoma / SSC in serrated  vs  MSI-H and MSS   (relabelled from L3-L6)

C1..C7 / E1..E7 are indexed by the normal cell type in the order ABS, GOB, CT, EE, STM, TUF, TAC.
Fitting reuses the early/late contrast functions, so these rows are directly
comparable with the existing E1-E3 / L1-L6 results.
"""
import importlib.util
import os
from pathlib import Path
import numpy as np
import pandas as pd

NEW_DE = r'L:\scRNA-seq\VUMC_HTAN\res\new_DE'
CT7 = ['ABS', 'GOB', 'CT', 'EE', 'STM', 'TUF', 'TAC']
EARLY_COH = ['DIS', 'VAL', 'PCA83']
MIN_CELLS = 10

_RUN_DE = Path(__file__).with_name('02_fit_early_late_contrasts.py')
spec = importlib.util.spec_from_file_location('run_de', _RUN_DE)
run_de = importlib.util.module_from_spec(spec)
spec.loader.exec_module(run_de)


def load_pb(path='pseudobulk_raw.npz'):
    z = np.load(path, allow_pickle=True)
    cohs = sorted({k.rsplit('_', 1)[0] for k in z.files})
    return {c: dict(units=z[f'{c}_units'], genes=z[f'{c}_genes'], mat=z[f'{c}_mat']) for c in cohs}


def arms(U):
    """Unit sets for every arm used below, all at the >= 10 cell threshold."""
    u = U[U.n_cells >= MIN_CELLS]
    out = {}
    for ct in CT7:
        out[f'Normal_{ct}'] = u[(u.cohort.isin(EARLY_COH)) & (u.group == 'Normal')
                                & (u.cell_type == ct)]
    out['MSI-H'] = u[(u.cohort == 'ABN') & (u.cell_type == 'CARCINOMA') & (u.polyp_class == 'MSI-H')]
    out['MSS'] = u[(u.cohort == 'ABN') & (u.cell_type == 'CARCINOMA') & (u.polyp_class == 'MSS')]
    # ASC / SSC drawn from all four files, restricted to the matching lesion type
    out['ASC'] = u[(u.cell_type == 'ASC') & (u.group == 'Adenoma')]
    out['SSC'] = u[(u.cell_type == 'SSC') & (u.group == 'Serrated')]
    return out


def run(PB, U, log=print):
    A = arms(U)
    res_all = []
    jobs = []
    for i, ct in enumerate(CT7, start=1):
        for test, suff in (('MSI-H', 'MSIH'), ('MSS', 'MSS')):
            jobs.append((f'C{i}_{suff}', test, f'Normal_{ct}', ct, 'normal_vs_carcinoma'))
    for i, ct in enumerate(CT7, start=1):
        for test in ('ASC', 'SSC'):
            jobs.append((f'E{i}_{test}', test, f'Normal_{ct}', ct, 'normal_vs_polyp_celltype'))

    for tag, test, ref, ct, stage in jobs:
        ut, ur = A[test], A[ref]
        if len(ut) < 2 or len(ur) < 3:
            log(f'{tag}: SKIPPED ({len(ut)} vs {len(ur)} specimens)')
            continue
        u = pd.concat([ut, ur])
        meta = pd.DataFrame({'group': [test] * len(ut) + [ref] * len(ur),
                             'cohort': u.cohort.values}, index=u.unit.values)
        counts = run_de.build_matrix(PB, list(u.unit.values))
        res, design, ng = run_de.fit(counts, meta, test, ref, design_extra=('cohort',))
        out = run_de.tidy(res, cell_type=ct, contrast=f'{tag}: {test} vs {ref}', stage=stage,
                          design=design, n_genes=ng, meta=meta, test=test, ref=ref, units=U)
        out.insert(3, 'contrast_tag', tag)
        res_all.append(out)
        log(f'{tag} {design} n={len(ut)}v{len(ur)} genes={ng} '
            f'padj<0.05={int((res.padj < 0.05).sum())}')
    return pd.concat(res_all, ignore_index=True)


RELABEL = {'L3': ('L1_MSIH', 'MSI-H', 'ASC'), 'L4': ('L1_MSS', 'MSS', 'ASC'),
           'L5': ('L2_MSIH', 'MSI-H', 'SSC'), 'L6': ('L2_MSS', 'MSS', 'SSC')}


def relabel_late(path=None):
    """Pull the already-fitted MSI-H/MSS vs ASC/SSC contrasts and give them the L1/L2 labels."""
    path = path or os.path.join(NEW_DE, 'de_all_contrasts.parquet')
    d = pd.read_parquet(path)
    d = d[d.contrast.str.split(':').str[0].isin(RELABEL)].copy()
    old = d.contrast.str.split(':').str[0]
    d.insert(3, 'contrast_tag', old.map(lambda t: RELABEL[t][0]))
    d['contrast'] = [f'{RELABEL[t][0]}: {RELABEL[t][1]} vs {RELABEL[t][2]}' for t in old]
    d['stage'] = 'polyp_celltype_vs_carcinoma'
    return d


def summarise(d):
    g = d.groupby(['contrast_tag', 'contrast', 'stage', 'cell_type', 'group1', 'group2',
                   'design', 'n_specimens_grp1', 'n_specimens_grp2', 'n_genes_tested'],
                  dropna=False)
    s = g.apply(lambda x: pd.Series({
        'n_padj_lt_0.05': int((x.padj < 0.05).sum()),
        'n_padj_lt_0.05_absLFC_gt_0.5': int(((x.padj < 0.05) & (x.log2FC.abs() > 0.5)).sum()),
        'n_up': int(((x.padj < 0.05) & (x.log2FC > 0.5)).sum()),
        'n_down': int(((x.padj < 0.05) & (x.log2FC < -0.5)).sum()),
        'n_padj_NA': int(x.padj.isna().sum())}), include_groups=False).reset_index()
    s['underpowered_flag'] = np.where(s[['n_specimens_grp1', 'n_specimens_grp2']].min(axis=1) < 3,
                                      'yes: fewer than 3 specimens in one arm', '')
    return s.sort_values('contrast_tag')
