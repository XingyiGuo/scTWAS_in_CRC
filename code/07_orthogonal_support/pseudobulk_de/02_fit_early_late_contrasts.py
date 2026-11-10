"""Specimen-level pseudobulk DE: early stage E1-E3, late stage L1-L6."""
import numpy as np, pandas as pd, warnings
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

FILENAME = {'DIS': 'VUMC_HTAN_DIS_EPI_V2.h5ad', 'VAL': 'VUMC_HTAN_VAL_EPI_V2.h5ad',
            'PCA83': 'PCA83samples.h5ad', 'ABN': 'VUMC_ABNORMALS_EPI_V2.h5ad'}
CT7 = ['ABS', 'GOB', 'CT', 'EE', 'STM', 'TUF', 'TAC']
EARLY_COH = ['DIS', 'VAL', 'PCA83']
MIN_CELLS = 10
MIN_FRAC_DETECTED = 0.10


def build_matrix(PB, unit_ids):
    """Counts for the requested units, on the gene intersection of contributing files."""
    cohs = sorted({u.split('|')[0] for u in unit_ids})
    genes = set(PB[cohs[0]]['genes'])
    for c in cohs[1:]:
        genes &= set(PB[c]['genes'])
    genes = np.array(sorted(genes))
    blocks = []
    for c in cohs:
        gi = pd.Index(PB[c]['genes']).get_indexer(genes)
        ui = pd.Index(PB[c]['units'])
        keep = [u for u in unit_ids if u.split('|')[0] == c]
        ri = ui.get_indexer(keep)
        if (ri < 0).any():
            missing = [unit for unit, row in zip(keep, ri) if row < 0]
            raise ValueError(f'Pseudobulk counts are missing units: {missing}')
        blocks.append(pd.DataFrame(PB[c]['mat'][np.ix_(ri, gi)], index=keep, columns=genes))
    return pd.concat(blocks).loc[list(unit_ids)]


def fit(counts, meta, test, ref, design_extra=('cohort',)):
    """DESeq2 Wald test of group=test vs group=ref, with apeglm-shrunk fold changes."""
    keep = (counts > 0).mean(axis=0) >= MIN_FRAC_DETECTED
    counts = counts.loc[:, keep]
    terms = ['group'] + [f for f in design_extra if meta[f].nunique() > 1]
    # a covariate level present in only one group arm makes the design rank-deficient
    for f in list(terms[1:]):
        ct = pd.crosstab(meta[f], meta['group'])
        if (ct.gt(0).sum(axis=1) == 1).all():
            terms.remove(f)
    design = '~' + '+'.join(terms)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        dds = DeseqDataSet(counts=counts.astype(int), metadata=meta, design=design,
                           refit_cooks=True, quiet=True, n_cpus=8)
        dds.deseq2()
        st = DeseqStats(dds, contrast=['group', test, ref], quiet=True, n_cpus=8)
        st.summary()
        res = st.results_df.copy()
        res['log2FC_MLE'] = res['log2FoldChange']
        # The contrast can be the negative of the fitted group coefficient.
        indices = np.flatnonzero(st.contrast_vector)
        if len(indices) != 1 or abs(st.contrast_vector[indices[0]]) != 1:
            raise ValueError(f'Expected a single group coefficient for {test} vs {ref}')
        idx = indices[0]
        st.lfc_shrink(coeff=dds.varm['LFC'].columns[idx])
        res['log2FoldChange'] = st.contrast_vector[idx] * st.results_df['log2FoldChange']
        res['lfcSE'] = st.results_df['lfcSE']
    return res, design, counts.shape[1]


def tidy(res, *, cell_type, contrast, stage, design, n_genes, meta, test, ref, units):
    nc = units.set_index('unit')['n_cells']
    g1 = meta.index[meta.group == test]
    g2 = meta.index[meta.group == ref]
    src = lambda ix: '; '.join(sorted({FILENAME[u.split('|')[0]] for u in ix}))
    out = res.rename_axis('gene').reset_index().rename(
        columns={'log2FoldChange': 'log2FC', 'lfcSE': 'SE'}
    )
    out.insert(1, 'cell_type', cell_type)
    out.insert(2, 'contrast', contrast)
    out.insert(3, 'stage', stage)
    out['n_specimens_grp1'] = len(g1)
    out['n_specimens_grp2'] = len(g2)
    out['mean_cells_grp1'] = round(float(nc.loc[g1].mean()), 1)
    out['mean_cells_grp2'] = round(float(nc.loc[g2].mean()), 1)
    out['source_grp1'] = src(g1)
    out['source_grp2'] = src(g2)
    out['group1'] = test
    out['group2'] = ref
    out['design'] = design
    out['n_genes_tested'] = n_genes
    return out[['gene', 'cell_type', 'contrast', 'stage', 'group1', 'group2',
                'log2FC', 'SE', 'stat', 'pvalue', 'padj', 'log2FC_MLE', 'baseMean',
                'n_specimens_grp1', 'n_specimens_grp2', 'mean_cells_grp1', 'mean_cells_grp2',
                'source_grp1', 'source_grp2', 'design', 'n_genes_tested']]


def run_all(PB, units, log=print):
    res_all, skipped = [], []
    # ---------- early stage ----------
    early = units[(units.cohort.isin(EARLY_COH)) & (units.group != 'Unclassified')]
    specs = [('E1', 'Polyp', ['Adenoma', 'Serrated'], 'Normal', MIN_CELLS),
             ('E2', 'Adenoma', ['Adenoma'], 'Normal', MIN_CELLS),
             ('E3', 'Serrated', ['Serrated'], 'Normal', MIN_CELLS),
             ('E1_EE5', 'Polyp', ['Adenoma', 'Serrated'], 'Normal', 5)]
    for tag, test, test_groups, ref, mincells in specs:
        cts = ['EE'] if tag == 'E1_EE5' else CT7
        for ct in cts:
            u = early[(early.cell_type == ct) & (early.n_cells >= mincells) &
                      (early.group.isin(test_groups + [ref]))]
            if u.group.nunique() < 2 or (u.group == ref).sum() < 3 or (u.group != ref).sum() < 3:
                skipped.append((ct, tag, 'too few specimens'))
                continue
            meta = pd.DataFrame({'group': np.where(u.group == ref, ref, test),
                                 'cohort': u.cohort.values,
                                 'participant': u.participant.values}, index=u.unit.values)
            counts = build_matrix(PB, list(u.unit.values))
            res, design, ng = fit(counts, meta, test, ref)
            res_all.append(tidy(res, cell_type=ct, contrast=f'{tag}: {test} vs {ref}',
                                stage='early', design=design, n_genes=ng, meta=meta,
                                test=test, ref=ref, units=units))
            log(f'{tag} {ct}: {design}  n={sum(meta.group==test)}v{sum(meta.group==ref)}  '
                f'genes={ng}  padj<0.05={int((res.padj<0.05).sum())}')
    # ---------- late stage ----------
    abn = units[units.cohort == 'ABN']
    arms = {
        'ASC_adenoma': abn[(abn.cell_type == 'ASC') & (abn.polyp_class.isin(['TA', 'TVA']))],
        'SSC_serrated': abn[(abn.cell_type == 'SSC') & (abn.polyp_class.isin(['SSL', 'TSA', 'HP']))],
        'Carcinoma': abn[abn.cell_type == 'CARCINOMA'],
        'MSI-H': abn[(abn.cell_type == 'CARCINOMA') & (abn.polyp_class == 'MSI-H')],
        'MSS': abn[(abn.cell_type == 'CARCINOMA') & (abn.polyp_class == 'MSS')]}
    arms = {k: v[v.n_cells >= MIN_CELLS] for k, v in arms.items()}
    late = [('L1', 'Carcinoma', 'ASC_adenoma', 'ASC'), ('L2', 'Carcinoma', 'SSC_serrated', 'SSC'),
            ('L3', 'MSI-H', 'ASC_adenoma', 'ASC'), ('L4', 'MSS', 'ASC_adenoma', 'ASC'),
            ('L5', 'MSI-H', 'SSC_serrated', 'SSC'), ('L6', 'MSS', 'SSC_serrated', 'SSC')]
    for tag, test, ref, ct_label in late:
        u = pd.concat([arms[test], arms[ref]])
        meta = pd.DataFrame({'group': [test] * len(arms[test]) + [ref] * len(arms[ref]),
                             'cohort': u.cohort.values}, index=u.unit.values)
        counts = build_matrix(PB, list(u.unit.values))
        res, design, ng = fit(counts, meta, test, ref, design_extra=())
        res_all.append(tidy(res, cell_type=ct_label, contrast=f'{tag}: {test} vs {ref}',
                            stage='late', design=design, n_genes=ng, meta=meta,
                            test=test, ref=ref, units=units))
        log(f'{tag} {ct_label}: {design}  n={len(arms[test])}v{len(arms[ref])}  '
            f'genes={ng}  padj<0.05={int((res.padj<0.05).sum())}')
    return pd.concat(res_all, ignore_index=True), skipped
