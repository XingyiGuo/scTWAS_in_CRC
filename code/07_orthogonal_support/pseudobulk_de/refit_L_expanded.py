"""Refit L1-L2 with the ASC / SSC arms expanded to all four files (ASC n=75, SSC n=46).

The current L1-L2 fits draw both arms from the abnormals object (ASC n=23, SSC n=14). Carcinoma
exists only there, so the group effect stays identified within that file either way; the extra
specimens enter through the dispersion trend and the cohort coefficients. `cohort` is retained
because ABN contributes to both arms, which keeps the design full rank.
"""
import os
import pandas as pd
import run_new_contrasts as R

MIN_CELLS = R.MIN_CELLS
JOBS = [('L1_MSIH_all', 'MSI-H', 'ASC', 'ASC'), ('L1_MSS_all', 'MSS', 'ASC', 'ASC'),
        ('L2_MSIH_all', 'MSI-H', 'SSC', 'SSC'), ('L2_MSS_all', 'MSS', 'SSC', 'SSC')]


def expanded_arms(U):
    u = U[U.n_cells >= MIN_CELLS]
    return {'MSI-H': u[(u.cohort == 'ABN') & (u.cell_type == 'CARCINOMA') & (u.polyp_class == 'MSI-H')],
            'MSS': u[(u.cohort == 'ABN') & (u.cell_type == 'CARCINOMA') & (u.polyp_class == 'MSS')],
            'ASC': u[(u.cell_type == 'ASC') & (u.group == 'Adenoma')],
            'SSC': u[(u.cell_type == 'SSC') & (u.group == 'Serrated')]}


def run(PB, U, log=print):
    A = expanded_arms(U)
    out_all = []
    for tag, test, ref, ct in JOBS:
        ut, ur = A[test], A[ref]
        u = pd.concat([ut, ur])
        meta = pd.DataFrame({'group': [test] * len(ut) + [ref] * len(ur),
                             'cohort': u.cohort.values}, index=u.unit.values)
        counts = R.run_de.build_matrix(PB, list(u.unit.values))
        res, design, ng = R.run_de.fit(counts, meta, test, ref, design_extra=('cohort',))
        out = R.run_de.tidy(res, cell_type=ct, contrast=f'{tag}: {test} vs {ref} (all files)',
                            stage='polyp_celltype_vs_carcinoma', design=design, n_genes=ng,
                            meta=meta, test=test, ref=ref, units=U)
        out.insert(3, 'contrast_tag', tag)
        out['arm_scope'] = 'all_four_files'
        out_all.append(out)
        log(f'{tag} {design} n={len(ut)}v{len(ur)} genes={ng} '
            f'padj<0.05={int((res.padj < 0.05).sum())}')
    return pd.concat(out_all, ignore_index=True)


def compare(new, old):
    """Old (abnormals-only) vs new (all files) for the same four comparisons."""
    keymap = {'L1_MSIH_all': 'L1_MSIH', 'L1_MSS_all': 'L1_MSS',
              'L2_MSIH_all': 'L2_MSIH', 'L2_MSS_all': 'L2_MSS'}
    rows = []
    for ntag, otag in keymap.items():
        n = new[new.contrast_tag == ntag].set_index('gene')
        o = old[old.contrast_tag == otag].set_index('gene')
        shared = n.index.intersection(o.index)
        ns, os_ = n.loc[shared], o.loc[shared]
        sig_n = (ns.padj < 0.05) & (ns.log2FC.abs() > 0.5)
        sig_o = (os_.padj < 0.05) & (os_.log2FC.abs() > 0.5)
        either = sig_n | sig_o
        rows.append(dict(
            comparison=otag,
            n_ref_specimens_old=int(o.n_specimens_grp2.iloc[0]),
            n_ref_specimens_new=int(n.n_specimens_grp2.iloc[0]),
            design_old=o.design.iloc[0], design_new=n.design.iloc[0],
            genes_tested_old=int(o.n_genes_tested.iloc[0]),
            genes_tested_new=int(n.n_genes_tested.iloc[0]),
            genes_shared=len(shared),
            sig_old=int(sig_o.sum()), sig_new=int(sig_n.sum()),
            sig_both=int((sig_n & sig_o).sum()),
            sig_new_only=int((sig_n & ~sig_o).sum()), sig_old_only=int((sig_o & ~sig_n).sum()),
            log2FC_pearson=round(float(ns.log2FC.corr(os_.log2FC)), 4),
            log2FC_spearman=round(float(ns.log2FC.corr(os_.log2FC, method='spearman')), 4),
            median_abs_log2FC_diff=round(float((ns.log2FC - os_.log2FC).abs().median()), 4),
            sign_agreement_sig_either=round(float(
                ((ns.log2FC[either] > 0) == (os_.log2FC[either] > 0)).mean()), 4),
            median_SE_ratio_new_over_old=round(float((ns.SE / os_.SE).median()), 4)))
    return pd.DataFrame(rows)
