"""Reconstruct per-panel source data for Figures 1-3 and Supplementary Figure 1
of the TF x risk-gene network analysis.

Logic is transcribed from the original figure scripts so the exported numbers are
the values actually rendered:
  makefig1b.py + netfig2.key_tfs  -> Figure 1 (a,b,c,d)
  makefig3.py  + spec.build       -> Figure 2 (a,b)
  makefig4.py  + spec.build       -> Figure 3 (a-f)
  makefig2.py                     -> Supplementary Figure 1 (a,b,c)
"""
import os
import numpy as np
import pandas as pd
from scipy import stats

BASE = r'F:\Dropbox\GuoLab\Claude\TFgenes\TF_network_analysis'
DATA = os.path.join(BASE, 'data')

# ---------------------------------------------------------------- inputs
df = pd.read_csv(os.path.join(DATA, 'S1_TF_gene_correlations_all.csv'))
occ = pd.read_csv(os.path.join(DATA, 'S2_TF_leadvariant_occupancy.csv'))
e = df[df.r > 0.5]

ORD_N = ['STM', 'TAC', 'ABS', 'GOB', 'EE', 'TUF', 'CT']
ORD_A = ['STM', 'TAC', 'ABS', 'GOB', 'EE', 'TUF', 'CT', 'ASC', 'SSC']
o_all = occ.groupby('TF').n_occupied.sum()
VMAX = o_all.max()

# ---------------------------------------------------------------- spec.build
CTS = {'Normal': ['STM', 'TAC', 'ABS', 'GOB', 'EE', 'TUF', 'CT'],
       'Adenoma': ['STM', 'TAC', 'ABS', 'GOB', 'ASC', 'SSC', 'EE', 'TUF', 'CT']}
MINGENE = 10


def tau(x):
    x = np.asarray(x, float)
    return np.nan if x.max() <= 0 else np.sum(1 - x / x.max()) / (len(x) - 1)


def bh(p):
    p = np.asarray(p, float)
    o = np.argsort(p)
    q = np.empty_like(p)
    m = len(p)
    prev = 1.0
    for i in range(m - 1, -1, -1):
        prev = min(prev, p[o[i]] * m / (i + 1))
        q[o[i]] = prev
    return q


def build(tis):
    cts = [c for c in CTS[tis]
           if df[(df.tissue == tis) & (df.cell_type == c)].gene.nunique() >= MINGENE]
    d = df[(df.tissue == tis) & (df.cell_type.isin(cts))]
    ng = d.groupby('cell_type').gene.nunique()[cts]
    deg = (e[(e.tissue == tis) & (e.cell_type.isin(cts))]
           .groupby(['TF', 'cell_type']).gene.nunique().unstack())
    deg = deg.reindex(index=sorted(set(df.TF)), columns=cts).fillna(0).astype(int)
    rate = deg.div(ng, axis=1)
    rows = []
    for tf in deg.index:
        for c in cts:
            a = deg.loc[tf, c]
            b = ng[c] - a
            cc = deg.loc[tf, [x for x in cts if x != c]].sum()
            dd = ng[[x for x in cts if x != c]].sum() - cc
            p = stats.fisher_exact([[a, b], [cc, dd]], alternative='greater')[1]
            rows.append(dict(tissue=tis, TF=tf, cell_type=c, n_linked=int(a),
                             n_genes=int(ng[c]), rate=a / ng[c],
                             rate_other=cc / (cc + dd) if (cc + dd) > 0 else 0, p=p))
    R = pd.DataFrame(rows)
    R['q'] = bh(R.p.values)
    T = rate.apply(tau, axis=1)
    S = pd.DataFrame({'tau': T, 'dominant': rate.idxmax(axis=1),
                      'max_rate': rate.max(axis=1), 'total_links': deg.sum(axis=1)})
    return cts, ng, deg, rate, R, S


BUILD = {t: build(t) for t in ('Normal', 'Adenoma')}

# ---------------------------------------------------------------- Figure 1
def key_tfs(ntop=20):
    out = {}
    for tis in ('Normal', 'Adenoma'):
        d = df[(df.tissue == tis) & (df.r > 0.5)]
        deg = d.groupby('TF').gene.nunique()
        m = pd.DataFrame({'deg': deg}).join(o_all.rename('occ')).fillna(0)
        m['sc'] = (m.deg.rank(pct=True) + m.occ.rank(pct=True)) / 2
        out[tis] = list(m.sort_values('sc', ascending=False).head(ntop).index)
    return out


K = key_tfs(20)
KEY = sorted(set(K['Normal']) | set(K['Adenoma']))

sheets = {}
panel_counts = {}

for tis in ('Normal', 'Adenoma'):
    d = df[(df.tissue == tis) & (df.r > 0.5) & (df.TF.isin(KEY))]
    pairs = (d.groupby(['TF', 'gene'])
             .agg(rmax=('r', 'max'), nct=('cell_type', 'nunique')).reset_index())
    gct = df[df.tissue == tis].groupby('gene').cell_type.apply(lambda s: sorted(set(s)))
    prim = {g: (gct[g][0] if len(gct[g]) == 1 else 'Multiple') for g in pairs.gene.unique()}
    ct_list = {g: ';'.join(gct[g]) for g in pairs.gene.unique()}
    edges = pairs.rename(columns={'rmax': 'max_r', 'nct': 'n_cell_types'})
    edges['gene_node_group'] = edges.gene.map(prim)
    edges['gene_cell_types_all_r'] = edges.gene.map(ct_list)
    edges['TF_node_fill_occupancy'] = edges.TF.map(o_all).fillna(0).astype(int)
    edges.insert(0, 'tissue', tis)
    lab = 'a' if tis == 'Normal' else 'b'
    sheets[f'Fig1{lab}_edges_{tis}'] = edges.sort_values(['TF', 'gene'])
    panel_counts[f'Figure 1{lab}'] = dict(
        panel=f'Figure 1{lab}',
        tissue=('Normal colon epithelium' if tis == 'Normal'
                else 'Adenoma and serrated lesions'),
        n_links=len(pairs), n_key_TFs=pairs.TF.nunique(), n_risk_genes=pairs.gene.nunique(),
        node_order=';'.join(ORD_N if tis == 'Normal' else ORD_A),
        TF_fill_colour_max=int(VMAX))

# Fig 1c
dn = df[(df.tissue == 'Normal') & (df.r > 0.5)].groupby('TF').gene.nunique()
da = df[(df.tissue == 'Adenoma') & (df.r > 0.5)].groupby('TF').gene.nunique()
rows = []
for tis, s in (('Normal', dn), ('Adenoma', da)):
    for tf, v in s.items():
        rows.append(dict(tissue=tis, TF=tf,
                         x_lead_variants_occupied=int(o_all.get(tf, 0)),
                         y_risk_genes_r_gt_0p5=int(v),
                         labelled_in_panel=tf in ['POLR2A', 'RAD21', 'SP1', 'MAX', 'CBX3',
                                                  'GATA6', 'SMC3', 'YY1', 'CEBPB', 'FOXA2',
                                                  'TFAP4', 'ELF1'] and tis == 'Normal'))
sheets['Fig1c_scatter'] = pd.DataFrame(rows).sort_values(['tissue', 'TF'])

# Fig 1d
m = pd.DataFrame({'Normal': dn, 'Adenoma': da, 'occ': o_all}).fillna(0).loc[KEY]
m['sc'] = m.Normal.rank() + m.Adenoma.rank() + m.occ.rank()
m = m.sort_values('sc')
fig1d = m.reset_index().rename(columns={'index': 'TF', 'occ': 'lead_variants_occupied',
                                        'Normal': 'coexpressed_genes_Normal',
                                        'Adenoma': 'coexpressed_genes_Adenoma',
                                        'sc': 'rank_sum_used_for_ordering'})
fig1d.insert(1, 'y_position_bottom_to_top', np.arange(len(fig1d)))
sheets['Fig1d_dotplot'] = fig1d

# ---------------------------------------------------------------- Figure 2
ORDER2 = {'Normal': ['STM', 'TAC', 'ABS', 'GOB', 'EE', 'TUF', 'CT'],
          'Adenoma': ['STM', 'TAC', 'ABS', 'GOB', 'ASC', 'SSC']}
for lab, tis in (('a', 'Normal'), ('b', 'Adenoma')):
    cts, ng, deg, rate, R, S = BUILD[tis]
    cts_shown = [c for c in ORDER2[tis] if c in cts]
    Rs = R[R.cell_type.isin(cts_shown)]
    hits = Rs[(Rs.p < 0.05) & (Rs.n_linked >= 3)]
    best = hits.sort_values('p').drop_duplicates('TF').set_index('TF')
    tfs = sorted(sorted(set(hits.TF)),
                 key=lambda t: (cts_shown.index(best.loc[t, 'cell_type']), best.loc[t, 'p']))
    out = Rs[Rs.TF.isin(tfs) & (Rs.n_linked > 0)].copy()
    out['dominant_cell_type_for_row_order'] = out.TF.map(best.cell_type)
    out['row_position_top_to_bottom'] = out.TF.map({t: i for i, t in enumerate(tfs)})
    out['col_position_left_to_right'] = out.cell_type.map(
        {c: i for i, c in enumerate(cts_shown)})
    out['marker_size_pt2'] = 8 + 330 * out.rate
    out['colour_neglog10P_capped6'] = np.minimum(
        -np.log10(out.p.clip(lower=1e-12)), 6.0)
    out['bold_outline_FDR_lt_0p05'] = (out.q < 0.05) & (out.n_linked >= 3)
    sheets[f'Fig2{lab}_dotmatrix_{tis}'] = out.sort_values(
        ['row_position_top_to_bottom', 'col_position_left_to_right'])
    panel_counts[f'Figure 2{lab}'] = dict(
        panel=f'Figure 2{lab}',
        tissue='Normal colon' if tis == 'Normal' else 'Adenoma and serrated lesions',
        n_links=int(out.n_linked.sum()), n_key_TFs=len(tfs), n_risk_genes=np.nan,
        node_order=';'.join(cts_shown), TF_fill_colour_max=6)

# ---------------------------------------------------------------- Figure 3
MODS = [('Normal', 'TAC', 'Transit-amplifying'), ('Normal', 'ABS', 'Absorptive'),
        ('Normal', 'GOB', 'Goblet'), ('Normal', 'CT', 'Crypt-top'),
        ('Adenoma', 'ASC', 'Adenoma-specific'), ('Adenoma', 'SSC', 'Serrated-specific')]
f3_nodes, f3_edges = [], []
for k, (tis, ct, full) in enumerate(MODS):
    lab = 'abcdef'[k]
    R = BUILD[tis][4]
    s = R[(R.cell_type == ct) & (R.p < 0.05) & (R.n_linked >= 3)].sort_values('p').head(6)
    tfs = list(s.TF)
    d = e[(e.tissue == tis) & (e.cell_type == ct) & (e.TF.isin(tfs))]
    cnt = d.groupby('gene').TF.nunique()
    genes = set(cnt[cnt >= 2].index)
    if len(genes) < 5:
        genes = set(cnt.index)
    for tf in tfs:
        sub = d[d.TF == tf]
        if not set(sub.gene) & genes:
            genes |= set(sub.sort_values('r', ascending=False).gene.head(2))
    genes = sorted(genes)
    gct = df[df.tissue == tis].groupby('gene').cell_type.nunique()
    si = s.set_index('TF')
    for i, tf in enumerate(tfs):
        f3_nodes.append(dict(panel=f'Figure 3{lab}', module=full, cell_type=ct, tissue=tis,
                             node_type='TF', node=tf, order_in_panel=i,
                             n_linked=int(si.loc[tf, 'n_linked']),
                             n_genes=int(si.loc[tf, 'n_genes']),
                             rate=si.loc[tf, 'rate'], p=si.loc[tf, 'p'], q=si.loc[tf, 'q'],
                             FDR_lt_0p05=bool(si.loc[tf, 'q'] < 0.05),
                             unique_to_this_cell_type=''))
    for i, g in enumerate(genes):
        f3_nodes.append(dict(panel=f'Figure 3{lab}', module=full, cell_type=ct, tissue=tis,
                             node_type='risk_gene', node=g, order_in_panel=i,
                             n_linked=np.nan, n_genes=np.nan, rate=np.nan, p=np.nan,
                             q=np.nan, FDR_lt_0p05='',
                             unique_to_this_cell_type=bool(gct[g] == 1)))
    p3 = d[d.gene.isin(genes)][['TF', 'gene', 'r']].copy()
    p3.insert(0, 'panel', f'Figure 3{lab}')
    p3.insert(1, 'module', full)
    p3.insert(2, 'cell_type', ct)
    p3.insert(3, 'tissue', tis)
    p3['edge_style'] = np.where(p3.r > 0.7, 'thick dark (r > 0.7)',
                                'thin grey (0.5 < r <= 0.7)')
    f3_edges.append(p3.sort_values(['TF', 'gene']))
sheets['Fig3_module_nodes'] = pd.DataFrame(f3_nodes)
sheets['Fig3_module_edges'] = pd.concat(f3_edges, ignore_index=True)

# ---------------------------------------------------------- Supplementary Figure 1
tfs_all = sorted(set(df.TF))
N = e[e.tissue == 'Normal'].groupby(['TF', 'cell_type']).gene.nunique().unstack(fill_value=0)
A = e[e.tissue == 'Adenoma'].groupby(['TF', 'cell_type']).gene.nunique().unstack(fill_value=0)
O = occ.pivot(index='TF', columns='cell_type', values='n_occupied').fillna(0)
N = N.reindex(index=tfs_all, columns=ORD_N).fillna(0)
A = A.reindex(index=tfs_all, columns=ORD_A).fillna(0)
O = O.reindex(index=tfs_all, columns=ORD_A).fillna(0)
order = (N.sum(axis=1).rank() + A.sum(axis=1).rank()
         + O.sum(axis=1).rank()).sort_values(ascending=False).index
N, A, O = N.loc[order], A.loc[order], O.loc[order]
rowpos = {t: i for i, t in enumerate(order)}
vmaxC = max(N.values.max(), A.values.max())
for lab, M, what, vmx in (('a', N, 'co-expressed risk genes, normal colon', vmaxC),
                          ('b', A, 'co-expressed risk genes, adenoma / serrated', vmaxC),
                          ('c', O, 'risk-gene lead variants occupied by TF', O.values.max())):
    long = M.stack().rename('value').reset_index()
    long.columns = ['TF', 'cell_type', 'value']
    long['value'] = long.value.astype(int)
    long.insert(0, 'panel', f'Supplementary Figure 1{lab}')
    long.insert(1, 'quantity', what)
    long['row_position_top_to_bottom'] = long.TF.map(rowpos)
    long['col_position_left_to_right'] = long.cell_type.map(
        {c: i for i, c in enumerate(M.columns)})
    long['colour_scale_max'] = int(vmx)
    sheets[f'FigS1{lab}_matrix'] = long.sort_values(
        ['row_position_top_to_bottom', 'col_position_left_to_right'])

sheets['Fig1_panel_summary'] = pd.DataFrame(panel_counts.values())

# ---------------------------------------------------------------- write
README = pd.DataFrame([
    ('Figure 1a', 'Fig1a_edges_Normal',
     'TF-risk gene links drawn in the normal-colon network (r > 0.5, key TFs). '
     'max_r = strongest correlation across cell types; n_cell_types = cell types supporting '
     'the link; gene_node_group = sector the gene node sits in; TF_node_fill_occupancy = TF fill value.'),
    ('Figure 1b', 'Fig1b_edges_Adenoma', 'As Figure 1a, adenoma and serrated lesions.'),
    ('Figure 1c', 'Fig1c_scatter',
     'One point per TF per tissue: x = risk-gene lead variants occupied, '
     'y = risk genes co-expressed at r > 0.5.'),
    ('Figure 1d', 'Fig1d_dotplot',
     'Key TFs, one row each; three plotted values per row plus the rank-sum used for row order.'),
    ('Figure 1 (all panels)', 'Fig1_panel_summary',
     'Link / TF / gene counts quoted in the panel subtitles, and axis orders.'),
    ('Figure 2a', 'Fig2a_dotmatrix_Normal',
     'One row per plotted dot. Dot area from rate (8 + 330*rate pt2), fill from '
     '-log10(P) capped at 6, bold outline when FDR < 0.05 and n_linked >= 3. '
     'P from one-sided Fisher exact test, cell type vs all other cell types of the same tissue; '
     'q = Benjamini-Hochberg across all TF x cell-type tests of that tissue.'),
    ('Figure 2b', 'Fig2b_dotmatrix_Adenoma', 'As Figure 2a, adenoma and serrated lesions.'),
    ('Figure 3a-f', 'Fig3_module_nodes',
     'Nodes of the six module panels: module TFs (up to 6, lowest Fisher P, n_linked >= 3) '
     'with their enrichment statistics, and the risk-gene nodes with their exclusivity flag.'),
    ('Figure 3a-f', 'Fig3_module_edges',
     'TF-risk gene edges drawn in each module panel, with r and the line style it maps to.'),
    ('Supplementary Figure 1a', 'FigS1a_matrix',
     'Heat-map cell values: risk genes co-expressed with each TF per cell type, normal colon. '
     'Row order is shared by all three panels (sum of the three row rank scores, descending).'),
    ('Supplementary Figure 1b', 'FigS1b_matrix', 'As panel a, adenoma and serrated lesions.'),
    ('Supplementary Figure 1c', 'FigS1c_matrix',
     'Heat-map cell values: risk-gene lead variants occupied by each TF per cell type.'),
], columns=['Figure panel', 'Sheet', 'Contents'])

src = pd.DataFrame([
    ('S1_TF_gene_correlations_all.csv', 'TF-risk gene Spearman/Pearson r per tissue and cell type'),
    ('S2_TF_leadvariant_occupancy.csv', 'ChIP-seq occupancy of risk-gene lead variants per TF and cell type'),
], columns=['Input file (data/)', 'Description'])

order_sheets = (['README', 'Fig1_panel_summary']
                + [k for k in sheets if k.startswith('Fig1') and k != 'Fig1_panel_summary']
                + [k for k in sheets if k.startswith('Fig2')]
                + [k for k in sheets if k.startswith('Fig3')]
                + [k for k in sheets if k.startswith('FigS1')])

OUT = 'SourceData_Figures1-3_S1.xlsx'
with pd.ExcelWriter(OUT, engine='openpyxl') as xw:
    README.to_excel(xw, sheet_name='README', index=False, startrow=0)
    src.to_excel(xw, sheet_name='README', index=False, startrow=len(README) + 3)
    for k in order_sheets:
        if k == 'README':
            continue
        sheets[k].to_excel(xw, sheet_name=k[:31], index=False)

os.makedirs('sourcedata_csv', exist_ok=True)
for k, v in sheets.items():
    v.to_csv(os.path.join('sourcedata_csv', k + '.csv'), index=False)

print('KEY TFs (n=%d):' % len(KEY), ','.join(KEY))
print('sheets:', [(k, sheets[k].shape) for k in order_sheets if k != 'README'])
