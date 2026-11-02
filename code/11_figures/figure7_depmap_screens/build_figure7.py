
import numpy as np, pandas as pd, matplotlib as mpl, matplotlib.pyplot as plt

ONC, ONC_M, ONC_L = '#2c6d5b', '#43a389', '#8fcdbc'
SUP, SUP_M = '#2c566d', '#4380a3'
GREY = '#8a8880'

def build_figure7(long24, order24, vals, ORDER2, PAL2, PW, shown, topC, cell, LINES,
                  out_png, out_pdf, figsize=(6.0, 6.0)):
    F = {'title':6.2, 'axis':5.2, 'tick':4.6, 'gene':4.3, 'cell':4.0, 'note':3.8, 'letter':7.5}
    fig = plt.figure(figsize=figsize)
    gs = fig.add_gridspec(2, 2, height_ratios=[5.55, 2.05], width_ratios=[1.0, 1.12],
                          left=0.098, right=0.995, top=0.955, bottom=0.115,
                          wspace=0.26, hspace=0.46)

    # ================= (a) gene-set essentiality, top left =================
    axa = fig.add_subplot(gs[0, 0])
    CAP, TOP = -2.0, 1.0
    data = [vals[k] for k in ORDER2]
    clipped = {k: int((vals[k] < CAP).sum()) for k in ORDER2}
    vp = axa.violinplot([np.clip(v, CAP, None) for v in data], positions=np.arange(4),
                        widths=0.80, showextrema=False, showmedians=False)
    for body, k in zip(vp['bodies'], ORDER2):
        body.set_facecolor(PAL2[k]); body.set_alpha(0.35); body.set_edgecolor('none')
    bx = axa.boxplot(data, positions=np.arange(4), widths=0.28, patch_artist=True, showfliers=False,
                     medianprops=dict(color='white', lw=0.9), whiskerprops=dict(color=GREY, lw=0.6),
                     capprops=dict(color=GREY, lw=0.6), boxprops=dict(lw=0.4, edgecolor=GREY))
    for p, k in zip(bx['boxes'], ORDER2): p.set_facecolor(PAL2[k])
    rng = np.random.default_rng(1)
    for i, k in enumerate(ORDER2):
        v = vals[k]; inside = v[v >= CAP]; below = v[v < CAP]
        axa.scatter(i + rng.uniform(-0.32, 0.32, len(inside)), inside, s=2.0, color=PAL2[k],
                    alpha=0.45, linewidths=0, zorder=1)
        if len(below):
            axa.scatter(i + rng.uniform(-0.13, 0.13, len(below)), np.full(len(below), CAP+0.03),
                        s=6, marker='v', color=PAL2[k], edgecolor='#52514e', linewidths=0.25, zorder=4)
    axa.axhline(0, color=GREY, lw=0.6)
    axa.axhline(-0.5, color='#8f6634', lw=0.7, ls='--')
    P = PW.set_index(['group_A','group_B'])
    base, step = 0.30, 0.215
    for lvl, (a_, b_) in enumerate(shown):
        i, j = ORDER2.index(a_), ORDER2.index(b_)
        y = base + step*lvl
        r = P.loc[(a_, b_)]
        col = ONC if r['significant_Holm_0.05'] else GREY
        axa.plot([i, i, j, j], [y-0.05, y, y, y-0.05], color=col, lw=0.7)
        lab = f"P = {r.P_two_sided_Holm:.0e}" if r.P_two_sided_Holm < 1e-3 else f"P = {r.P_two_sided_Holm:.3f}"
        axa.text((i+j)/2, y+0.012, lab, ha='center', va='bottom', fontsize=F['note'], color=col)
    axa.set_ylim(CAP, TOP)
    axa.set_xticks(np.arange(4))
    axa.set_xticklabels([f"{lab}\nn={len(vals[k])}" for lab, k in
                         zip(['Random\nbackground','CRC\nconsensus','scTWAS\nrisk genes','CRC\ndriver'], ORDER2)],
                        fontsize=F['tick'])
    axa.set_ylabel('DepMap gene effect, median across 63 lines', fontsize=F['axis'])
    axa.set_title('Gene-set dependency enrichment', loc='left', fontsize=F['title'], pad=11)
    axa.tick_params(labelsize=F['tick'])
    axa.spines[['top','right']].set_visible(False)
    axa.text(0.0, -0.155, 'Wilcoxon rank-sum, Holm-adjusted over all 6 pairs; axis capped at −2 ('
             + '\u25be = ' + str(sum(clipped.values())) + ' genes below)',
             transform=axa.transAxes, fontsize=F['note'], va='top', color='#52514e')

    # ================= (c) CRISPR matrix, top right =================
    axc = fig.add_subplot(gs[0, 1])
    DIV = mpl.colors.LinearSegmentedColormap.from_list('div', [ONC, ONC_M, '#f2f1ed', SUP_M, SUP])
    genes, N, nc, lim = list(topC.index), len(topC), len(LINES), 25.0
    M = np.array([[cell[(g, L)][0] for L in LINES] for g in genes], float)
    S = np.array([[cell[(g, L)][1] for L in LINES] for g in genes], bool)
    im = axc.imshow(M, cmap=DIV, norm=mpl.colors.TwoSlopeNorm(vmin=-lim, vcenter=0, vmax=lim),
                    aspect='auto', extent=(-0.5, nc-0.5, N-0.5, -0.5))
    for i in range(N):
        for j in range(nc):
            if S[i, j]:
                axc.add_patch(mpl.patches.Rectangle((j-0.5, i-0.5), 1, 1, fill=False,
                              edgecolor='black', lw=0.9, zorder=3))
                axc.text(j, i, f"{M[i,j]:.0f}", ha='center', va='center', fontsize=F['cell'],
                         color='white' if abs(M[i,j]) > lim*0.6 else '#0b0b0b', zorder=4)
    x0 = nc - 0.5
    for i, g in enumerate(genes):
        z = float(topC.loc[g,'rep_Z'])
        axc.add_patch(mpl.patches.Rectangle((x0+0.10, i-0.5), 0.80, 1, facecolor=ONC_M if z > 0 else SUP_M,
                      alpha=0.8, edgecolor='white', lw=0.4, clip_on=False, zorder=5))
        axc.text(x0+0.50, i, f"{z:+.1f}", ha='center', va='center', fontsize=F['cell'], zorder=6)
    axc.set_xlim(-0.5, nc+0.95); axc.set_ylim(N-0.5, -0.5)
    axc.set_xticks(range(nc)); axc.set_xticklabels(LINES, fontsize=F['tick'])
    axc.set_yticks(range(N)); axc.set_yticklabels(genes, style='italic', fontsize=F['gene'])
    axc.text(x0+0.50, -0.72, 'CRC risk\nscore', ha='center', va='bottom', fontsize=F['note'],
             color='#52514e', linespacing=0.95)
    axc.set_title('CRISPR-Cas9 knockout of the 24 risk-increasing genes', loc='left', fontsize=F['title'], pad=11)
    for s in axc.spines.values(): s.set_visible(False)
    axc.tick_params(length=0, labelsize=F['tick'])
    cax = axc.inset_axes([0.0, -0.075, 0.62, 0.014])
    cax.imshow(np.linspace(-lim, lim, 256).reshape(1, -1), cmap=DIV, norm=im.norm, aspect='auto')
    cax.set_yticks([]); cax.set_xticks([0, 127.5, 255])
    cax.set_xticklabels(['−25', '0', '+25'], fontsize=F['note'])
    cax.tick_params(length=1.2, pad=1.0)
    for sp in cax.spines.values(): sp.set_linewidth(0.4)
    cax.set_xlabel('Δ viability vs control (points); outline P < 0.05\n'
                   'green = oncogene-like, blue = suppressor-like',
                   fontsize=F['note'], labelpad=1.2, linespacing=1.15)

    # ================= (b) 24 dependency genes, bottom =================
    axb = fig.add_subplot(gs[1, :])
    d24 = [long24.loc[long24.gene == g, 'gene_effect'].values for g in order24]
    bp = axb.boxplot(d24, positions=np.arange(len(order24)), widths=0.66, patch_artist=True,
                     showfliers=False, medianprops=dict(color='black', lw=0.8),
                     whiskerprops=dict(color=GREY, lw=0.5), capprops=dict(color=GREY, lw=0.5),
                     boxprops=dict(lw=0.35, edgecolor=GREY))
    for p in bp['boxes']: p.set_facecolor(ONC_M); p.set_alpha(0.55)
    rng2 = np.random.default_rng(0)
    for i, v in enumerate(d24):
        axb.scatter(i + rng2.uniform(-0.19, 0.19, len(v)), v, s=1.4, color=ONC_M,
                    alpha=0.5, linewidths=0, zorder=3)
    axb.axhline(0, color=GREY, lw=0.6)
    axb.axhline(-0.5, color='#8f6634', lw=0.7, ls='--')
    axb.set_xticks(np.arange(len(order24)))
    axb.set_xticklabels(order24, rotation=90, style='italic', fontsize=F['gene'])
    axb.set_xlim(-0.7, len(order24)-0.3)
    axb.set_ylabel('DepMap gene effect', fontsize=F['axis'])
    axb.set_title('24 scTWAS genes are pan-essential dependencies in 63 colorectal lines',
                  loc='left', fontsize=F['title'])
    axb.tick_params(labelsize=F['tick'])
    axb.spines[['top','right']].set_visible(False)

    for ax, letter, dx, dy in [(axa,'a',-0.125,1.045), (axc,'c',-0.155,1.045), (axb,'b',-0.062,1.06)]:
        ax.text(dx, dy, letter, transform=ax.transAxes, fontsize=F['letter'],
                fontweight='bold', va='bottom', ha='left')
    fig.savefig(out_png, dpi=450)
    fig.savefig(out_pdf)
    return fig, clipped
