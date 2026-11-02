#!/usr/bin/env python
"""
scTWAS analytic workflow -- TWO figures, one script:
    Fig_scTWAS_workflow_NormalColon.pdf   (7 epithelial cell types)
    Fig_scTWAS_workflow_ColonPolyps.pdf   (9 epithelial populations)

Each figure runs one lineage end to end, from the scRNA-seq / genotype
resources and the cell-type scATAC-seq peak sets down to the susceptibility
genes, and names every cell type explicitly.

Style follows 'Supplementary Figure 2.pdf' (the sTF-TWAS workflow) with two
deliberate departures: NO filled boxes, and short labels only.
    dashed grey outline  = data resource
    solid teal outline   = analysis step   (navy + thicker = core model step)
    thin black arrows down a central spine
Arial 11 pt throughout.

Content: MethodsResults.docx (Results section 1; Methods).
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyArrow

# ---------------------------------------------------------------- style
TEAL    = (0.082353, 0.376471, 0.509804)
NAVY    = (0.015686, 0.141176, 0.200000)
MAGENTA = (0.847059, 0.431373, 0.800000)
GREY    = (0.30, 0.30, 0.30)

for fam in ("Arial", "Helvetica", "DejaVu Sans"):
    if any(fam in f.name for f in matplotlib.font_manager.fontManager.ttflist):
        FONT = fam
        break
plt.rcParams.update({"font.family": FONT, "pdf.fonttype": 42, "ps.fonttype": 42})

W, H = 792.0, 612.0      # points (11 x 8.5 in), as in the reference PDF
FS   = 11.0              # body and header text: Arial 11
LH   = 15.0              # line pitch

SPEC = {
    "NormalColon": dict(
        title="Normal colon",
        cts=["ABS", "CT", "EE", "GOB", "STM", "TAC", "TUF"],
        rna_head="scRNA-seq + genotype   (N = 381)",
        rna=["VUMC HTAN + Gut Cell Atlas (n = 103)",
             "Wellcome Sanger IBDverse (n = 278)"],
        atac=["VUMC  (n = 31 normal colon)",
              "Stanford HTAN  (n = 8 normal colon)"],
        models="7 cell-type-specific models",
        # Normal colon draws on two independent scRNA-seq + genotype panels
        # (VUMC n = 103, Wellcome Sanger n = 278; 103 + 278 = 381). Each panel
        # yields its own TWAS P per gene x cell type; the Cauchy combination
        # collapses the pair before Bonferroni thresholding.
        assoc_extra=["Cauchy combination of $P$ values:  "
                     "VUMC HTAN + Gut Cell Atlas  and  Wellcome Sanger"],
    ),
    "ColonPolyps": dict(
        title="Colorectal polyps",
        cts=["ABS", "CT", "EE", "GOB", "STM", "TAC", "TUF", "ASC", "SSC"],
        rna_head="scRNA-seq + genotype   (N = 84)",
        rna=["VUMC HTAN COLON MAP",
             "adenoma and serrated polyps"],
        atac=["VUMC  (n = 23 polyps)",
              "Stanford HTAN  (n = 48 polyps)"],
        models="9 cell-type-specific models",
    ),
}


def build(key, spec):
    fig = plt.figure(figsize=(W / 72.0, H / 72.0), facecolor="white")
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, W); ax.set_ylim(H, 0)          # y downward, PDF-like
    ax.axis("off")

    # ---- primitives ----------------------------------------------------
    def data_box(x0, y0, x1, header, bullets, hc=TEAL, pad=7.0):
        hh = LH + 4.0
        h = pad + hh + LH * len(bullets) + pad
        ax.add_patch(Rectangle((x0, y0), x1 - x0, h, fill=False, ec=GREY,
                               lw=0.8, ls=(0, (3, 2)), zorder=2))
        ax.text((x0 + x1) / 2, y0 + pad + hh / 2 - 1.0, header, ha="center",
                va="center", color=hc, fontsize=FS, fontweight="bold",
                zorder=4)
        ax.plot([x0 + 6, x1 - 6], [y0 + pad + hh, y0 + pad + hh], color=hc,
                lw=0.9, zorder=3)
        ty = y0 + pad + hh + LH / 2 + 2.0
        for b in bullets:
            ax.text(x0 + 12, ty, "\u2022  " + b, ha="left", va="center",
                    fontsize=FS, color="black", zorder=4)
            ty += LH
        return y0 + h

    def band(x0, y0, x1, lines, emph=False, pad=8.0, heights=None):
        hs = heights or [LH] * len(lines)
        h = 2 * pad + sum(hs)
        col = NAVY if emph else TEAL
        ax.add_patch(Rectangle((x0, y0), x1 - x0, h, fill=False, ec=col,
                               lw=1.6 if emph else 0.9, zorder=3))
        ty = y0 + pad
        for i, (t, lh) in enumerate(zip(lines, hs)):
            ax.text((x0 + x1) / 2, ty + lh / 2, t, ha="center", va="center",
                    color=col if i == 0 else "black", fontsize=FS,
                    fontweight="bold" if i == 0 else "normal", zorder=4)
            ty += lh
        return y0 + h

    def arrow(x, y0, y1, width=0.9):
        ax.add_patch(FancyArrow(x, y0, 0, y1 - y0 - 5.0, width=width,
                                head_width=6.5, head_length=5.6,
                                length_includes_head=False, fc="black",
                                ec="black", lw=0, zorder=2))

    def elbow(x_from, y_from, x_to, y_join, y_to):
        ax.plot([x_from, x_from], [y_from, y_join], color="black", lw=0.9,
                zorder=2)
        ax.plot([x_from, x_to], [y_join, y_join], color="black", lw=0.9,
                zorder=2)
        arrow(x_to, y_join, y_to)

    # ---- title ---------------------------------------------------------
    ax.text(40, 26, "scTWAS \u2014 " + spec["title"], ha="left", va="center",
            fontsize=14.5, fontweight="bold", color=NAVY)

    # ---- row 1: the two input resources --------------------------------
    b = data_box(40, 44, 392, spec["rna_head"], spec["rna"])
    data_box(400, 44, 752, "scATAC-seq open chromatin", spec["atac"],
             hc=MAGENTA)

    # ---- row 2: the cell types, one tile each --------------------------
    y = b + 32.0
    rail = y - 11.0
    ax.plot([216, 216], [b, rail], color="black", lw=0.9, zorder=2)
    ax.plot([576, 576], [b, rail], color="black", lw=0.9, zorder=2)
    cts = spec["cts"]
    x0, x1, gap = 40.0, 752.0, 7.0
    tw = (x1 - x0 - gap * (len(cts) - 1)) / len(cts)
    th = 28.0
    ax.plot([x0 + tw / 2, x1 - tw / 2], [rail, rail], color="black", lw=0.9,
            zorder=2)
    for i, ct in enumerate(cts):
        bx = x0 + i * (tw + gap)
        arrow(bx + tw / 2, rail, y + 0.4, width=0.7)
        ax.add_patch(Rectangle((bx, y), tw, th, fill=False, ec=GREY, lw=0.8,
                               ls=(0, (3, 2)), zorder=2))
        ax.text(bx + tw / 2, y + th / 2, ct, ha="center", va="center",
                fontsize=FS, fontweight="bold", color=TEAL, zorder=4)
    y_end = y + th

    # ---- row 3: variant prioritization ---------------------------------
    y = y_end + 24.0
    arrow(396, y_end, y)
    y_end = band(112, y, 680, [
        "Variants in cell-type open chromatin",
        "MAF > 0.05    \u00b7    HWE $P$ > 1 \u00d7 10$^{-6}$    \u00b7    \u00b1 1 Mb of TSS",
    ])

    # ---- row 4: GReX model construction --------------------------------
    y = y_end + 22.0
    arrow(396, y_end, y)
    y_end = band(112, y, 680, [
        "Elastic-net GReX models   (sTF-TWAS, \u03b1 = 0.5)",
        spec["models"],
    ], emph=True)

    # ---- row 5: the two GWAS resources ---------------------------------
    y_model = y_end
    y = y_end + 21.0
    data_box(104, y, 384, "CRC GWAS \u2014 European",
             ["99,260 cases / 1,119,512 controls"])
    y_end = data_box(408, y, 688, "CRC GWAS \u2014 multi-ancestry",
                     ["132,261 cases / 1,752,398 controls"])

    # ---- row 6: association test ---------------------------------------
    y = y_end + 26.0
    # GReX models feed the association step down the central spine, which
    # passes between the two GWAS boxes; each GWAS box joins that spine.
    for xb in (244, 548):
        ax.plot([xb, xb], [y_end, y_end + 13.0], color="black", lw=0.9,
                zorder=2)
        ax.plot([xb, 396], [y_end + 13.0, y_end + 13.0], color="black",
                lw=0.9, zorder=2)
    arrow(396, y_model, y)
    extra = spec.get("assoc_extra", [])
    y_end = band(112, y, 680,
                 ["scTWAS association",
                  r"$Z_g \;\approx \sum_{l\,\in\, \mathrm{Model}_g} \omega_{lg}\,"
                  r"\frac{\sigma_l}{\sigma_g}\,\frac{\beta_l}{se(\beta_l)}$"]
                 + extra
                 + ["Bonferroni-corrected $P$ < 0.05"],
                 heights=[LH, 34.0] + [LH] * (len(extra) + 1))

    # ---- row 7: shared-signal evidence ---------------------------------
    y = y_end + 18.0
    arrow(396, y_end, y)
    y_end = band(178, y, 614, [
        "Colocalization  +  SMR\u2013HEIDI",
    ])

    # ---- row 8: output --------------------------------------------------
    y = y_end + 18.0
    arrow(396, y_end, y)
    y_end = band(178, y, 614, [
        "CRC susceptibility genes",
    ])

    out = "Fig_scTWAS_workflow_%s" % key
    fig.savefig(out + ".pdf")
    fig.savefig(out + ".png", dpi=400)
    plt.close(fig)
    return out, round(y_end, 1)


for k, s in SPEC.items():
    print(*build(k, s))
