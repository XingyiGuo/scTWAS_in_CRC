#!/usr/bin/env python3
"""11 - Emit the numbers the report needs, so the document never restates a value by hand."""
import pandas as pd, numpy as np, json

risk = pd.read_csv("results/risk_genes_harmonized.csv")
g = pd.read_csv("results/integrated_gene_priority.csv")
drg = pd.read_csv("results/gene_druggability.csv")
rev = pd.read_csv("results/signature_reversal_all.csv.gz")
summ = pd.read_csv("results/signature_reversal_summary.csv")
rec = pd.read_csv("results/signature_reversal_recurrent.csv")
prox = pd.read_csv("results/network_proximity_significant.csv")
proxall = pd.read_csv("results/network_proximity.csv.gz")
gd = pd.read_csv("results/gene_drug_candidates.csv")
enr = pd.read_csv("results/pathway_enrichment.csv")
pc = pd.read_csv("results/reversal_positive_control.csv").reset_index(drop=True)
pck = pd.read_csv("results/reversal_positive_control_knowndrugs.csv")
tra = pd.read_csv("results/gene_tractability_score.csv")
mod = pd.read_csv("results/module_knockout_impact.csv")

D = {}
D["n_assoc"] = len(risk)
D["n_genes"] = int(risk["gene"].nunique())
D["n_mapped"] = int(g["gene"].isin(tra["gene"]).sum())
D["n_normal_assoc"] = int((risk.tissue == "Normal").sum())
D["n_adenoma_assoc"] = int((risk.tissue == "Adenoma").sum())
D["n_normal_ct"] = int(risk[risk.tissue == "Normal"].cell_type.nunique())
D["n_adenoma_ct"] = int(risk[risk.tissue == "Adenoma"].cell_type.nunique())
D["n_shared_genes"] = int((risk.groupby("gene")["tissue"].nunique() == 2).sum())
D["n_risk_incr"] = int((risk.direction == "risk-increasing").sum())
D["n_protective"] = int((risk.direction == "protective").sum())
D["n_in_network"] = int(tra.shape[0])
D["n_lincs_kd"] = 90
D["n_lincs_oe"] = 49

for k in ["ko_impact_validation", "ko_validation_controlled", "ko_validation_summary",
          "dti_cv_performance", "gdsc_proximity_validation", "tractability_validation",
          "ko_signature_validation"]:
    D[k] = json.load(open(f"results/{k}.json"))

D["tiers"] = drg["tier"].value_counts().to_dict()
D["tier1_genes"] = drg[drg.tier.str.startswith("T1")][
    ["gene", "direction", "n_ct", "min_p", "matched_approved_drugs"]].to_dict("records")
D["n_screens"] = int(summ.shape[0])
D["n_candidates"] = int(rev["candidate"].sum())
D["cand_per_ct"] = summ.sort_values("n_candidates", ascending=False).to_dict("records")
D["n_pc_compounds"] = int(len(pc))
D["pc_known"] = pck[["rank", "drug_name", "connectivity_score", "fdr_query"]].to_dict("records")
D["n_prox_pairs"] = int(len(prox))
D["n_prox_compounds"] = int(prox.drug_name.nunique())
D["n_prox_tested"] = int(len(proxall))
ap = prox[prox.groups.astype(str).str.contains("approved")]
D["n_prox_approved"] = int(ap.drug_name.nunique())
D["prox_top"] = (ap.groupby(["drug_name", "targets"])
                 .agg(n_celltypes=("cell_type", "size"), mean_z=("z_proximity", "mean"))
                 .reset_index().sort_values(["n_celltypes", "mean_z"])
                 .sort_values("n_celltypes", ascending=False).head(10).to_dict("records"))
D["rec_top"] = rec[rec.groups.astype(str).str.contains("approved")].head(10)[
    ["drug_name", "groups", "n_celltypes", "mean_reversal", "targets"]].to_dict("records")
D["top_targets"] = g.head(15)[["gene", "direction", "n_ct", "min_p", "tier", "ko_impact",
                               "tractability_z", "priority_score",
                               "matched_approved_drugs"]].to_dict("records")
D["gene_drug"] = gd[gd.evidence_layers >= 3][
    ["gene", "direction", "drug_name", "groups", "action_class", "reversal_candidate",
     "proximal", "n_ct", "min_p"]].to_dict("records")
D["n_gene_drug_pairs"] = int(len(gd))
D["n_direction_matched"] = int(gd.direction_match.sum())
D["n_pathway_sig"] = int((enr.fdr < 0.1).sum())
D["pathways_top"] = (enr[enr.fdr < 0.1].sort_values("p").groupby(["tissue", "cell_type"])
                     .head(2).sort_values(["tissue", "cell_type", "p"])
                     [["tissue", "cell_type", "pathway", "source", "overlap", "fdr", "genes"]]
                     .head(18).to_dict("records"))
D["erbb2"] = g[g.gene == "ERBB2"].iloc[0].to_dict()
D["erbb2"]["min_p"] = float(drg.loc[drg.gene == "ERBB2", "min_p"].iat[0])
D["mod_top"] = mod.sort_values("tv_z", ascending=False).drop_duplicates("gene").head(10)[
    ["gene", "tissue", "cell_type", "tv_z"]].to_dict("records")
D["n_module_tests"] = int(len(mod))


# ---------------------------------------------------------------- extended knockout
spec = pd.read_csv("results/ko_impact_spectrum_full.csv")
dsj = json.load(open("results/ko_driver_spectrum.json"))
dj = json.load(open("results/ko_driver_analysis.json"))
conv = json.load(open("results/ko_convergence.json"))
mc = pd.read_csv("results/ko_module_convergence.csv")
ctrlj = json.load(open("results/ko_downstream_control.json"))
recur_pw = pd.read_csv("results/ko_downstream_pathways_recurrent.csv")
enr_ko = pd.read_csv("results/ko_downstream_pathways.csv")
lj = json.load(open("results/lineage_network.json"))
lz = pd.read_csv("results/lineage_network_proximity.csv")
sej = json.load(open("results/stage_epistasis.json"))
ecj = json.load(open("results/epistasis_control.json"))
mpp = pd.read_csv("results/module_pathway_proximity.csv")
dprox = pd.read_csv("results/risk_gene_driver_proximity.csv")
drvp = pd.read_csv("results/driver_perturbation_by_ko.csv")
stg = pd.read_csv("results/stage_contrast.csv")

D["spectrum"] = dsj["spectrum"]
D["spectrum_cmp"] = dsj["comparisons"]
D["driver_overlap"] = dsj["driver_overlap"]
D["driver_neighbours"] = dsj["direct_neighbours"]
D["driver_proximity_overall"] = dj["driver_proximity"]
D["driver_closest"] = (dprox[~dprox["is_driver"]].nsmallest(10, "driver_proximity_z")
                       [["gene", "direction", "tissues", "n_celltype_assoc",
                         "dist_to_nearest_driver", "driver_proximity_z"]].to_dict("records"))
D["driver_perturb_top"] = drvp.head(6)[["driver", "frac_risk", "frac_bg", "odds_ratio", "p", "fdr"]].to_dict("records")
D["n_driver_sig"] = int((drvp["fdr"] < 0.1).sum())
D["n_drivers_tested"] = int(len(drvp))

mpiv = mpp.pivot_table(index="tissue", columns="pathway", values="z")
D["tissue_pathway"] = mpiv.round(3).to_dict()
D["module_pathway_sig"] = (mpp[mpp["z"] < -1.645].sort_values("z")
                           [["tissue", "cell_type", "pathway", "z"]].head(12).to_dict("records"))
D["n_module_pathway_sig"] = int((mpp["z"] < -1.645).sum())
D["n_module_pathway_tests"] = int(len(mpp))

D["ko_downstream"] = {
    "n_enrichments": int(len(enr_ko)), "n_genes": int(enr_ko["gene"].nunique()),
    "top": recur_pw.head(8)[["pathway", "source", "n_genes"]].to_dict("records"),
    "control": ctrlj,
}
D["convergence"] = conv
D["module_convergence"] = mc.to_dict("records")
D["lineage"] = lj
D["lineage_sig"] = (lz[lz["z"] < -1.645].sort_values("z")
                    [["tissue", "cell_type", "programme", "z", "is_own"]].to_dict("records"))
D["stage"] = {k: v for k, v in sej.items() if k not in ("epistasis", "interacting_pairs")}
D["stage_counts"] = stg["stage"].value_counts().to_dict()
D["epistasis"] = ecj


# ---------------------------------------------------------------- spatial omics
sde = pd.read_csv("results/spatial_risk_gene_DE.csv")
scon = pd.read_csv("results/spatial_direction_concordance.csv")
svb = pd.read_csv("results/spatial_risk_vs_background.csv")
srep = pd.read_csv("results/spatial_replication.csv")
scos = json.load(open("results/cosmx_summary.json"))
sct = pd.read_csv("results/cosmx_celltype_trend.csv")

epi = sde[(sde.experiment == "exp2") & (sde.compartment == "PanCK+")]
# significance among the tested risk genes, not genome-wide (see 20_spatial_de.py)
sig = epi[epi["sig"]]
sig_gw = epi[epi["FDR"] < 0.05]
D["spatial"] = {
    "n_risk_measured_exp2": int(epi.gene.nunique()),
    "n_risk_measured_exp1": int(sde[(sde.experiment == "exp1")].gene.nunique()),
    "n_sig_any": int(sig.gene.nunique()),
    "n_sig_any_genomewide": int(sig_gw.gene.nunique()),
    "sig_per_transition_genomewide": {t: int((sig_gw.transition == t).sum())
                                      for t in sig_gw.transition.unique()},
    "sig_per_transition": {t: int((sig.transition == t).sum()) for t in sig.transition.unique()},
    "verdicts": sig.groupby(["transition", "verdict"]).size().unstack(fill_value=0).to_dict(),
    "concordance_epi": scon[(scon.experiment == "exp2") & (scon.compartment == "PanCK+")]
        [["transition", "n_genes", "n_concordant", "frac_concordant", "p_binomial",
          "p_direction_permutation", "p_vs_background"]].to_dict("records"),
    "vs_background_epi": svb[(svb.experiment == "exp2") & (svb.compartment == "PanCK+")]
        [["transition", "frac_risk_sig", "frac_other_sig", "odds_ratio", "p_fisher"]].to_dict("records"),
    "replication": srep[srep.compartment == "PanCK+"]
        [["transition", "rho_all", "rho_risk", "n_risk_shared", "sign_agreement_risk"]].to_dict("records"),
    "top_genes": {t: sig[sig.transition == t].nsmallest(10, "FDR")
                  [["gene", "twas_direction", "log2FC", "FDR", "verdict", "cell_types"]]
                  .to_dict("records")
                  for t in ["NOR->ADE", "ADE->CORE", "NOR->CORE"]},
    "cosmx": {"n_cells": scos["n_cells"], "n_patients": scos["n_patients"],
              "stages": scos["stages"],
              "top_trend": sct.head(6)[["cellType", "mean_prop_im", "mean_prop_T1sm1",
                                        "mean_prop_T1sm3", "spearman_rho", "p_trend",
                                        "fdr_trend"]].to_dict("records"),
              "celltypes": scos["cosmx_celltypes"]},
}

# ---------------------------------------------------------------- validated gene set
sval = json.load(open("results/validated_summary.json"))
svf = json.load(open("results/validated_vs_failed.json"))
spc = json.load(open("results/cosmx_percell_summary.json"))
vpri = pd.read_csv("results/validated_priority.csv")
vct = pd.read_csv("results/validated_by_celltype.csv")
vall = pd.read_csv("results/validated_risk_genes_all.csv")
cxb = pd.read_csv("results/cosmx_crc_vs_normal_epi.csv")

D["validated"] = {
    "rule": sval["rule"],
    "n_measured": sval["n_measured"], "n_validated": sval["n_validated"],
    "frac_validated": sval["frac_validated"],
    "n_oncogene_like": sval["n_oncogene_like"], "n_tsg_like": sval["n_tsg_like"],
    "n_monotone": sval["n_monotone"],
    "n_never_significant": sval["n_never_significant"],
    "n_sig_opposite_only": sval["n_sig_opposite_only"],
    "n_excluded_by_contradiction": sval["n_excluded_by_contradiction"],
    "permutation": sval["permutation"],
    "lead_transition_counts": sval["lead_transition_counts"],
    "by_celltype": vct.to_dict("records"),
    "top": vpri.head(20).to_dict("records"),
    "vs_failed": svf,
    "erbb2": vall[vall["gene"] == "ERBB2"][
        ["gene", "validated", "log2FC_NOR_CORE", "FDR_NOR_CORE", "twas_direction"]
    ].to_dict("records"),
    "znrf3": vall[vall["gene"] == "ZNRF3"][
        ["gene", "validated", "log2FC_NOR_CORE", "FDR_NOR_CORE"]].to_dict("records"),
    "ctnnb1": vall[vall["gene"] == "CTNNB1"][
        ["gene", "validated", "log2FC_NOR_CORE", "FDR_NOR_CORE"]].to_dict("records"),
}

D["cosmx_percell"] = {
    "n_genes_on_panel": spc["n_genes_on_panel"],
    "n_validated_on_panel": spc["crc_vs_normal"]["n_validated_on_panel"],
    "n_patients": spc["n_patients"],
    "crc_vs_normal": spc["crc_vs_normal"],
    "replication": spc["replication"],
    "stage_trend": {"n_tests": spc["stage_trend"]["n_tests"],
                    "n_fdr05": spc["stage_trend"]["n_fdr05"]},
    "top_up": cxb.sort_values("delta", ascending=False).head(8)[
        ["gene", "twas_direction", "mean_CRC", "mean_normal_epi", "delta", "p"]
    ].to_dict("records"),
}


# ---------------------------------------------------------------- summary statistics
sstat = json.load(open("results/summary_statistics_overview.json"))
SW = pd.read_csv("results/risk_gene_summary_statistics.csv")
dprox = pd.read_csv("results/risk_gene_driver_proximity.csv")
SW2 = SW.merge(dprox[["gene", "is_driver", "dist_to_nearest_driver"]], on="gene", how="left")
crc_rel = SW2[(SW2["dist_to_nearest_driver"] <= 1) | (SW2["is_driver"] == True)]
D["summary_stats"] = {
    **sstat,
    "n_crc_relevant": int(len(crc_rel)),
    "n_crc_relevant_sig": int(((crc_rel["lead_FDR"] < 0.05)).sum()),
    "figure18_genes": ["SMARCD1", "CSNK2B", "PTGES3", "GTF2I", "UBE2V1", "UBE2B",
                       "PSMD3", "LAMC1", "SMAD9", "BMP2", "WNT4", "HLA-G"],
    "figure18_rows": SW2[SW2["gene"].isin(["SMARCD1", "CSNK2B", "PTGES3", "GTF2I",
                                           "UBE2V1", "UBE2B", "PSMD3", "LAMC1",
                                           "SMAD9", "BMP2", "WNT4", "HLA-G"])]
        [["gene", "class", "lead_comparison", "lead_log2FC", "lead_P_nominal",
          "lead_FDR", "twas_cell_types", "cosmx_delta_CRC_minus_normal_epi",
          "cosmx_P_paired"]].sort_values("lead_FDR").to_dict("records"),
    "top_by_lead_fdr": SW.head(16)[["gene", "twas_direction", "twas_cell_types",
                                    "lead_comparison", "lead_log2FC", "lead_P_nominal",
                                    "lead_FDR", "validated", "class"]].to_dict("records"),
}


# ---------------------------------------------------------------- functional screens
fq = json.load(open("results/func_qc.json"))
fs = json.load(open("results/func_summary.json"))
frg = pd.read_csv("results/func_risk_gene_results.csv")
fsi = pd.read_csv("results/func_sirna_gene_stats.csv")
fcr = pd.read_csv("results/func_crispr_gene_stats.csv")
D["screens"] = {
    "qc": fq["per_line"], "excluded": fq["excluded_cell_line"],
    "retained": fs["retained"],
    "n_genes_sirna": fq["n_genes_sirna"], "n_genes_crispr": fq["n_genes_crispr"],
    "n_detected_sirna": int(fsi["detected"].sum()),
    "n_strong_sirna": int(fsi["strong"].sum()),
    "n_detected_crispr": int(fcr["detected"].sum()),
    "n_strong_crispr": int(fcr["strong"].sum()),
    "sirna": fs["siRNA"], "crispr": fs["CRISPR"],
    "calibration": fs["calibration"], "concordance": fs["concordance"],
    "sensitivity": fs["sensitivity_three_lines"],
    "risk_hits": frg[frg["detected"]].sort_values("mean_z")[
        ["gene", "twas_class", "is_validated", "mean_z", "mean_pct_NT", "fdr", "strong",
         "crispr_mean_z", "crispr_hit"]].to_dict("records"),
    "n_risk_assayed": int(len(frg)),
    "zprime_sirna": {r["cell_line"]: r["z_prime"] for r in fq["per_line"]
                     if r["assay"] == "siRNA"},
}
_b = frg[frg["detected"]].dropna(subset=["crispr_mean_z"]).copy()
_b["agree"] = np.sign(_b["mean_z"]) == np.sign(_b["crispr_mean_z"])
_mm = pd.read_csv("results/func_sirna_vs_crispr.csv")
_contra = _mm[(_mm.detected_si | _mm.detected_cr) & (~_mm.agree_direction) & _mm.is_risk_si]
D["screens"]["corroborated"] = _b[(_b.agree) & (_b.mean_z < 0)][
    ["gene", "twas_class", "mean_pct_NT", "mean_z", "crispr_mean_z"]].to_dict("records")
D["screens"]["contradicted"] = _contra[
    ["gene", "mean_z_si", "fdr_si", "detected_si", "mean_z_cr", "fdr_cr", "detected_cr"]
].to_dict("records")
D["screens"]["n_corroborated"] = int(((_b.agree) & (_b.mean_z < 0)).sum())
_cl = json.load(open("results/crispr_by_cellline.json"))
_w = pd.read_csv("results/crispr_risk36_by_cellline.csv")
D["screens"]["by_cellline"] = {
    **{k: v for k, v in _cl.items() if k != "table"},
    "n_by_lines_nominal": {str(k): v for k, v in _cl["n_by_lines_nominal"].items()},
    "n_by_lines": {str(k): v for k, v in _cl["n_by_lines"].items()},
    "all_three": _w[_w["n_lines_supported"] == 3]["gene"].tolist(),
    "two_lines": _w[_w["n_lines_supported"] == 2]["gene"].tolist(),
    "top": _w.sort_values(["n_lines_supported", "mean_z_all"],
                          ascending=[False, True]).head(10)[
        ["gene", "class", "z_HCT116", "z_RKO", "z_SW480", "n_lines_supported",
         "supported_direction"]].to_dict("records"),
}
D["screens"]["n_contradicted"] = int(len(_contra))

_rq = json.load(open("results/crispr_replicate_qc.json"))
_rqc = pd.read_csv("results/crispr_replicate_qc.csv")
_cons = pd.DataFrame(_rq["consequence"])
_sw = _cons[_cons.cell_line == "SW480"].set_index("set")
# the pre-exclusion concordance, recomputed so the report can quote both readings
_pre = json.load(open("results/crispr_concordance_all_plates.json"))
D["screens"]["replicate_qc"] = {
    "excluded": _rq["excluded"], "n_excluded": _rq["n_excluded"],
    "by_plate": _rqc.to_dict("records"),
    "min_z_prime": float(_rqc["z_prime"].min()),
    "t_crit_df2": _rq["t_crit_df2"], "t_crit_df1": _rq["t_crit_df1"],
    "sw480_sig_all": int(_sw.loc["all 3", "sig_panel"]),
    "sw480_sig_qc": int(_sw.loc["after QC", "sig_panel"]),
    "sw480_effect_all": float(_sw.loc["all 3", "median_abs_effect"]),
    "sw480_effect_qc": float(_sw.loc["after QC", "median_abs_effect"]),
    "concordance_before": _pre["concordance_calls"],
    "concordance_before_genes": _pre["concordance_genes"],
    "concordance_after_qc": _pre.get("concordance_after_qc", _pre["concordance_calls"]),
}

def _within(path):
    c = pd.read_csv(path)
    w = c[c.same_cell_line]
    out = {"within_median": float(w.pearson_r.median()),
           "across_median": float(c[~c.same_cell_line].pearson_r.median()),
           "spearman_within_median": float(w.spearman_rho.median())}
    for ln in ["HCT116", "RKO", "SW480"]:
        out[ln] = sorted(w[w.plate_a.str.startswith(ln)].pearson_r.round(4).tolist())
    return out

D["screens"]["plate_cor"] = {"crispr": _within("results/crispr_plate_correlations.csv"),
                             "sirna": _within("results/sirna_plate_correlations.csv")}

_cq = pd.read_csv("results/cosmx_requested_genes.csv")
D["screens"]["cosmx_requested"] = {
    "requested": ["DDX39B", "CSNK2B", "CLIC1", "CTNNB1", "UBE2V1", "RPS5", "TRIM28", "LIMA1"],
    "drawn": _cq.loc[_cq.on_cosmx_panel, "gene"].tolist(),
    "off_panel": _cq.loc[~_cq.on_cosmx_panel, "gene"].tolist(),
    "table": _cq.to_dict("records"),
}

_hc = json.load(open("results/highconf_replication.json"))
_hcw = pd.read_csv("results/highconf_replication.csv")
D["screens"]["highconf"] = {
    **{k: v for k, v in _hc.items() if k != "odds_ratio"},
    "pct_same_screened": 100.0 * _hc["rho_screened"][2] / _hc["rho_screened"][3],
    "pct_same_other": 100.0 * _hc["rho_other"][2] / _hc["rho_other"][3],
    "pct_same_all": 100.0 * _hc["rho_all"][2] / _hc["rho_all"][3],
    "n_cosmx_sig": int((_hcw.p < 0.05).sum()),
    "n_same_sign": int(_hcw.replicates.sum()),
    "func_called_genes": sorted(_hcw.loc[_hcw.func_called, "gene"].tolist()),
    "discordant": _hcw.loc[~_hcw.replicates].reindex(
        _hcw.loc[~_hcw.replicates, "delta"].abs().sort_values(ascending=False).index)[
        ["gene", "geomx_log2FC", "delta", "p"]].to_dict("records"),
    "n_discordant_sig": int(((~_hcw.replicates) & (_hcw.p < 0.05)).sum()),
    "n_discordant_small": int(((~_hcw.replicates) & (_hcw.delta.abs() < 0.05)).sum()),
    "top": _hcw.reindex(_hcw.delta.abs().sort_values(ascending=False).index).head(10)[
        ["gene", "geomx_log2FC", "delta", "p", "geomx_validated", "func_called"]
    ].to_dict("records"),
    "table": _hcw[["gene", "geomx_log2FC", "delta", "p", "fdr", "twas_direction",
                   "geomx_validated", "replicates", "func_assayed", "func_called"]]
    .to_dict("records"),
}

_cs = json.load(open("results/cosmx_selected_genes.json"))
D["screens"]["cosmx_selected"] = {
    **{k: v for k, v in _cs.items() if k != "stats"},
    "stats": _cs["stats"],
    "genes": _cs["picked"],
}

D["screens"]["risk_increasing"] = json.load(open("results/risk_increasing_summary.json"))
D["screens"]["ri_crispr"] = json.load(open("results/risk_increasing_crispr.json"))

_sh = pd.read_csv("results/shared_genes_both_assays.csv")
from scipy.stats import spearmanr as _sp
_r, _p = _sp(_sh["cr_mean"], _sh["si_mean"])
_both = _sh[_sh.called_in == "both modalities"]
D["screens"]["shared"] = {
    "n_shared": int(len(_sh)),
    "only_crispr": sorted(set(pd.read_csv("results/crispr_risk36_by_cellline.csv")["gene"])
                          - set(_sh["gene"])),
    "counts": {k: int(v) for k, v in _sh["called_in"].value_counts().items()},
    "both_genes": _both.sort_values("si_mean")["gene"].tolist(),
    "both_same_sign": int((np.sign(_both.cr_mean) == np.sign(_both.si_mean)).sum()),
    "same_sign": int(_sh["same_sign"].sum()),
    "rho": float(_r), "p": float(_p),
    "table": _sh[["gene", "twas_direction", "cr_mean", "si_mean", "cr_verdict",
                  "si_verdict", "called_in"]].to_dict("records"),
}

_nt = json.load(open("results/nontwas_comparison.json"))
_ntw = pd.read_csv("results/nontwas_by_cellline.csv")
D["screens"]["nontwas"] = {
    **_nt,
    "n_by_lines_other": {str(k): v for k, v in _nt["n_by_lines_other"].items()},
    "n_by_lines_risk": {str(k): v for k, v in _nt["n_by_lines_risk"].items()},
    "top": _ntw.nsmallest(10, "mean_delta")[
        ["gene", "core", "delta_HCT116", "delta_RKO", "delta_SW480", "n_lines_sig"]
    ].to_dict("records"),
}

_sl = json.load(open("results/sirna_by_cellline.json"))
_sw = pd.read_csv("results/sirna_risk_by_cellline.csv")
D["screens"]["sirna_by_cellline"] = {
    **{k: v for k, v in _sl.items() if k != "table"},
    "n_by_lines": {crit: {str(k): v for k, v in d.items()}
                   for crit, d in _sl["n_by_lines"].items()},
    "all_three": _sw[_sw["n_lines_supported"] == 3]["gene"].tolist(),
    "two_lines": _sw[_sw["n_lines_supported"] == 2]["gene"].tolist(),
    "one_line": _sw[_sw["n_lines_supported"] == 1]["gene"].tolist(),
    "top": _sw.sort_values(["n_lines_supported", "mean_z_all"],
                           ascending=[False, True]).head(12)[
        ["gene", "twas_direction", "class", "z_HCT116", "z_RKO", "z_SW480",
         "n_lines_supported", "n_lines_nominal", "supported_direction"]].to_dict("records"),
}


# ---------------------------------------------------------------- pathway enrichment
_pj = json.load(open("results/pathway_celltype.json"))
_ps = pd.read_csv("results/pathway_celltype_significant.csv")
_pn = pd.read_csv("results/pathway_celltype_null.csv")
_pg = pd.read_csv("results/pathway_leading_genes.csv")
_pp = pd.read_csv("results/pathway_pooled_enrichment.csv")
_pt = _pn[_pn.tested].copy()
_surv = _pt[_pt.p_perm_matched_fdr < 0.05].sort_values("p_perm_matched_fdr")
CTL = {"ABS": "absorptive", "ASC": "absorptive stem-like", "CT": "crypt top",
       "EE": "enteroendocrine", "GOB": "goblet", "SSC": "secretory stem-like",
       "STM": "stem", "TAC": "transit-amplifying", "TUF": "tuft"}
D["pathways_ct"] = {
    **{k: v for k, v in _pj.items()
       if k not in ("null", "families", "pooled_top", "novel_top", "leading_genes",
                    "gene_core", "per_celltype")},
    "n_celltypes_with_hits": int(_ps.groupby(["tissue", "cell_type"]).ngroups),
    "n_core_genes": int(len({g for x in _ps.genes for g in x.split("|")})),
    "surviving": [{"label": f"{r.tissue.lower()} {CTL.get(r.cell_type, r.cell_type)}",
                   "tissue": r.tissue, "cell_type": r.cell_type, "n_mapped": int(r.n_mapped),
                   "observed": int(r.observed_fdr),
                   "matched_mean": float(r.matched_fdr_mean),
                   "p_matched": float(r.p_perm_matched_fdr),
                   "p_size": float(r.p_perm_fdr)} for r in _surv.itertuples(index=False)],
    "failed": [{"label": f"{r.tissue.lower()} {CTL.get(r.cell_type, r.cell_type)}",
                "observed": int(r.observed_fdr), "matched_mean": float(r.matched_fdr_mean),
                "p_matched": float(r.p_perm_matched_fdr)}
               for r in _pt[(_pt.observed_fdr >= 5) &
                            (_pt.p_perm_matched_fdr >= 0.05)].itertuples(index=False)],
    "families": _pj["families"],
    "families_missed": [f for f in ["WNT", "TGF-β / BMP", "RTK - MAPK", "PI3K - AKT", "p53",
                                    "Mismatch repair", "NOTCH", "Cell cycle", "Apoptosis",
                                    "Chromatin remodelling"]
                        if f not in set(_ps[_ps.established].family)],
    "novel_top": _pj["novel_top"],
    "gene_core": _pj["gene_core"],
    "leading": _pg.head(12).to_dict("records"),
    "pooled_top": _pp[_pp.fdr < 0.1][
        ["pathway", "source", "family", "n_pathway", "overlap", "p", "fdr"]].to_dict("records"),
    "n_blocks": int(_ps.groupby("pathway").agg(
        core=("genes", lambda s: "|".join(sorted({g for x in s for g in x.split("|")}))))
        .core.nunique()),
    "by_celltype": [{"label": f"{t.lower()} {CTL.get(c, c)}", "tissue": t, "cell_type": c,
                     "n": int(len(g)), "n_est": int(g.established.sum()),
                     "core": sorted({x for y in g.genes for x in y.split("|")})}
                    for (t, c), g in _ps.groupby(["tissue", "cell_type"])],
}

_mj = json.load(open("results/pathway_merged.json"))
_ms = pd.read_csv("results/pathway_merged_significant.csv")
_mn = pd.read_csv("results/pathway_merged_null.csv")
_mc = pd.read_csv("results/pathway_merged_vs_split.csv")
D["pathways_merged"] = {
    **{k: v for k, v in _mj.items()
       if k not in ("modules", "compare", "families", "top", "established_top",
                    "surviving")},
    "surviving": _mj["surviving"],
    "modules": _mj["modules"],
    "compare": _mj["compare"],
    "families": _mj["families"],
    "established": _mj["established_top"],
    "by_module": [{"label": l,
                   "n": int(_mn.loc[_mn.label == l, "n_mapped"].iat[0]),
                   "n_normal": int(_mn.loc[_mn.label == l, "n_normal"].iat[0]),
                   "n_adenoma": int(_mn.loc[_mn.label == l, "n_adenoma"].iat[0]),
                   "n_shared": int(_mn.loc[_mn.label == l, "n_shared"].iat[0]),
                   "sets": int(len(g)), "established": int(g.established.sum()),
                   "both": int((g.tissues_contributing == "normal+adenoma").sum()),
                   "p_matched": float(_mn.loc[_mn.label == l, "p_perm_matched_fdr"].iat[0]),
                   "core": sorted({x for y in g.genes for x in y.split("|")})}
                  for l, g in _ms.groupby("label")],
    "stem_top": _ms[_ms.label == "stem"].nsmallest(5, "p")[
        ["pathway", "overlap", "fdr_merged", "tissues_contributing", "genes"]
    ].to_dict("records"),
    "goblet_top": _ms[_ms.label == "goblet"].nsmallest(6, "p")[
        ["pathway", "family", "overlap", "fdr_merged", "tissues_contributing", "genes"]
    ].to_dict("records"),
}

_ej = json.load(open("results/all347_etiology.json"))
_em = pd.read_csv("results/canonical_membership_summary.csv")
_et = pd.read_csv("results/etiology_themes.csv")
_ea = pd.read_csv("results/etiology_gene_assignment.csv")
_emg = pd.read_csv("results/canonical_membership.csv")
D["etiology"] = {
    **{k: v for k, v in _ej.items() if k not in ("all_top", "up_top", "down_top",
                                                 "families", "themes", "proximity",
                                                 "nulls", "status")},
    "nulls": _ej["nulls"], "status": _ej["status"], "proximity": _ej["proximity"],
    "enriched": [r for r in _ej["all_top"] if float(r["fdr"]) < 0.1],
    "families": _em.to_dict("records"),
    "n_in_any_family": int(_emg.gene.nunique()),
    "n_multi_family": int((_emg.groupby("gene").family.nunique() >= 3).sum()),
    "multi_family_genes": sorted(_emg.groupby("gene").family.nunique()
                                 .pipe(lambda s: s[s >= 3]).index.tolist()),
    "themes": _et.to_dict("records"),
    "themes_sig": _et[_et.fdr < 0.1].to_dict("records"),
    "theme_zero": _et[_et.n_risk_genes == 0].to_dict("records"),
    "n_annotated": int(_ea.in_pathway_universe.sum()),
    "n_no_entrez": int(D["n_genes"] - _ej["nulls"][0]["n_genes"]),
}

_dj = json.load(open("results/driver_overlap.json"))
_dv = pd.read_csv("results/driver_overlap.csv")
_ord = {"established": 0, "probable": 1, "contextual": 2}
_dv = (_dv.assign(_k=_dv.confidence.map(_ord), _c=~_dv.crc_driver_list)
       .sort_values(["_c", "_k", "p_twas"]).drop(columns=["_k", "_c"]))
D["drivers"] = {
    **{k: v for k, v in _dj.items() if k != "table"},
    "table": _dv.to_dict("records"),
    "n_total": int(len(_dv)),
    "n_established": int((_dv.confidence == "established").sum()),
    "n_probable": int((_dv.confidence == "probable").sum()),
    "n_contextual": int((_dv.confidence == "contextual").sum()),
    "n_consistent": int((_dv.direction_vs_role == "consistent").sum()),
    "n_opposite": int((_dv.direction_vs_role == "opposite").sum()),
    "n_testable_direction": int(_dv.direction_vs_role.notna().sum()),
    "opposite_genes": _dv.loc[_dv.direction_vs_role == "opposite", "gene"].tolist(),
    "n_validated": int(_dv.spatially_validated.sum()),
    "validated_genes": _dv.loc[_dv.spatially_validated, "gene"].tolist(),
    "n_screened": int(_dv.screen_effect.sum()),
    "screened_genes": _dv.loc[_dv.screen_effect, "gene"].tolist(),
}


def clean(o):
    if isinstance(o, dict):
        return {k: clean(v) for k, v in o.items()}
    if isinstance(o, list):
        return [clean(v) for v in o]
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return None if (np.isnan(o) or np.isinf(o)) else float(o)
    if isinstance(o, (np.bool_,)):
        return bool(o)
    if isinstance(o, float) and (np.isnan(o) or np.isinf(o)):
        return None
    return o


json.dump(clean(D), open("results/report_data.json", "w"), indent=1)
print("report data written;", len(D), "top-level keys")
print(json.dumps(clean(D["erbb2"]), indent=1)[:600])
