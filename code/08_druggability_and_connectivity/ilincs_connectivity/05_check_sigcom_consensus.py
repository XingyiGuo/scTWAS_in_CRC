#!/usr/bin/env python
"""Consensus reversal query against the full LINCS collection via SigCom LINCS.

The local scan (scripts 01–04) covers the 16,767 large-intestine iLINCS
signatures. SigCom LINCS indexes the whole L1000 collection and scores a tagged
up/down query against it server-side with its own two-sided rank statistic,
returning the strongest mimickers and reversers with p, FDR and Bonferroni.

Two libraries are queried per cell-type model:
  l1000_cp       693,215 individual chemical perturbagen signatures
  l1000_mean_cp  consensus signature per compound x cell line, i.e. the
                 aggregate view that averages over dose and replicate

A reverser is a signature whose expression change opposes the query: our
up-tagged genes (most positive TWAS Z) go down and vice versa - the
therapeutically interesting direction.

Outputs (under <ROOT>/iLINCS_reversal):
  table/sigcom_reversers_top.csv        top reversers per model and library, annotated
  table/sigcom_model_summary.csv        per model: library size, best FDR, counts
  table/sigcom_recurrent_reversers.csv  recurrent perturbagens across models
"""
import json
import re
import time
import urllib.error
import urllib.parse
import urllib.request

import pandas as pd

OUT = r"W:\VUMC_projects\scTWAS\iLINCS_reversal"
META_API = "https://maayanlab.cloud/sigcom-lincs/metadata-api"
ENRICH = "https://data-api.sigcom-lincs.k8s.maayanlab.cloud/v1/api/v1/enrich/ranktwosided"
LIBS = ["l1000_cp", "l1000_mean_cp"]
TOP = 100


def post(url, payload, timeout=300, tries=3):
    for a in range(tries):
        try:
            req = urllib.request.Request(url, data=json.dumps(payload).encode(),
                                         headers={"Content-Type": "application/json"})
            with urllib.request.urlopen(req, timeout=timeout) as f:
                return json.load(f)
        except (urllib.error.HTTPError, urllib.error.URLError, TimeoutError):
            if a == tries - 1:
                raise
            time.sleep(3 * (a + 1))


def entity_ids(symbols):
    out = {}
    for i in range(0, len(symbols), 200):
        r = post(f"{META_API}/entities/find",
                 {"filter": {"where": {"meta.symbol": {"inq": symbols[i:i + 200]}}, "limit": 500}})
        for e in r:
            out[e["meta"]["symbol"]] = e["id"]
    return out


def signature_meta(uuids):
    """uuid -> perturbagen, cell line, dose, time."""
    rows = []
    for i in range(0, len(uuids), 100):
        r = post(f"{META_API}/signatures/find",
                 {"filter": {"where": {"id": {"inq": uuids[i:i + 100]}}, "limit": 200}})
        for s in r:
            m = s.get("meta", {})
            rows.append({"uuid": s["id"],
                         "perturbagen": m.get("pert_name") or m.get("pert_id"),
                         "pert_id": m.get("pert_id"),
                         "cell_line": m.get("cell_line") or m.get("cell_id"),
                         "dose": m.get("pert_dose"), "timepoint": m.get("pert_time"),
                         "moa": m.get("moa"), "target": m.get("target")})
    return pd.DataFrame(rows)


def main():
    qdf = pd.read_csv(rf"{OUT}\table\ilincs_query_genes.csv")
    ent = entity_ids(sorted(qdf.gene_name.unique()))
    print("query genes mapped to SigCom entities:", len(ent))

    rows, summary = [], []
    for model, q in qdf.groupby("model"):
        up = [ent[g] for g in q.loc[q.side.eq("positive"), "gene_name"] if g in ent]
        dn = [ent[g] for g in q.loc[q.side.eq("negative"), "gene_name"] if g in ent]
        for lib in LIBS:
            r = post(ENRICH, {"up_entities": up, "down_entities": dn,
                              "limit": TOP, "database": lib})
            res = pd.DataFrame(r["results"])
            rev = res[res.type.eq("reversers")].copy()
            rev["model"], rev["library"] = model, lib
            rows.append(rev)
            best = rev.loc[rev["z-sum"].idxmin()] if len(rev) else None
            summary.append({
                "model": model, "library": lib, "library_signatures": r["maxRank"],
                "n_reversers": r["reversers"], "n_mimickers": r["mimickers"],
                "best_z_sum": None if best is None else best["z-sum"],
                "best_fdr_up": None if best is None else best["fdr-up"],
                "best_fdr_down": None if best is None else best["fdr-down"],
                "best_logp_fisher": None if best is None else best["logp-fisher"],
                "n_reversers_fdr05": int(((rev["fdr-up"] < 0.05) & (rev["fdr-down"] < 0.05)).sum()),
                "n_reversers_bonferroni05": int(((rev["p-up-bonferroni"] < 0.05) &
                                                 (rev["p-down-bonferroni"] < 0.05)).sum()),
                "query_seconds": r["queryTimeSec"]})
            print(f"{model:30s} {lib:14s} sigs {r['maxRank']:7d} | best z-sum "
                  f"{summary[-1]['best_z_sum']:+.2f} fdr {min(summary[-1]['best_fdr_up'], summary[-1]['best_fdr_down']):.3g}"
                  f" | FDR<0.05 {summary[-1]['n_reversers_fdr05']}")

    allrev = pd.concat(rows, ignore_index=True)
    meta = signature_meta(sorted(allrev.uuid.unique()))
    allrev = allrev.merge(meta, on="uuid", how="left")
    summ = pd.DataFrame(summary)
    allrev.to_csv(rf"{OUT}\table\sigcom_reversers_top.csv", index=False)
    summ.to_csv(rf"{OUT}\table\sigcom_model_summary.csv", index=False)

    rec = (allrev.groupby(["library", "perturbagen"])
           .agg(n_models=("model", "nunique"), n_signatures=("uuid", "nunique"),
                best_z=("z-sum", "min"), best_fdr=("fdr-up", "min"))
           .reset_index().sort_values(["library", "n_models", "best_z"], ascending=[True, False, True]))
    rec.to_csv(rf"{OUT}\table\sigcom_recurrent_reversers.csv", index=False)
    print("\n", summ.groupby("library")[["library_signatures", "n_reversers_fdr05",
                                         "n_reversers_bonferroni05"]].agg(["min", "max"]).to_string())
    print("\nmost recurrent reversers per library:")
    for lib in LIBS:
        print(rec[rec.library.eq(lib)].head(8).to_string(index=False))
    return allrev, summ, rec


if __name__ == "__main__":
    allrev, summ, rec = main()
