## 08_summarize_gene_level_support.R --------------------------------
## Collapse the 546 scTWAS gene x cell-type pairs to the 347 unique genes, fix
## one direction per gene, and check that direction against every available
## differential-expression resource.
##
## Why this exists
##   The support tables for scRNA-seq and spatial expression are
##   per-pair: a gene nominated in five cell-type models occupies five rows and
##   is counted five times. Every external DE resource used here is gene level,
##   so the per-pair layout multiplies one measurement across rows and makes the
##   gene appear to be supported repeatedly. This script produces the gene-level
##   view instead: one row per gene, one direction, one support call per source.
##
## Direction rule (as specified for this analysis)
##   A gene nominated in several cell-type models can carry Z scores of opposite
##   sign in different models. Where all models agree the common sign is used.
##   Where they disagree the representative direction is taken from the model
##   with the smallest scTWAS P value. 10 of the 347 genes are discordant and
##   are flagged, so any conclusion resting on them can be inspected.
##
## Input
##   data/scTWAS_ADSER_NL_collected_results.csv                546 pairs, 16 models
##   data/risk_gene_all_comparisons.csv                        GeoMx contrasts
##   data/TCGA_CPTAC_support.csv                               header on row 2
##   <Visium repo>/Visium/results/tables/de_senescence_low_grade_vs_normal.csv   M1
##   <Visium repo>/Visium/results/tables/de_stemness_high_vs_low.csv             M3
##
## Output
##   table/scTWAS_347genes_model_Zscore_matrix.csv
##        347 x 16 Z-score matrix, one column per Dataset_CellType model
##   table/scTWAS_347genes_direction_summary.csv
##        per gene: models, up/down counts, concordance, representative Z/P/model
##   table/scTWAS_347genes_multisource_DE_support.csv
##        the summary plus one support block per DE resource
##   table/scTWAS_347genes_support_counts.csv
##        label counts per resource, with denominators
##
## Support rule, applied identically to all seven resources so the columns are
## comparable: FDR < 0.05 AND |log2FC| > 0.5 AND sign(log2FC) == sign(rep Z).
##   Up             supported, gene higher in the more advanced state
##   Down           supported, gene lower in the more advanced state
##   DE_discordant  passes FDR and the effect floor but opposes the scTWAS Z
##   NS             tested, does not pass FDR or the effect floor
##   not_tested     gene absent from that resource
## The published TCGA_DE_support / CPTAC_DE_support flags are carried through
## verbatim as well; they use an FDR cut only, with no effect-size floor, so
## they are looser than the columns above by construction.
##
## One-sidedness caveat, Visium only
##   The deposited Visium tables contain only genes enriched in each contrast
##   (adjusted P < 0.05, log fold-change > 0.2, detection difference > 10%), so
##   every log2FC there is positive. A gene whose representative Z is negative
##   can never be direction-concordant against M1 or M3 under any threshold.
##   The Visium_assessable column marks which genes can be tested at all.
## -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
})

## paths -----------------------------------------------------------------------
PROJ <- "W:/VUMC_projects/scTWAS/VUMC_HTAN_scTWAS_genes_DE_support"
if (!dir.exists(file.path(PROJ, "data"))) PROJ <- "."

DIR_DATA <- file.path(PROJ, "data")
DIR_TAB  <- file.path(PROJ, "table")
dir.create(DIR_TAB, recursive = TRUE, showWarnings = FALSE)

TWAS_CSV <- file.path(DIR_DATA, "scTWAS_ADSER_NL_collected_results.csv")
SPAT_CSV <- file.path(DIR_DATA, "risk_gene_all_comparisons.csv")
if (!file.exists(SPAT_CSV)) {
  SPAT_CSV <- file.path(DIR_DATA, "spatialDE_CosMx_PMID41986711_risk_gene_all_comparisons.csv")
}
BULK_CSV <- file.path(DIR_DATA, "TCGA_CPTAC_support.csv")

## The Visium repository has been checked out under either project root.
.vis <- c("W:/VUMC_projects/scTWAS/CRC_SpatilTrans_PMID42620287/Spatial-multi-omics/Visium/results/tables",
          "L:/SpatialExpr/CRC_SpatilTrans_PMID42620287/Spatial-multi-omics/Visium/results/tables")
VIS_DIR <- .vis[dir.exists(.vis)][1]
if (is.na(VIS_DIR)) stop("Visium results/tables not found. Checked:\n  ",
                         paste(.vis, collapse = "\n  "))

## settings --------------------------------------------------------------------
Q_CUT      <- 0.05
LFC_CUT    <- 0.5
FDR_COL    <- "FDR_risk_restricted"
PANEL_MAIN <- "exp2"          # WTA, carries all 271 spatial genes
CMP_MAIN   <- "Epithelium"    # all nine scTWAS cell types are epithelial

CT_ORDER <- c("STM", "TAC", "ABS", "CT", "GOB", "EE", "TUF", "ASC", "SSC")

## =============================================================================
## Functions
## =============================================================================

## Support label, identical semantics to 08_sctwas_spatialDE_CosMx_support.R.
label_support <- function(log2FC, fdr, z, q_cut = Q_CUT, lfc_cut = LFC_CUT) {
  out  <- rep("NS", length(log2FC))
  sig  <- !is.na(fdr) & fdr < q_cut & !is.na(log2FC) & abs(log2FC) > lfc_cut
  conc <- sign(log2FC) == sign(z)
  out[sig &  conc & log2FC > 0] <- "Up"
  out[sig &  conc & log2FC < 0] <- "Down"
  out[sig & !conc]              <- "DE_discordant"
  out[is.na(log2FC)]            <- "not_tested"
  out
}


## 347 x 16 matrix of scTWAS Z scores. A cell holds the Z score when the gene
## reached significance in that Dataset_CellType model and is empty otherwise.
build_model_z_matrix <- function(twas, ct_order = CT_ORDER) {
  d <- copy(twas)
  d[, model := paste(Dataset, CellType, sep = "_")]
  stopifnot(!anyDuplicated(d, by = c("gene", "model")))

  w <- dcast(d, gene + gene_name ~ model, value.var = "Z_score")

  ## column order: NL models first, then ADSER, each in crypt-axis order
  lev <- unlist(lapply(c("NL", "ADSER"),
                       function(ds) paste(ds, ct_order, sep = "_")))
  keep <- intersect(lev, names(w))
  setcolorder(w, c("gene", "gene_name", keep))
  attr(w, "model_cols") <- keep
  w[]
}


## One row per gene: how many models nominated it, how the signs split, whether
## they agree, and the representative Z / P / model under the smallest-P rule.
summarize_gene_direction <- function(twas) {
  d <- copy(twas)
  d[, model := paste(Dataset, CellType, sep = "_")]
  setorder(d, gene, P_value)          # smallest scTWAS P first within each gene

  s <- d[, {
    up <- sum(Z_score > 0); dn <- sum(Z_score < 0)
    list(
      n_models      = .N,
      models        = paste(model, collapse = ","),
      cell_types    = paste(sort(unique(CellType)), collapse = ","),
      datasets      = paste(sort(unique(Dataset)),  collapse = ","),
      n_up          = up,
      n_down        = dn,
      ## Single_model, Concordant when every model agrees in sign, otherwise
      ## Discordant; the representative below resolves the Discordant genes.
      direction_status = if (.N == 1L) "Single_model"
                         else if (up == 0L || dn == 0L) "Concordant"
                         else "Discordant",
      rep_model     = model[1],       # smallest P after the setorder above
      rep_CellType  = CellType[1],
      rep_Dataset   = Dataset[1],
      rep_Z         = Z_score[1],
      rep_P         = P_value[1],
      rep_Sig_resource = Sig_resource[1],
      Z_min         = min(Z_score),
      Z_max         = max(Z_score)
    )
  }, by = .(gene, gene_name)]

  s[, rep_direction := ifelse(rep_Z > 0, "risk-increasing", "risk-decreasing")]
  s[, Visium_assessable := rep_Z > 0]   # see the one-sidedness note in the header
  setorder(s, -n_models, rep_P)
  s[]
}


## Attach one support block (log2FC, padj, DEsupport) for a gene-level DE table.
## `de` must hold one row per gene with columns gene_symbol, log2FC, fdr.
add_support_block <- function(summary_dt, de, prefix, z_col = "rep_Z") {
  stopifnot(all(c("gene_symbol", "log2FC", "fdr") %in% names(de)))
  stopifnot(!anyDuplicated(de, by = "gene_symbol"))
  m <- de[.(toupper(summary_dt$gene_name)), on = .(gene_symbol)]
  summary_dt[[paste0(prefix, "_log2FC")]] <- m$log2FC
  summary_dt[[paste0(prefix, "_padj")]]   <- m$fdr
  summary_dt[[paste0(prefix, "_DEsupport")]] <-
    label_support(m$log2FC, m$fdr, summary_dt[[z_col]])
  summary_dt
}

## =============================================================================
## 1. scTWAS input, Z matrix, direction summary
## =============================================================================

twas <- fread(TWAS_CSV)
stopifnot(all(c("Dataset", "CellType", "gene", "gene_name",
                "Sig_resource", "Z_score", "P_value") %in% names(twas)))
stopifnot(nrow(twas) == 546L, uniqueN(twas$gene) == 347L,
          uniqueN(twas$gene_name) == 347L,
          twas[, uniqueN(gene_name), by = gene][, max(V1)] == 1L)

zmat <- build_model_z_matrix(twas)
MODEL_COLS <- attr(zmat, "model_cols")
stopifnot(nrow(zmat) == 347L, length(MODEL_COLS) == 16L)
## every non-missing cell must correspond to a real scTWAS pair
stopifnot(sum(!is.na(as.matrix(zmat[, ..MODEL_COLS]))) == nrow(twas))

gsum <- summarize_gene_direction(twas)
stopifnot(nrow(gsum) == 347L, sum(gsum$n_models) == nrow(twas))

## write the matrix with the direction summary appended on the right
zout <- merge(zmat, gsum[, .(gene, n_models, n_up, n_down, direction_status,
                             rep_model, rep_Z, rep_P, rep_direction)],
              by = "gene", sort = FALSE)
setcolorder(zout, c("gene", "gene_name", MODEL_COLS))
fwrite(zout, file.path(DIR_TAB, "scTWAS_347genes_model_Zscore_matrix.csv"), na = "")
fwrite(gsum, file.path(DIR_TAB, "scTWAS_347genes_direction_summary.csv"), na = "")

## =============================================================================
## 2. DE resources, each reduced to one row per gene
## =============================================================================

sup <- copy(gsum)

## --- GeoMx: NOR->CORE, ADE->CORE, NOR->ADE (epithelium, WTA) ----------------
## Keep the historical CosMx_ column prefix for existing plot scripts.
spat <- fread(SPAT_CSV)
setnames(spat, FDR_COL, "fdr")
spat[, `:=`(log2FC = as.numeric(log2FC), fdr = as.numeric(fdr))]
## guard: the source file's own flags must still match their definitions
stopifnot(identical(spat$sig, spat$fdr < Q_CUT))

for (tr in c("NOR->CORE", "ADE->CORE", "NOR->ADE")) {
  sub <- spat[compartment_name == CMP_MAIN & experiment == PANEL_MAIN &
                transition == tr]
  if (!nrow(sub)) stop("no GeoMx rows for ", tr)
  de <- sub[, .(gene_symbol = toupper(gene), log2FC, fdr)]
  setkey(de, gene_symbol)
  sup <- add_support_block(sup, de,
                           paste0("CosMx_", gsub("->", "_", tr, fixed = TRUE)))
}

## CTA measurements are displayed alongside WTA in Table 20.
cta <- spat[compartment_name == CMP_MAIN & experiment == "exp1" &
              transition == "NOR->CORE",
            .(gene_symbol = toupper(gene), log2FC, fdr)]
stopifnot(nrow(cta) > 0L, !anyDuplicated(cta, by = "gene_symbol"))
cta <- cta[.(toupper(sup$gene_name)), on = .(gene_symbol)]
sup[, `:=`(GeoMx_CTA_NOR_CORE_log2FC = cta$log2FC,
           GeoMx_CTA_NOR_CORE_padj = cta$fdr)]
cta_support <- label_support(cta$log2FC, cta$fdr, sup$rep_Z)
## The combined display column stays separate from the WTA resource roll-up.
sup[, GeoMx_NOR_CORE_DEsupport_either := fcase(
  CosMx_NOR_CORE_DEsupport %chin% c("Up", "Down"), CosMx_NOR_CORE_DEsupport,
  cta_support %chin% c("Up", "Down"), cta_support,
  CosMx_NOR_CORE_DEsupport == "DE_discordant" | cta_support == "DE_discordant", "DE_discordant",
  CosMx_NOR_CORE_DEsupport == "not_tested" & cta_support == "not_tested", "not_tested",
  default = "NS"
)]

## --- Visium modules M1 and M3 ----------------------------------------------
read_visium <- function(f) {
  d <- fread(file.path(VIS_DIR, f))
  stopifnot(all(c("gene", "avg_log2FC", "p_val_adj") %in% names(d)))
  ## these tables are one-sided by construction; assert it so the caveat in the
  ## header cannot silently stop being true
  stopifnot(all(d$avg_log2FC > 0))
  de <- d[, .(gene_symbol = toupper(gene), log2FC = avg_log2FC,
              fdr = p_val_adj)]
  de <- de[order(fdr)][!duplicated(gene_symbol)]
  setkey(de, gene_symbol)
  de
}
sup <- add_support_block(sup, read_visium("de_senescence_low_grade_vs_normal.csv"),
                         "Visium_M1")
sup <- add_support_block(sup, read_visium("de_stemness_high_vs_low.csv"),
                         "Visium_M3")

## --- TCGA and CPTAC bulk ----------------------------------------------------
## Header sits on the second line of the file. log2FC/FDR are gene level and
## constant across the pairs of a gene, which is asserted below.
bulk <- fread(BULK_CSV, skip = 1L)
for (cc in c("Z_score", "TCGA_log2FC", "TCGA_FDR", "CPTAC_log2FC", "CPTAC_FDR"))
  bulk[[cc]] <- as.numeric(bulk[[cc]])
stopifnot(nrow(bulk) == 546L)

for (coh in c("TCGA", "CPTAC")) {
  lf <- paste0(coh, "_log2FC"); fd <- paste0(coh, "_FDR")
  chk <- bulk[, .(n = uniqueN(na.omit(get(lf))) + uniqueN(na.omit(get(fd)))),
              by = gene_name]
  stopifnot(chk[, max(n)] <= 2L)          # at most one distinct value each
  de <- unique(bulk[, .(gene_symbol = toupper(gene_name),
                        log2FC = get(lf), fdr = get(fd))],
               by = "gene_symbol")
  setkey(de, gene_symbol)
  sup <- add_support_block(sup, de, coh)
  ## carry the published per-pair flag through, taken from the representative model
  pub <- bulk[, .(Dataset, CellType, gene_name,
                  flag = get(paste0(coh, "_DE_support")))]
  pub[, model := paste(Dataset, CellType, sep = "_")]
  m <- pub[.(sup$gene_name, sup$rep_model), on = .(gene_name, model)]
  sup[[paste0(coh, "_DE_support_published")]] <- m$flag
}

## The published flags use an FDR cut with no effect-size floor. Where a gene
## clears the floor the two definitions must agree; that is the guard.
for (coh in c("TCGA", "CPTAC")) {
  lbl <- sup[[paste0(coh, "_DEsupport")]]
  pub <- sup[[paste0(coh, "_DE_support_published")]]
  ok  <- lbl %chin% c("Up", "Down")
  stopifnot(all(pub[ok] == "Support"))
  ok2 <- lbl == "DE_discordant"
  stopifnot(all(pub[ok2] == "Opposite_direction"))
}

## =============================================================================
## 3. Roll-up across resources
## =============================================================================

SOURCES <- c("CosMx_NOR_CORE", "CosMx_ADE_CORE", "CosMx_NOR_ADE",
             "Visium_M1", "Visium_M3", "TCGA", "CPTAC")
sup_cols <- paste0(SOURCES, "_DEsupport")

lbl_mat <- as.matrix(sup[, ..sup_cols])
## `%in%` on a matrix returns a plain vector, so the mask has to be rebuilt
## with the original dimensions before rowSums.
is_sup <- matrix(lbl_mat %chin% c("Up", "Down"),
                 nrow = nrow(lbl_mat), dimnames = dimnames(lbl_mat))
sup[, n_sources_tested     := rowSums(lbl_mat != "not_tested")]
sup[, n_sources_supporting := rowSums(is_sup)]
sup[, n_sources_discordant := rowSums(lbl_mat == "DE_discordant")]
sup[, supporting_sources := apply(is_sup, 1, function(r)
  paste(SOURCES[r], collapse = ","))]

## Table 20 counts either GeoMx cohort as one source, alongside TCGA and CPTAC.
table20_sources <- c("GeoMx_NOR_CORE", "TCGA", "CPTAC")
table20_labels <- as.matrix(sup[, .(GeoMx_NOR_CORE_DEsupport_either,
                                    TCGA_DEsupport, CPTAC_DEsupport)])
sup[, bulk_DE_supports := apply(table20_labels, 1L, function(labels)
  paste(table20_sources[labels %chin% c("Up", "Down")], collapse = ","))]

ord <- c("gene", "gene_name", "n_models", "models", "cell_types", "datasets",
         "n_up", "n_down", "direction_status",
         "rep_model", "rep_Dataset", "rep_CellType", "rep_Z", "rep_P",
         "rep_Sig_resource", "rep_direction", "Z_min", "Z_max",
         "Visium_assessable",
         unlist(lapply(SOURCES, function(s)
           paste0(s, c("_log2FC", "_padj", "_DEsupport")))),
         "TCGA_DE_support_published", "CPTAC_DE_support_published",
         "n_sources_tested", "n_sources_supporting", "n_sources_discordant",
         "supporting_sources", "bulk_DE_supports")
ord <- append(ord, c("GeoMx_CTA_NOR_CORE_log2FC", "GeoMx_CTA_NOR_CORE_padj",
                     "GeoMx_NOR_CORE_DEsupport_either"),
              after = match("CosMx_NOR_CORE_DEsupport", ord))
sup <- sup[, ..ord]
setorder(sup, -n_sources_supporting, -n_models, rep_P)
fwrite(sup, file.path(DIR_TAB, "scTWAS_347genes_multisource_DE_support.csv"),
       na = "")

## label counts per resource --------------------------------------------------
lv <- c("Up", "Down", "DE_discordant", "NS", "not_tested")
cnt <- rbindlist(lapply(SOURCES, function(s) {
  v <- factor(sup[[paste0(s, "_DEsupport")]], levels = lv)
  t <- as.list(table(v))
  c(list(resource = s, n_genes = nrow(sup)), t,
    list(supported = t$Up + t$Down,
         pct_of_tested = round(100 * (t$Up + t$Down) /
                                 max(1L, nrow(sup) - t$not_tested), 1)))
}))
fwrite(cnt, file.path(DIR_TAB, "scTWAS_347genes_support_counts.csv"), na = "")

## =============================================================================
## 4. Console summary
## =============================================================================

cat("\n== direction across the 16 cell-type models (347 genes) ==\n")
print(gsum[, .N, by = direction_status][order(-N)])
cat("\nmodels per gene:\n")
print(gsum[, .N, by = n_models][order(n_models)])
cat("\ndiscordant-direction genes (sign differs between models):\n")
print(gsum[direction_status == "Discordant",
           .(gene_name, n_models, n_up, n_down, rep_model,
             rep_Z = round(rep_Z, 2), rep_P, Z_min = round(Z_min, 2),
             Z_max = round(Z_max, 2))])

cat("\n== DE support per resource (347 genes) ==\n")
print(cnt)

cat("\n== genes by number of supporting resources ==\n")
print(sup[, .N, by = n_sources_supporting][order(-n_sources_supporting)])
cat("\ntop genes by converging support:\n")
print(head(sup[n_sources_supporting > 0,
               .(gene_name, n_models, direction_status, rep_model,
                 rep_Z = round(rep_Z, 2), n_sources_supporting,
                 supporting_sources)], 20))

cat("\nwrote 4 tables to ", DIR_TAB, "\n", sep = "")
