#install.packages("remotes")
#install.packages("data.table")
#remotes::install_github("wjawaid/enrichR")

# module load StdEnv/2023
# module load r/4.5.0

library(data.table)
library("remotes")
library("readxl")
library("enrichR")
library("dplyr")
listEnrichrSites()
setEnrichrSite("Enrichr")

format_output <- function(enriched, ouptut_prefix, write_file = TRUE){
  header <- c("Term", "Overlap", "P.value", "Adjusted.P.value",
              "Old.P.value", "Old.Adjusted.P.value", "Odds.Ratio",
              "Combined.Score", "Genes")

  filtered_list <- list()

  for (sheet_name in names(enriched)) {
    df <- enriched[[sheet_name]]

    if ("Adjusted.P.value" %in% colnames(df)) {
      df_sig <- df[df$Adjusted.P.value < 0.1, , drop = FALSE]

      if (nrow(df_sig) > 0) {
        df_sig$Source <- sheet_name
        filtered_list[[sheet_name]] <- df_sig
      }
    }
  }

  if (length(filtered_list) > 0) {
    combined_df <- do.call(rbind, filtered_list)
    combined_df <- combined_df[, c(header, "Source"), drop = FALSE]

    if (write_file) {
      fwrite(
        combined_df,
        file = paste0(
          "/data/l2_bioinfo1/liq17/scTWAS/Functional_res/EnrichR/",
          ouptut_prefix, "_Pathway_enriched_genes.tsv"
        ),
        sep = "\t"
      )
    }

    return(combined_df)
  }

  message("No significant results (Adjusted.P.value < 0.1) found in any enrichment set.")
  return(NULL)
}

dbs_lib <- c("Reactome_Pathways_2024", "KEGG_2021_Human", "MSigDB_Hallmark_2020")
#genes_df <- read_excel("D:/PostDoc/Projects/CRC_scTWAS/ms_figures_tables/scTWAS_Supplementary_Tables.xlsx", sheet="S1", skip=1)
wk="/data/l2_bioinfo1/liq17/scTWAS/TWAS_HTAN_GSE266546_WellcomeSanger_res/SPrediXcan_peakSNP_association/"

test="ADSER"
if(test=="ADSER"){	
	genes_df <- data.frame(fread(paste0(wk, "VUMC_HTAN_ADSER_EUR_Multi_sig/VUMC_HTAN_ADSER_9CTs_peakSNP.EUR_Multi.sig.csv")))
	model_cts <- c("ABS","CT","EE","GOB","STM","TAC","ASC","SSC")
}else if(test=="NL"){
	genes_df <-data.frame(fread(paste0(wk, "Cauchy_Meta_NL_WellcomeSanger_NL_EUR_Multi/HTAN_NL_WellcomeSanger_NL_7celltypes_peakSNP.EUR_Multi.CauchyMeta.sig.csv")))
	model_cts <- c("ABS","CT","EE","GOB","STM","TAC")
}

gene_cols <- grep("^gene_name_", names(genes_df), value = TRUE)
tmp <- lapply(genes_df[gene_cols],function(x) {
    x <- trimws(x)
    x[x == ""] <- NA
    x
  }
)

genes_df$gene_name <- do.call(dplyr::coalesce, tmp)
input <- unique(na.omit(genes_df$gene_name))
print(length(input))
enriched <- enrichr(input, dbs_lib)
format_output(enriched, test)

genes_df %>% group_by(CellType) %>% summarise(gene_name = n())

# Run Enrichr separately for each cell type and collect all significant pathways
celltype_enrichr_list <- list()

for(model in unique(genes_df$CellType)){
  model_df <- genes_df[genes_df$CellType == model, ]
  input <- unique(na.omit(model_df$gene_name))
  input <- input[input != ""]

  message("Running Enrichr for ", model, ": ", length(input), " unique genes")

  if (length(input) == 0) {
    message("Skipping ", model, " because no valid gene names were found.")
    next
  }

  enriched <- enrichr(input, dbs_lib)

  # Keep the original per-cell-type output
  sig_df <- format_output(enriched, paste0(test, "_", model))

  # Also collect significant pathways into one combined TSV
  if (!is.null(sig_df) && nrow(sig_df) > 0) {
    sig_df$CellType <- model
    sig_df$Dataset <- test
    celltype_enrichr_list[[model]] <- sig_df
  }
}

if (length(celltype_enrichr_list) > 0) {
  celltype_enrichr_df <- rbindlist(celltype_enrichr_list, use.names = TRUE, fill = TRUE)

  # Put identifiers first
  setcolorder(
    celltype_enrichr_df,
    c("Dataset", "CellType", "Source", "Term", "Overlap", "P.value",
      "Adjusted.P.value", "Old.P.value", "Old.Adjusted.P.value",
      "Odds.Ratio", "Combined.Score", "Genes")
  )

  fwrite(
    celltype_enrichr_df,
    file = paste0(
      "/data/l2_bioinfo1/liq17/scTWAS/Functional_res/EnrichR/",
      test, "_CellTypes_sig_EnrichR_pathways.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    na = ""
  )

  message(
    "Combined significant cell-type Enrichr pathways written to: ",
    paste0(
      "/data/l2_bioinfo1/liq17/scTWAS/Functional_res/EnrichR/",
      test, "_CellTypes_sig_EnrichR_pathways.tsv"
    )
  )
} else {
  message("No significant cell-type Enrichr pathways were found.")
}




