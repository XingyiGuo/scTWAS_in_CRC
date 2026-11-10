#!/usr/bin/env Rscript
# Export the full Enrichr results and the pathways with adjusted P < 0.05.

args <- commandArgs(trailingOnly = TRUE)
usage <- paste(
  "Usage: Rscript 03_enrich_druggable_genes_enrichr.R",
  "<63_druggable_genes_summary.csv> <output.tsv> [library]"
)
if (length(args) == 1L && args[1L] %in% c("-h", "--help")) {
  cat(usage, "\n")
  quit(save = "no")
}
if (length(args) < 2L || length(args) > 3L) {
  stop(usage, call. = FALSE)
}

summary <- read.csv(args[1L], colClasses = "character", check.names = FALSE)
if (!"gene_name" %in% names(summary)) {
  stop("The input CSV must contain a gene_name column.", call. = FALSE)
}
genes <- unique(trimws(summary$gene_name))
genes <- genes[!is.na(genes) & nzchar(genes)]
if (!length(genes)) {
  stop("No gene symbols found in the input CSV.", call. = FALSE)
}
database <- if (length(args) == 3L) args[3L] else "KEGG_2021_Human"

enrichR::setEnrichrSite("Enrichr")
databases <- enrichR::listEnrichrDbs()
if (!database %in% databases$libraryName) {
  stop("Enrichr library unavailable: ", database, call. = FALSE)
}
results <- enrichR::enrichr(genes, databases = database)[[database]]
if (!is.data.frame(results) || !all(c("Term", "Adjusted.P.value") %in% names(results))) {
  stop("Enrichr did not return an enrichment table for ", database, call. = FALSE)
}

results$Library <- rep(database, nrow(results))
results <- results[order(results$Adjusted.P.value, results$Term), , drop = FALSE]
significant <- results[which(results$Adjusted.P.value < 0.05), , drop = FALSE]
output <- args[2L]
significant_output <- paste0(tools::file_path_sans_ext(output), "_significant.tsv")
dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
write.table(results, file = output, sep = "\t", row.names = FALSE, quote = FALSE, na = "")
write.table(significant, file = significant_output, sep = "\t", row.names = FALSE,
            quote = FALSE, na = "")
message("Submitted ", length(genes), " gene symbols to ", database,
        " and wrote ", nrow(results), " terms to ", output, ".")
message("Wrote ", nrow(significant), " pathways with Enrichr adjusted P < 0.05 to ",
        significant_output, ".")
