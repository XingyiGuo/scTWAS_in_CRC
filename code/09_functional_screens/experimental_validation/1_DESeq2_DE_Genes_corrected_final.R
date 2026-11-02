# Differential expression analysis for multiple cell lines and KD genes
# DESeq2: KD vs NT

library(readxl)
library(DESeq2)
library(dplyr)
library(data.table)

# Set the random seed
set.seed(20260909)

source("/data/l2_bioinfo1/liq17/LQ_frameworks/GenoPhenoBasic.R")
# ============================================================
# Input
# ============================================================

META_FILE <- "/data/l2_bioinfo1/liq17/scTWAS/data/Experiments/15174-XG_project_summary.xlsx"
COUNT_FILE <- "/data/l2_bioinfo1/liq17/scTWAS/data/Experiments/salmon.merged.gene_counts.tsv"
OUT_DIR <- "/data/l2_bioinfo1/liq17/scTWAS/Experiments_res/RNA_seq_DE"
GTF_FILE <- "/data/l2_bioinfo1/liq17/LQ_references/gencode/gencode.v26.annotation.gtf"

gtf_v26 <- get_gene_annotation_all(GTF_FILE)

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

sample_meta_df <- as.data.frame(read_excel(META_FILE, sheet = "Sheet1"))
count_df <- as.data.frame(fread(COUNT_FILE))

CellTypes <- c("HCT116", "SW480")
Genes <- c("LAMC2", "METTL1", "NOS3", "PC", "PILRB", "PPP2R1B", "PYGL")

# ENSG -> gene symbol annotation
gene_anno <- count_df %>%
    select(gene_id, gene_name) %>%
    distinct(gene_id, .keep_all = TRUE)

# Store significant results and summary from all comparisons
sig_res_list <- list()
sig_summary_list <- list()

# ============================================================
# DESeq2 analysis
# ============================================================

for (ct in CellTypes) {
    for (gene in Genes) {

        message("Running: ", ct, " / ", gene)

        # ----------------------------------------------------
        # 1. Select NT and KD samples
        # ----------------------------------------------------
        pairs_meta_df <- sample_meta_df %>%
            filter(Cell == ct, Gene %in% c("NT", gene))

        # Make NT the reference group
        pairs_meta_df$Gene <- factor(pairs_meta_df$Gene, levels = c("NT", gene))

        # Print the exact samples used in this comparison
        message("  NT samples: ", paste(pairs_meta_df$Filename[pairs_meta_df$Gene == "NT"], collapse = ", "))
        message("  KD samples: ", paste(pairs_meta_df$Filename[pairs_meta_df$Gene == gene], collapse = ", "))

        # ----------------------------------------------------
        # 2. Prepare raw count matrix
        # ----------------------------------------------------
        count_pairs_df <- count_df[, c("gene_id", pairs_meta_df$Filename)]

        rownames(count_pairs_df) <- count_pairs_df$gene_id
        count_pairs_df$gene_id <- NULL

		count_pairs_df <- as.matrix(count_pairs_df)
		storage.mode(count_pairs_df) <- "numeric"

		# Check for invalid values
		stopifnot(!anyNA(count_pairs_df),all(count_pairs_df >= 0))

		# Convert Salmon estimated counts to integer counts for DESeq2
		count_pairs_df <- round(count_pairs_df)
		storage.mode(count_pairs_df) <- "integer"

        # ----------------------------------------------------
        # 3. Prepare sample metadata
        # ----------------------------------------------------
        coldata <- data.frame(
            Gene = pairs_meta_df$Gene,
            row.names = pairs_meta_df$Filename
        )

        # Count columns and metadata rows must have identical order
        stopifnot(identical(rownames(coldata), colnames(count_pairs_df)))

        # ----------------------------------------------------
        # 4. Construct DESeq2 object and filter low-count genes
        # ----------------------------------------------------
        dds <- DESeqDataSetFromMatrix(
            countData = count_pairs_df,
            colData = coldata,
            design = ~ Gene
        )

        # Keep genes with >=5 counts in at least 2 samples (2 replicates per group)
        keep <- rowSums(counts(dds) >= 5) >= 2
        dds <- dds[keep, ]

        # ----------------------------------------------------
        # 5. Fit DESeq2 model
        # ----------------------------------------------------
        dds2 <- DESeq(dds)

        # Dispersion plot
        pdf(
            file.path(OUT_DIR, paste0(ct, "_", gene, "_dispersion_plot.pdf")),
            width = 6,
            height = 6
        )
        plotDispEsts(dds2)
        dev.off()

        # ----------------------------------------------------
        # 6. Normalized counts and group means
        # ----------------------------------------------------
        norm_counts <- counts(dds2, normalized = TRUE)

        NT_mean <- rowMeans(norm_counts[, dds2$Gene == "NT", drop = FALSE])

        KD_mean <- rowMeans(norm_counts[, dds2$Gene == gene, drop = FALSE])

        Dispersion <- mcols(dds2)$dispersion

        count_info <- cbind(
            norm_counts,
            NT = NT_mean,
            KD = KD_mean,
            Dispersion = Dispersion
        )

        # ----------------------------------------------------
        # 7. Differential expression: KD vs NT
        #    Positive log2FC = higher expression in KD
        #    Negative log2FC = lower expression in KD
        # ----------------------------------------------------
        res_gene <- results(
            dds2,
            contrast = c("Gene", gene, "NT")
        )

        res_gene <- cbind(count_info, as.data.frame(res_gene))
        res_gene <- as.data.frame(res_gene)

        # Add ENSG and gene symbol explicitly
        res_gene$gene_id <- rownames(res_gene)
        res_gene <- res_gene %>%
            left_join(gene_anno, by = "gene_id") %>%
            relocate(gene_id, gene_name) %>%
            arrange(is.na(padj), padj)

        # Save all DESeq2 results
        write.csv(
            res_gene,
            file.path(OUT_DIR, paste0(ct, "_", gene, "_DESeq2_res.csv")),
            row.names = FALSE
        )

        # ----------------------------------------------------
        # 8. Collect significant DE genes
        #
        # Significance rule:
        #   1) If padj is available, use padj < 0.05.
        #   2) If padj is NA, fall back to:
        #          pvalue < 0.05 / number of tested genes
        #
        # Effect-size requirement:
        #   abs(log2FoldChange) > 0.3
        # ----------------------------------------------------
        n_tests <- nrow(res_gene)
        bonferroni_cutoff <- 0.05 / n_tests

        sig_res <- res_gene %>%
            filter(
                (
                    !is.na(padj) &
                    padj < 0.05
                ) |
                (
                    is.na(padj) &
                    !is.na(pvalue) &
                    pvalue < bonferroni_cutoff
                ),
                !is.na(log2FoldChange),
                abs(log2FoldChange) > 0.3
            ) %>%
            mutate(
                Cell = ct,
                KD_Gene = gene,
                Significance_method = dplyr::if_else(
                    !is.na(padj),
                    "BH_FDR_padj",
                    "Bonferroni_rawP"
                ),
                Direction = dplyr::if_else(
                    log2FoldChange > 0.3,
                    "Up",
                    "Down"
                ),
                N_tests = n_tests,
                Bonferroni_cutoff = bonferroni_cutoff
            ) %>%
            # IMPORTANT:
            # Keep only summary/statistical columns here.
            # Individual sample normalized-count columns are intentionally
            # excluded from the combined significant-results table.
            dplyr::select(
                Cell,
                KD_Gene,
                gene_id,
                gene_name,
                Significance_method,
                Direction,
                NT,
                KD,
                Dispersion,
                baseMean,
                log2FoldChange,
                lfcSE,
                stat,
                pvalue,
                padj,
                N_tests,
                Bonferroni_cutoff
            )

        key <- paste(ct, gene, sep = "_")

        # Only add non-empty significant results to the combined-results list.
        # Zero-hit comparisons are still retained in the summary below.
        if (nrow(sig_res) > 0) {
            sig_res_list[[key]] <- sig_res
        }

        # Record summary even when there are zero significant genes
        sig_summary_list[[key]] <- data.frame(
            Cell = ct,
            KD_Gene = gene,
            N_tests = n_tests,
            Bonferroni_cutoff = bonferroni_cutoff,
            N_sig = nrow(sig_res),
            N_up = sum(sig_res$Direction == "Up"),
            N_down = sum(sig_res$Direction == "Down")
        )

        message(
            "Finished ", ct, " / ", gene,
            " | significant DE genes = ", nrow(sig_res)
        )
    }
}

# ============================================================
# Combine significant results from all comparisons
# ============================================================

sig_res_all <- bind_rows(sig_res_list)

# Final safety check: the combined significant-results file must not contain
# individual sample normalized-count columns.
sig_res_all <- sig_res_all %>%
    dplyr::select(
        Cell,
        KD_Gene,
        gene_id,
        gene_name,
        Significance_method,
        Direction,
        NT,
        KD,
        Dispersion,
        baseMean,
        log2FoldChange,
        lfcSE,
        stat,
        pvalue,
        padj,
        N_tests,
        Bonferroni_cutoff
    )


# If gene_name still contains an ENSG identifier, annotate it with gtf_v26
if (nrow(sig_res_all) > 0L) {
    idx_ensg <- !is.na(sig_res_all$gene_name) &
        grepl("^ENSG", sig_res_all$gene_name)

    if (any(idx_ensg)) {
        ensg <- sub("\\.[0-9]+$", "", sig_res_all$gene_name[idx_ensg])
        new_name <- gtf_v26$genename[match(ensg, gtf_v26$geneid)]

        current_name <- sig_res_all$gene_name[idx_ensg]
        replace_idx <- !is.na(new_name) & nzchar(new_name)
        current_name[replace_idx] <- new_name[replace_idx]
        sig_res_all$gene_name[idx_ensg] <- current_name
    }
}

write.csv(
    sig_res_all,
    file.path(OUT_DIR, "HCT116_SW480_7KD_DESeq2_significant_results.csv"),
    row.names = FALSE
)

# Summary: number of significant up/down genes per comparison
# Includes comparisons with zero significant genes.
sig_summary <- bind_rows(sig_summary_list)

write.csv(
    sig_summary,
    file.path(OUT_DIR, "HCT116_SW480_7KD_DESeq2_significant_summary.csv"),
    row.names = FALSE
)

print(sig_summary)
message("All DESeq2 analyses finished.")
