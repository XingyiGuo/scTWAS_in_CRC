#!/usr/bin/env Rscript

# ==============================================================================
# CRC 51 TF motif logos from JASPAR2024
#
# Compatible with:
#   R 4.5.x
#   Bioconductor 3.22
#
# Outputs:
#   1. CRC51TF_JASPAR2024_selected_motifs.csv
#   2. CRC51TF_JASPAR2024_all_candidate_motifs.csv
#   3. CRC51TF_no_direct_JASPAR_motif.txt
#   4. CRC51TF_JASPAR2024_all_51_TF_motif_logos.pdf
#   5. CRC51TF_JASPAR2024_direct_motifs_only.pdf
#   6. individual_logos/*.pdf and *.png
#
# Important fix:
#   JASPAR2024 must be opened through its SQLite database connection.
#   Do NOT use getMatrixSet(JASPAR2024, opts), because JASPAR2024 is a function.
# ==============================================================================


# ==============================================================================
# 0. Install packages if needed
# ==============================================================================

if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

bioc_pkgs <- c(
    "JASPAR2024",
    "TFBSTools"
)

for (pkg in bioc_pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        BiocManager::install(
            pkg,
            ask = FALSE,
            update = FALSE
        )
    }
}

cran_pkgs <- c(
    "RSQLite",
    "ggplot2",
    "ggseqlogo",
    "gridExtra"
)

for (pkg in cran_pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        install.packages(pkg)
    }
}


# ==============================================================================
# 1. Load packages
# ==============================================================================

suppressPackageStartupMessages({
    library(JASPAR2024)
    library(TFBSTools)
    library(RSQLite)
    library(ggplot2)
    library(ggseqlogo)
    library(gridExtra)
})


# ==============================================================================
# 2. Define the 51 TFs
# ==============================================================================

TFs <- c(
    "ATF2",
    "ATF3",
    "CBX3",
    "CDX2",
    "CEBPB",
    "CNOT3",
    "E2F3",
    "E2F8",
    "EGR1",
    "ELF1",
    "EP300",
    "ETS2",
    "ETV5",
    "FOSL1",
    "FOXA2",
    "GATA6",
    "HNF4A",
    "HOXA6",
    "ICE1",
    "IRF1",
    "JUN",
    "JUND",
    "KLF5",
    "KMT2B",
    "MAX",
    "MED12",
    "MYC",
    "NIPBL",
    "NR1H2",
    "NR1H3",
    "PAF1",
    "POLR2A",
    "PPARG",
    "PROX1",
    "RAD21",
    "RELA",
    "REPIN1",
    "SIN3A",
    "SMC3",
    "SOX2",
    "SOX9",
    "SP1",
    "SRF",
    "TCF7L1",
    "TCF7L2",
    "TFAP4",
    "USF1",
    "VDR",
    "YY1",
    "ZBTB7B",
    "ZNF83"
)

stopifnot(length(TFs) == 51)

cat("Number of requested TFs:", length(TFs), "\n")


# ==============================================================================
# 3. Output directories
# ==============================================================================

OUTDIR <- "CRC51TF_JASPAR2024_motif_logos"
INDIVIDUAL_DIR <- file.path(OUTDIR, "individual_logos")

dir.create(
    OUTDIR,
    recursive = TRUE,
    showWarnings = FALSE
)

dir.create(
    INDIVIDUAL_DIR,
    recursive = TRUE,
    showWarnings = FALSE
)


# ==============================================================================
# 4. Open JASPAR2024 through its SQLite database
# ==============================================================================

jaspar <- JASPAR2024::JASPAR2024()

jaspar_db_path <- JASPAR2024::db(jaspar)

cat("\nJASPAR2024 SQLite database:\n")
cat(jaspar_db_path, "\n\n")

if (!file.exists(jaspar_db_path)) {
    stop(
        "Could not find the JASPAR2024 SQLite database at: ",
        jaspar_db_path
    )
}

jaspar_con <- RSQLite::dbConnect(
    RSQLite::SQLite(),
    jaspar_db_path
)


# ==============================================================================
# 5. Retrieve current human CORE PFMs
# ==============================================================================

opts <- list(
    collection   = "CORE",
    tax_group    = "vertebrates",
    species      = "Homo sapiens",
    matrixtype   = "PFM",
    all_versions = FALSE
)

cat("Retrieving JASPAR2024 human CORE motifs...\n")

motifs <- TFBSTools::getMatrixSet(
    x = jaspar_con,
    opts = opts
)

cat(
    "Number of human CORE motifs retrieved:",
    length(motifs),
    "\n\n"
)

if (length(motifs) == 0) {
    RSQLite::dbDisconnect(jaspar_con)
    stop("No JASPAR motifs were retrieved.")
}


# ==============================================================================
# 6. Extract motif metadata
# ==============================================================================

motif_info <- data.frame(
    motif_index = seq_along(motifs),

    motif_id = vapply(
        motifs,
        function(x) TFBSTools::ID(x),
        character(1)
    ),

    motif_name = vapply(
        motifs,
        function(x) TFBSTools::name(x),
        character(1)
    ),

    stringsAsFactors = FALSE
)

motif_info$motif_name_upper <- toupper(motif_info$motif_name)

cat("First retrieved motifs:\n")
print(head(motif_info, 10))


# ==============================================================================
# 7. Calculate total motif information content
#
# Used only to choose one representative motif if JASPAR contains more than
# one current exact-name motif for the same TF.
# ==============================================================================

calculate_total_ic <- function(pfm) {

    pfm <- as.matrix(pfm)
    required_bases <- c("A", "C", "G", "T")

    if (!all(required_bases %in% rownames(pfm))) {
        return(NA_real_)
    }

    pfm <- pfm[required_bases, , drop = FALSE]

    column_sums <- colSums(pfm)
    valid <- column_sums > 0

    pfm <- pfm[, valid, drop = FALSE]
    column_sums <- column_sums[valid]

    if (ncol(pfm) == 0) {
        return(NA_real_)
    }

    prob <- sweep(
        pfm,
        2,
        column_sums,
        "/"
    )

    position_ic <- apply(
        prob,
        2,
        function(p) {
            p <- p[p > 0]
            2 + sum(p * log2(p))
        }
    )

    sum(position_ic)
}


motif_info$total_information_bits <- vapply(
    seq_along(motifs),
    function(i) {
        pfm <- TFBSTools::Matrix(motifs[[i]])
        calculate_total_ic(pfm)
    },
    numeric(1)
)


# ==============================================================================
# 8. Exact-match the 51 factors to JASPAR motif names
#
# We intentionally do not force a motif from another TF or heterodimer onto
# factors without their own exact JASPAR motif.
# ==============================================================================

candidate_list <- lapply(
    TFs,
    function(tf) {

        idx <- which(
            motif_info$motif_name_upper == toupper(tf)
        )

        if (length(idx) == 0) {
            return(NULL)
        }

        result <- motif_info[
            idx,
            ,
            drop = FALSE
        ]

        result$TF <- tf
        result
    }
)

candidate_list <- Filter(Negate(is.null), candidate_list)

if (length(candidate_list) > 0) {

    candidate_df <- do.call(
        rbind,
        candidate_list
    )

} else {

    candidate_df <- data.frame(
        motif_index = integer(),
        motif_id = character(),
        motif_name = character(),
        motif_name_upper = character(),
        total_information_bits = numeric(),
        TF = character(),
        stringsAsFactors = FALSE
    )
}

write.csv(
    candidate_df,
    file.path(
        OUTDIR,
        "CRC51TF_JASPAR2024_all_candidate_motifs.csv"
    ),
    row.names = FALSE
)


# ==============================================================================
# 9. Select one representative motif per TF
#
# Selection rule:
#   - exact TF-name match
#   - current JASPAR versions only (all_versions = FALSE)
#   - if multiple exact motifs remain, choose highest total information content
# ==============================================================================

selection_list <- lapply(
    TFs,
    function(tf) {

        candidates <- candidate_df[
            candidate_df$TF == tf,
            ,
            drop = FALSE
        ]

        if (nrow(candidates) == 0) {

            return(
                data.frame(
                    TF = tf,
                    status = "NO_DIRECT_JASPAR_MOTIF",
                    n_candidates = 0L,
                    selected_motif_id = NA_character_,
                    selected_motif_name = NA_character_,
                    selected_motif_index = NA_integer_,
                    total_information_bits = NA_real_,
                    stringsAsFactors = FALSE
                )
            )
        }

        candidates <- candidates[
            order(
                -candidates$total_information_bits,
                candidates$motif_id,
                na.last = TRUE
            ),
            ,
            drop = FALSE
        ]

        selected <- candidates[1, ]

        data.frame(
            TF = tf,
            status = "FOUND",
            n_candidates = nrow(candidates),
            selected_motif_id = selected$motif_id,
            selected_motif_name = selected$motif_name,
            selected_motif_index = as.integer(selected$motif_index),
            total_information_bits = selected$total_information_bits,
            stringsAsFactors = FALSE
        )
    }
)

selection_df <- do.call(
    rbind,
    selection_list
)

write.csv(
    selection_df,
    file.path(
        OUTDIR,
        "CRC51TF_JASPAR2024_selected_motifs.csv"
    ),
    row.names = FALSE
)


# ==============================================================================
# 10. Print matching summary
# ==============================================================================

found_TFs <- selection_df$TF[
    selection_df$status == "FOUND"
]

missing_TFs <- selection_df$TF[
    selection_df$status != "FOUND"
]

cat("\n============================================================\n")
cat("JASPAR matching summary\n")
cat("============================================================\n")

cat("Total requested TFs:", length(TFs), "\n")
cat("Direct JASPAR motifs found:", length(found_TFs), "\n")
cat("No direct JASPAR motif:", length(missing_TFs), "\n\n")

if (length(missing_TFs) > 0) {
    cat("TFs without a direct JASPAR2024 motif:\n")
    print(missing_TFs)
}

writeLines(
    missing_TFs,
    file.path(
        OUTDIR,
        "CRC51TF_no_direct_JASPAR_motif.txt"
    )
)


# ==============================================================================
# 11. Save TFs with multiple exact candidate motifs
# ==============================================================================

multiple_df <- selection_df[
    selection_df$n_candidates > 1,
    ,
    drop = FALSE
]

if (nrow(multiple_df) > 0) {

    cat("\nTFs with multiple exact JASPAR motifs:\n")

    print(
        multiple_df[
            ,
            c(
                "TF",
                "n_candidates",
                "selected_motif_id"
            )
        ]
    )

    write.csv(
        multiple_df,
        file.path(
            OUTDIR,
            "CRC51TF_multiple_motif_TFs.csv"
        ),
        row.names = FALSE
    )
}


# ==============================================================================
# 12. Sequence-logo plotting function
#
# method = "bits" creates the same type of logo as the example figure:
# stacked A/C/G/T letters with information content on a 0-2 bit scale.
# ==============================================================================

make_logo_plot <- function(
    motif_object,
    tf_name,
    motif_id
) {

    pfm <- as.matrix(
        TFBSTools::Matrix(motif_object)
    )

    required_bases <- c("A", "C", "G", "T")

    if (!all(required_bases %in% rownames(pfm))) {
        stop(
            "Motif ",
            motif_id,
            " does not contain standard A/C/G/T row names."
        )
    }

    pfm <- pfm[
        required_bases,
        ,
        drop = FALSE
    ]

    ggseqlogo::ggseqlogo(
        pfm,
        method = "bits",
        seq_type = "dna"
    ) +
        labs(
            title = tf_name,
            subtitle = motif_id,
            x = "Position",
            y = "Bits"
        ) +
        scale_y_continuous(
            limits = c(0, 2),
            breaks = c(0, 0.5, 1, 1.5, 2),
            expand = c(0, 0)
        ) +
        theme_classic(
            base_size = 10
        ) +
        theme(
            plot.title = element_text(
                face = "bold",
                size = 12,
                hjust = 0.5
            ),
            plot.subtitle = element_text(
                size = 8,
                hjust = 0.5
            ),
            axis.title.x = element_text(size = 9),
            axis.title.y = element_text(size = 9),
            axis.text.x = element_text(size = 7),
            axis.text.y = element_text(size = 7),
            plot.margin = margin(5, 5, 5, 5)
        )
}


# ==============================================================================
# 13. Placeholder panel for factors without a direct JASPAR motif
# ==============================================================================

make_missing_plot <- function(tf_name) {

    ggplot() +
        annotate(
            "text",
            x = 0.5,
            y = 0.60,
            label = tf_name,
            fontface = "bold",
            size = 4
        ) +
        annotate(
            "text",
            x = 0.5,
            y = 0.40,
            label = "No direct JASPAR2024 motif",
            size = 3
        ) +
        xlim(0, 1) +
        ylim(0, 1) +
        theme_void() +
        theme(
            panel.border = element_rect(
                fill = NA,
                linewidth = 0.3
            )
        )
}


# ==============================================================================
# 14. Generate individual motif logos
# ==============================================================================

plot_list <- vector(
    "list",
    length(TFs)
)

names(plot_list) <- TFs


for (tf in TFs) {

    cat("Processing:", tf, "\n")

    selected <- selection_df[
        selection_df$TF == tf,
        ,
        drop = FALSE
    ]

    if (
        nrow(selected) == 0 ||
        selected$status[1] != "FOUND"
    ) {

        plot_list[[tf]] <- make_missing_plot(tf)
        next
    }

    motif_index <- selected$selected_motif_index[1]
    motif_id <- selected$selected_motif_id[1]

    motif_object <- motifs[[motif_index]]

    p <- make_logo_plot(
        motif_object = motif_object,
        tf_name = tf,
        motif_id = motif_id
    )

    plot_list[[tf]] <- p

    ggsave(
        filename = file.path(
            INDIVIDUAL_DIR,
            paste0(
                tf,
                "_",
                motif_id,
                "_motif_logo.pdf"
            )
        ),
        plot = p,
        width = 5,
        height = 3,
        units = "in"
    )

    ggsave(
        filename = file.path(
            INDIVIDUAL_DIR,
            paste0(
                tf,
                "_",
                motif_id,
                "_motif_logo.png"
            )
        ),
        plot = p,
        width = 5,
        height = 3,
        units = "in",
        dpi = 600
    )
}


# ==============================================================================
# 15. Combined multi-page PDF with all 51 factors
#
# 3 columns x 4 rows = 12 panels/page
# ==============================================================================

combined_pages <- gridExtra::marrangeGrob(
    grobs = plot_list,
    nrow = 4,
    ncol = 3,
    top = NULL
)

ggsave(
    filename = file.path(
        OUTDIR,
        "CRC51TF_JASPAR2024_all_51_TF_motif_logos.pdf"
    ),
    plot = combined_pages,
    width = 12,
    height = 10,
    units = "in"
)


# ==============================================================================
# 16. Combined PDF containing only factors with direct JASPAR motifs
# ==============================================================================

found_plot_list <- plot_list[
    found_TFs
]

if (length(found_plot_list) > 0) {

    found_pages <- gridExtra::marrangeGrob(
        grobs = found_plot_list,
        nrow = 4,
        ncol = 3,
        top = NULL
    )

    ggsave(
        filename = file.path(
            OUTDIR,
            "CRC51TF_JASPAR2024_direct_motifs_only.pdf"
        ),
        plot = found_pages,
        width = 12,
        height = 10,
        units = "in"
    )
}


# ==============================================================================
# 17. Close database connection
# ==============================================================================

RSQLite::dbDisconnect(jaspar_con)


# ==============================================================================
# 18. Final report
# ==============================================================================

cat("\n")
cat("============================================================\n")
cat("Finished\n")
cat("============================================================\n")

cat("Requested TFs:   ", length(TFs), "\n")
cat("Motifs found:    ", length(found_TFs), "\n")
cat("No direct motif: ", length(missing_TFs), "\n")

cat(
    "Output directory:",
    normalizePath(OUTDIR),
    "\n"
)

cat("\nMain outputs:\n")
cat("  CRC51TF_JASPAR2024_selected_motifs.csv\n")
cat("  CRC51TF_JASPAR2024_all_candidate_motifs.csv\n")
cat("  CRC51TF_no_direct_JASPAR_motif.txt\n")
cat("  CRC51TF_JASPAR2024_all_51_TF_motif_logos.pdf\n")
cat("  CRC51TF_JASPAR2024_direct_motifs_only.pdf\n")
cat("  individual_logos/\n")

cat("\nDone.\n")
