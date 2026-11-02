

suppressPackageStartupMessages({
  library(Seurat)        # >= 4.3 (works with v5 too)
  library(Matrix)
  library(dplyr)
  library(stringr)
  library(harmony)
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(R.utils)
  library(future)
})


set.seed(123)

# ---- paths ------------------------------------------------------------------
# Locate script directory robustly (works under Rscript and source())
.find_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  hit  <- grep("--file=", args, value = TRUE)
  if (length(hit) > 0) return(dirname(normalizePath(sub("--file=", "", hit[1]))))
  if (!is.null(sys.frames()) && length(sys.frames()) > 0) {
    f <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
    if (!is.null(f)) return(dirname(normalizePath(f)))
  }
  getwd()
}
SCRIPT_DIR <- .find_script_dir()

ROOT <- "/data/l2_bioinfo1/lyul1/single_cell_data/IBD/"
RESULTS_DIR <- file.path(ROOT, "Results")
dir.create(RESULTS_DIR, showWarnings = FALSE, recursive = TRUE)

GSE266_DIR <- "/data/l2_bioinfo1/lyul1/single_cell_data/IBD/raw_counts"
HTAN_DIR   <- "/data/l2_bioinfo1/HTAN_scRNA" # "/data/l2_bioinfo1/HTAN_scRNA"

message("ROOT       = ", ROOT)
message("GSE266_DIR = ", GSE266_DIR)
message("HTAN_DIR   = ", HTAN_DIR)


#============================================================================
# 1) Load GSE266546 -- ASCENDING COLON, Control + Normal_CD only
# ============================================================================
message("\n[1/5] Loading GSE266546 (Vanderbilt CD/Control, AC)...")

load_gse266 <- function(counts_gz, annot_gz, source_tag) {
  ann <- fread(annot_gz, sep = "\t", header = TRUE, data.table = FALSE)
  
  bc_col <- intersect(c("Barcodes", "Barcode", "barcode", "Cell", "cell", "Cell_Name", "V1"),
                      colnames(ann))[1]
  if (is.na(bc_col)) bc_col <- colnames(ann)[1]
  
  ann$.barcode <- ann[[bc_col]]
  rownames(ann) <- ann$.barcode
  
  cnt <- fread(counts_gz, sep = "\t", header = FALSE, data.table = FALSE)
  
  rownames(cnt) <- cnt[[1]]
  cnt[[1]] <- NULL
  
  if (ncol(cnt) != nrow(ann)) {
    stop(
      "After removing gene column, count columns (", ncol(cnt),
      ") do not match annotation rows (", nrow(ann), ") for ",
      basename(counts_gz)
    )
  }
  
  colnames(cnt) <- ann$.barcode
  cnt <- as(as.matrix(cnt), "dgCMatrix")
  
  ann$source <- source_tag
  list(counts = cnt, meta = ann)
}


ac_endo <- load_gse266(
  file.path(GSE266_DIR, "GSE266546_Raw_counts_AC_Endoscopy.txt.gz"),
  file.path(GSE266_DIR, "GSE266546_Annotation_AC_Endoscopy.txt.gz"),
  source_tag = "AC_Endoscopy"
)
ac_surg <- load_gse266( 
  file.path(GSE266_DIR, "GSE266546_Raw_counts_AC_Surgical.txt.gz"),
  file.path(GSE266_DIR, "GSE266546_Annotation_AC_Surgical.txt.gz"),
  source_tag = "AC_Surgical"
)

ti_endo <- load_gse266(
  file.path(GSE266_DIR, "GSE266546_Raw_counts_TI_Endoscopy.txt.gz"),
  file.path(GSE266_DIR, "GSE266546_Annotation_TI_Endoscopy.txt.gz"),
  source_tag = "TI_Endoscopy"
)
ti_surg <- load_gse266( 
  file.path(GSE266_DIR, "GSE266546_Raw_counts_TI_Surgical.txt.gz"),
  file.path(GSE266_DIR, "GSE266546_Annotation_TI_Surgical.txt.gz"),
  source_tag = "TI_Surgical"
)

# Stack the two AC compartments
# Check if the genes are the same.
gene_lists <- list(
  ac_endo = rownames(ac_endo$counts),
  ac_surg = rownames(ac_surg$counts),
  ti_endo = rownames(ti_endo$counts),
  ti_surg = rownames(ti_surg$counts)
)

# Number of genes in each
sapply(gene_lists, length)

# Are gene sets identical, ignoring order?
sapply(gene_lists, function(x) setequal(gene_lists$ac_endo, x))

#common_genes  <- intersect(rownames(ac_endo$counts), rownames(ac_surg$counts))
common_genes <- Reduce(intersect, list(
  rownames(ac_endo$counts),
  rownames(ac_surg$counts),
  rownames(ti_endo$counts),
  rownames(ti_surg$counts)
))

length(common_genes)

# Keep union of all genes across GSE266546 matrices
all_genes <- Reduce(union, list(
  rownames(ac_endo$counts),
  rownames(ac_surg$counts),
  rownames(ti_endo$counts),
  rownames(ti_surg$counts)
))

align_to_genes <- function(mat, genes) {
  missing_genes <- setdiff(genes, rownames(mat))
  
  if (length(missing_genes) > 0) {
    zero_mat <- Matrix::Matrix(
      0,
      nrow = length(missing_genes),
      ncol = ncol(mat),
      sparse = TRUE
    )
    rownames(zero_mat) <- missing_genes
    colnames(zero_mat) <- colnames(mat)
    
    mat <- rbind(mat, zero_mat)
  }
  
  mat[genes, , drop = FALSE]
}

ac_endo_mat <- align_to_genes(ac_endo$counts, all_genes)
ac_surg_mat <- align_to_genes(ac_surg$counts, all_genes)
ti_endo_mat <- align_to_genes(ti_endo$counts, all_genes)
ti_surg_mat <- align_to_genes(ti_surg$counts, all_genes)

gse266_counts <- cbind(
  ac_endo_mat,
  ac_surg_mat,
  ti_endo_mat,
  ti_surg_mat
)

gse266_meta <- bind_rows(
  ac_endo$meta,
  ac_surg$meta,
  ti_endo$meta,
  ti_surg$meta
)

rownames(gse266_meta) <- gse266_meta$.barcode

gse266_meta <- gse266_meta[colnames(gse266_counts), , drop = FALSE]

stopifnot(
  ncol(gse266_counts) == nrow(gse266_meta),
  all(colnames(gse266_counts) == rownames(gse266_meta))
)
length(all_genes)
dim(gse266_counts)

message(
  "  GSE266546 combined matrix: ",
  nrow(gse266_counts), " genes x ",
  ncol(gse266_counts), " cells"
)

message(
  "  GSE266546 metadata: ",
  length(unique(gse266_meta$PatientID)), " patients across ",
  length(unique(gse266_meta$SampleID)), " samples"
)

message("  Cells by tissue:")
print(table(gse266_meta$Tissue, useNA = "ifany"))

message("  Cells by histology:")
print(table(gse266_meta$Histology, useNA = "ifany"))

message("  Cells by source:")
print(table(gse266_meta$source, useNA = "ifany"))


# gse266_counts <- cbind(ac_endo$counts[common_genes, ],
#                        ac_surg$counts[common_genes, ],
#                        ti_endo$counts[common_genes, ],
#                        ti_surg$counts[common_genes, ]
#                        )
# gse266_meta   <- bind_rows(ac_endo$meta, ac_surg$meta, ti_endo$meta, ti_surg$meta)
# bind_rows drops rownames -- restore from the explicit .barcode column
# rownames(gse266_meta) <- gse266_meta$.barcode
# Make sure order matches the count matrix exactly
# gse266_meta   <- gse266_meta[colnames(gse266_counts), , drop = FALSE]
# stopifnot(ncol(gse266_counts) == nrow(gse266_meta),
#           all(colnames(gse266_counts) == rownames(gse266_meta)))


# Pull patient + condition out of the cell/sample name if not already present.
# Sample names look like: GCA017-Control-AC-Endoscopy ; cells often inherit
# "<sample>_<barcode>". The annotation file usually has a `Sample` column.

gse266_meta$sample_condition  <- gse266_meta$Histology
samp_col <- "SampleID"

# Filter: colon (AC) controls and "Normal_CD" (CD patients but the sampled
# tissue is histologically normal -- the user's "Control and Normal tissues").
keep_cells <- gse266_meta$sample_condition %in% c("Control", "Normal_CD")
gse266_counts <- gse266_counts[, keep_cells, drop = FALSE]
gse266_meta   <- gse266_meta[keep_cells, , drop = FALSE]
gse266_meta$CellType_cb<- gse266_meta$Celltypes

message("  GSE266546 retained cells after normal/control filter: ", ncol(gse266_counts),
        " across ", length(unique(gse266_meta[[samp_col]])), " samples ", length(unique(gse266_meta[["PatientID"]])), " patients")

message(
  "    ", nrow(gse266_counts), " genes x ",
  ncol(gse266_counts), " cells"
)
message(
  "    ", length(unique(gse266_meta[[samp_col]])), " samples, ",
  length(unique(gse266_meta$PatientID)), " patients"
)
message("    Cells by sample_condition:")
print(table(gse266_meta$sample_condition, useNA = "ifany"))
message("    Cells by tissue:")
print(table(gse266_meta$Tissue, useNA = "ifany"))
message("    Cells by source:")
print(table(gse266_meta$source, useNA = "ifany"))


obj266 <- CreateSeuratObject(
  counts    = gse266_counts,
  meta.data = gse266_meta,
  project   = "GSE266546",
  min.cells = 3, min.features = 200
)
obj266$study <- "GSE266546"
obj266$batch <- "GSE266546"

saveRDS(obj266, file.path(RESULTS_DIR, "obj266_preqc.rds"))
write.csv(obj266@meta.data, file.path(RESULTS_DIR, "obj266_metadata_preqc.csv"), row.names = F)

rm(ac_endo, ac_surg, gse266_counts, ti_endo, ti_surg, gse266_meta); gc()


# 3) filter to preferred site AC
preferred_site <- "Ascending Colon"

patient_site <- unique(obj266@meta.data[, c("PatientID", "Tissue")])

patients_with_preferred <- unique(
  patient_site$PatientID[patient_site$Tissue == preferred_site]
)

obj266 <- subset(
  obj266,
  subset = Tissue == preferred_site | !(PatientID %in% patients_with_preferred)
)

table(obj266$Tissue, useNA = "ifany")
length(unique(obj266$PatientID))

saveRDS(obj266, file.path(RESULTS_DIR, "obj266_preqc_pref_AC.rds"))

write.csv(obj266@meta.data, file.path(RESULTS_DIR, "obj266_metadata_preqc_pref_AC.csv"), row.names = F)

message("  GSE266546 after preferred-site filter:")
message(
  "    ", nrow(obj266), " genes x ",
  ncol(obj266), " cells"
)
message(
  "    ", length(unique(obj266$SampleID)), " samples, ",
  length(unique(obj266$PatientID)), " patients"
)
message("    Cells by tissue:")
print(table(obj266$Tissue, useNA = "ifany"))

site_per_patient_after <- tapply(
  obj266$Tissue,
  obj266$PatientID,
  function(x) length(unique(x))
)
message("    Tissue sites per patient after:")
print(table(site_per_patient_after, useNA = "ifany"))


# ============================================================================
# 2) Load HTAN Islam 2024 (Vanderbilt Colon MAP)
#    User has manually downloaded files from
#    https://humantumoratlas.org/publications/hta11_2024_nature_mirazul-islam
#    into HTAN_DIR. Two layouts are auto-detected:
#       (a) per-sample 10x folders containing matrix.mtx.gz etc.
#       (b) one or more *.h5ad files (Level 4 HTAN)
# ============================================================================
message("\n[2/5] Loading HTAN Islam 2024 (Vanderbilt Colon MAP)...")

# Discover what's there
h5ad_files <- list.files(HTAN_DIR, pattern = "\\.h5ad$",
                         recursive = TRUE, full.names = TRUE)

message("  found ", length(h5ad_files))
if (length(h5ad_files) == 0 ) {
  stop("HTAN_DIR is empty. See 01_download_data.sh for what to put there.")
}



# ---- loader (b): h5ad files -----------------------------------------------
# Uses the CRAN `anndata` package (which wraps reticulate + python anndata).
# Alternative: zellkonverter (Bioc) + as.Seurat() on the SingleCellExperiment.

load_one_h5ad <- function(f) {
  #py <- "/home/lyul1/venvs/r-h5ad/bin/python"
  #py_path <- paste(
  #  "/home/lyul1/venvs/r-h5ad/lib/python3.12/site-packages",
  #  "/cvmfs/soft.computecanada.ca/easybuild/software/2023/x86-64-v3/Compiler/gcccore/scipy-stack/2025a/lib/python3.12/site-packages",
  #  sep = ":"
  #)
  
  py <- "/mnt_home/venvs/r-h5ad/bin/python"
  
  py_path <- paste(
    "/mnt_home/venvs/r-h5ad/lib/python3.12/site-packages",
    "/cvmfs/soft.computecanada.ca/easybuild/software/2023/x86-64-v3/Compiler/gcccore/scipy-stack/2025a/lib/python3.12/site-packages",
    sep = ":"
  )
  outdir <- tempfile("h5ad_")
  dir.create(outdir)
  
  script <- tempfile(fileext = ".py")
  writeLines(c(
    "import sys",
    "import anndata",
    "import pandas as pd",
    "import scipy.sparse as sp",
    "import scipy.io as sio",
    "",
    "f = sys.argv[1]",
    "outdir = sys.argv[2]",
    "ad = anndata.read_h5ad(f)",
    "",
    "if ad.raw is not None:",
    "    X = ad.raw.X",
    "    genes = ad.raw.var_names",
    "else:",
    "    X = ad.X",
    "    genes = ad.var_names",
    "",
    "# anndata is cells x genes; Seurat wants genes x cells",
    "X = X.T",
    "if not sp.issparse(X):",
    "    X = sp.csr_matrix(X)",
    "",
    "sio.mmwrite(outdir + '/matrix.mtx', X)",
    "pd.Series(genes).to_csv(outdir + '/features.tsv', sep='\\t', index=False, header=False)",
    "pd.Series(ad.obs_names).to_csv(outdir + '/barcodes.tsv', sep='\\t', index=False, header=False)",
    "ad.obs.to_csv(outdir + '/obs.tsv', sep='\\t')"
  ), script)
  
  status <- system2(
    py,
    c(script, f, outdir),
    env = c(
      "PYTHONHOME=",
      paste0("PYTHONPATH=", py_path)
    )
  )
  
  if (status != 0) {
    stop("Python h5ad conversion failed for ", basename(f))
  }
  
  mat <- Matrix::readMM(file.path(outdir, "matrix.mtx"))
  
  genes <- read.delim(
    file.path(outdir, "features.tsv"),
    header = FALSE,
    stringsAsFactors = FALSE
  )[[1]]
  
  barcodes <- read.delim(
    file.path(outdir, "barcodes.tsv"),
    header = FALSE,
    stringsAsFactors = FALSE
  )[[1]]
  
  rownames(mat) <- make.unique(as.character(genes))
  colnames(mat) <- as.character(barcodes)
  mat <- as(mat, "CsparseMatrix")
  
  obj <- CreateSeuratObject(
    counts = mat,
    project = tools::file_path_sans_ext(basename(f)),
    min.cells = 3,
    min.features = 200
  )
  
  obs <- read.delim(
    file.path(outdir, "obs.tsv"),
    row.names = 1,
    check.names = FALSE
  )
  
  obs <- obs[colnames(obj), , drop = FALSE]
  obj <- AddMetaData(obj, metadata = obs)
  
  obj
}

first_meta_col <- function(md, candidates) {
  hit <- intersect(candidates, colnames(md))[1]
  if (is.na(hit)) return(NULL)
  hit
}

standardize_htan_meta <- function(o, nm) {
  md <- o@meta.data
  
  o$study <- "HTAN_islam2024"
  
  if (nm == "PCA83samples") {
    o$batch <- "HTAN_PCA83samples"
  } else {
    o$batch <- "HTAN_VUMC"
  }
  
  if (!"sample" %in% colnames(md)) {
    o$sample <- nm
  }
  
  if ("STUDYID.htan" %in% colnames(md)) {
    o$PatientID <- as.character(md[["STUDYID.htan"]])
  } else if ("HTAN Parent Data File ID" %in% colnames(md)) {
    o$PatientID <- sub("^(([^_]+_[^_]+)).*$", "\\1", md[["HTAN Parent Data File ID"]])
  } else if ("HTAN Specimen ID" %in% colnames(md)) {
    o$PatientID <- sub("^(([^_]+_[^_]+)).*$", "\\1", md[["HTAN Specimen ID"]])
  } else {
    o$PatientID <- nm
  }
  
  if ("SPECID.htan" %in% colnames(md)) {
    o$SampleID <- gsub(
      "^(HTA[0-9]+_[0-9]+_[0-9]+)_([0-9]+)$",
      "\\1\\2",
      as.character(md[["SPECID.htan"]])
    )
  } else if ("HTAN Specimen ID" %in% colnames(md)) {
    o$SampleID <- as.character(md[["HTAN Specimen ID"]])
  } else {
    o$SampleID <- NA_character_
  }
  
  o$CellType_cb <- NA_character_
  
  if ("Cell_Type" %in% colnames(md)) {
    idx <- is.na(o$CellType_cb) & !is.na(md[["Cell_Type"]])
    o$CellType_cb[idx] <- as.character(md[["Cell_Type"]][idx])
  }
  
  if ("Cell Type" %in% colnames(md)) {
    idx <- is.na(o$CellType_cb) & !is.na(md[["Cell Type"]])
    o$CellType_cb[idx] <- as.character(md[["Cell Type"]][idx])
  }
  
  type_col <- first_meta_col(md, c("HPTYPE", "Tumor_Type", "Polyp_Type", "Type"))
  if (!is.null(type_col)) {
    o$Type <- as.character(md[[type_col]])
  } else {
    o$Type <- NA_character_
  }
  
  
  cond_col <- first_meta_col(md, c(
    "TISSUE",
    "Sample_Classification",
    "sample_condition",
    "condition",
    "Status",
    "Histology"
  ))
  
  if (!is.null(cond_col)) {
    o$sample_condition <- as.character(md[[cond_col]])
  } else {
    o$sample_condition <- NA_character_
  }
  
  o
}

summarize_loaded_h5ad <- function(o, nm) {
  message("  HTAN h5ad loaded: ", nm)
  message("    genes: ", nrow(o))
  message("    cells: ", ncol(o))
  
  if ("PatientID" %in% colnames(o@meta.data)) {
    message("    patients: ", length(unique(o$PatientID)))
  } else {
    message("    patients: PatientID column not found")
  }
  
  if ("SampleID" %in% colnames(o@meta.data)) {
    message("    samples: ", length(unique(o$SampleID)))
  } else {
    message("    samples: SampleID column not found")
  }
  
  if ("CellType_cb" %in% colnames(o@meta.data)) {
    message("    CellType_cb:")
    print(table(o$CellType_cb, useNA = "ifany"))
  }
  
  if ("sample_condition" %in% colnames(o@meta.data)) {
    message("    sample_condition:")
    print(table(o$sample_condition, useNA = "ifany"))
  }
  
  if ("Type" %in% colnames(o@meta.data)) {
    message("    Type:")
    print(table(o$Type, useNA = "ifany"))
  }
  
  invisible(NULL)
}

htan_list <- list()

for (f in h5ad_files) {
  nm <- tools::file_path_sans_ext(basename(f))
  
  o <- tryCatch(load_one_h5ad(f), error = function(e) {
    message("    [warn] h5ad load failed ", nm, ": ", conditionMessage(e))
    NULL
  })
  
  if (!is.null(o)) {
    o <- standardize_htan_meta(o, nm)
    
    summarize_loaded_h5ad(o, nm)
    
    htan_list[[paste0("h5ad_", nm)]] <- o
  }
}

if (length(htan_list) == 0) stop("All HTAN files failed to load.")

objH <- if (length(htan_list) == 1) htan_list[[1]] else
  merge(htan_list[[1]], y = htan_list[-1],
        add.cell.ids = names(htan_list))

message("  HTAN before normal/NL filter:")
message(
  "    ", nrow(objH), " genes x ",
  ncol(objH), " cells"
)
message(
  "    ", length(unique(objH$SampleID)), " samples, ",
  length(unique(objH$PatientID)), " patients"
)
message("    Cells by sample_condition:")
print(table(objH$sample_condition, useNA = "ifany"))
message("    Cells by Type:")
print(table(objH$Type, useNA = "ifany"))
message("    Cells by batch:")
print(table(objH$batch, useNA = "ifany"))
message("    Cells by study:")
print(table(objH$study, useNA = "ifany"))
message("    Cells by cell type:")
print(table(objH$CellType_cb, useNA = "ifany"))


# Filter HTAN to normal/NL only
objH <- subset(
  objH,
  subset = sample_condition %in% c("Normal", "NL")
)


### message printing
count_by_file <- function(obj, group_col = "orig.ident") {
  obj@meta.data %>%
    group_by(.data[[group_col]]) %>%
    summarise(
      genes = nrow(obj),
      cells = n(),
      patients = if ("PatientID" %in% colnames(obj@meta.data)) n_distinct(PatientID) else NA_integer_,
      samples = if ("SampleID" %in% colnames(obj@meta.data)) n_distinct(SampleID) else NA_integer_,
      .groups = "drop"
    )
}
htan_by_file <- count_by_file(objH, "orig.ident")
htan_by_file


message("HTAN normal/NL cells retained: ", ncol(objH),
        " across ", length(unique(objH$PatientID)), " patient(s) ", length(unique(objH$SampleID)), " sample(s)")

message(
  "    ", nrow(objH), " genes x ",
  ncol(objH), " cells"
)
message(
  "    ", length(unique(objH$SampleID)), " samples, ",
  length(unique(objH$PatientID)), " patients"
)
message("    Cells by sample_condition:")
print(table(objH$sample_condition, useNA = "ifany"))
message("    Cells by Type:")
print(table(objH$Type, useNA = "ifany"))
message("    Cells by batch:")
print(table(objH$batch, useNA = "ifany"))
message("    Cells by study:")
print(table(objH$study, useNA = "ifany"))
message("    Cells by cell type:")
print(table(objH$CellType_cb, useNA = "ifany"))

saveRDS(objH, file.path(RESULTS_DIR, "objH_preqc.rds"))

write.csv(objH@meta.data, file.path(RESULTS_DIR, "objH_metadata_preqc.csv"), row.names = F)

rm(htan_list); gc()

