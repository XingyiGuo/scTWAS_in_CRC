## 05_sctwas_multisource_DEsupport_diverging.R --------------------------------
## Panel C: how many of the 347 scTWAS genes are supported, and in which
## direction, by each external differential-expression resource, broken down by
## the cell-type model that nominated the gene.
##
## Three platforms, one contrast: normal epithelium against colorectal
## carcinoma, left to right
##   CosMx_NOR_CORE   spatial single cell, RNA (carcinoma core)
##   TCGA             bulk RNA
##   CPTAC            bulk protein
##
## Excluded by decision: the CosMx adenoma -> core and normal -> adenoma
## transitions (2 and 5 supported genes of 347) and the two Visium modules
## (a different contrast, and one-sided published tables).
##
## Input
##   table/DE_key_results/c_scTWAS_347genes_multisource_DE_support.csv
##     one row per gene (347), with a *_DEsupport call per resource assigned
##     upstream in 10_sctwas_gene_level_direction_and_support.R under the house
##     rule padj < 0.05 AND |log2FC| > 0.5 AND sign(log2FC) matching the gene's
##     representative scTWAS Z.
##
## Counting
##   One gene, one row. Each gene is placed in the row of its representative
##   cell-type model (rep_Dataset x rep_CellType) - the model with the smallest
##   scTWAS P value, which is also the model whose Z fixed the gene's direction.
##   A gene nominated in several models therefore contributes to exactly one
##   row, so bars are counts of distinct genes and sum to the gene-level totals
##   reported elsewhere (CosMx normal->core 58, Visium M1 16, Visium M3 15,
##   TCGA 79, CPTAC 15). The consequence to keep in mind is the converse: a row
##   understates how many models implicate that cell type.
##
## Shared x scale
##   All five panels share one x range rather than each being scaled to its own
##   maximum. The resources differ several-fold in how much support they carry
##   (largest single bar: TCGA 10 genes, CPTAC 2), and per panel scaling would
##   draw a two-gene bar at the same length as a ten-gene bar. The thin panels
##   are meant to read as thin.
##
## Output
##   figure/Fig_scTWAS_multisource_DEsupport_diverging.{png,pdf}
##   figure/Fig_scTWAS_multisource_DEsupport_diverging_sourcedata.csv
##   table/DE_key_results/scTWAS_multisource_DEsupport_by_model.csv
## -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

PROJ <- "W:/VUMC_projects/scTWAS/VUMC_HTAN_scTWAS_genes_DE_support"
if (!dir.exists(file.path(PROJ, "table"))) PROJ <- "."
DIR_KEY <- file.path(PROJ, "table", "DE_key_results")
DIR_FIG <- file.path(PROJ, "figure")
dir.create(DIR_FIG, recursive = TRUE, showWarnings = FALSE)

GENE_CSV <- file.path(DIR_KEY,
  "c_scTWAS_347genes_multisource_DE_support.csv")

## Resource order and strip labels. Spatial first, then bulk.
##
## All three resources now measure the SAME biological contrast - normal
## epithelium against colorectal carcinoma - and differ only in platform, so
## the panels are directly comparable and the x axis can name the direction in
## carcinoma terms rather than as a generic "later state".
##
## Excluded by decision: the CosMx adenoma -> core and normal -> adenoma
## transitions (2 and 5 supported genes of 347), and the two Visium modules
## (a different contrast, and their published tables are up-genes only).
##
## one_sided marks a resource whose deposited table contains only up-regulated
## genes. None of the three qualify now, so no panel is shaded; the mechanism
## is kept so a one-sided resource cannot be added silently.
RES <- data.table(
  resource  = c("CosMx_NOR_CORE", "TCGA", "CPTAC"),
  short     = c("CosMx", "TCGA", "CPTAC"),
  label = c("CosMx\nnormal \u2192 colorectal carcinoma\n(spatial single cell, RNA)",
            "TCGA\ncolorectal carcinoma vs normal\n(bulk RNA)",
            "CPTAC\ncolorectal carcinoma vs normal\n(bulk protein)"),
  one_sided = c(FALSE, FALSE, FALSE))

## One colour pair per platform, so a panel is identifiable by colour alone.
##
## Palette constraints applied here:
##   - no red / green pairing anywhere, which journals in this family will not
##     accept and which red-green colour-vision deficiency cannot resolve;
##   - low saturation throughout (roughly 35-45%), the register these journals
##     print in, rather than the fully saturated ColorBrewer endpoints used in
##     the earlier draft of this figure;
##   - warm for up and cool for down in every pair, so direction reads the same
##     way in all three panels and only the hue distinguishes the platform;
##   - the three hue pairs are well separated: terracotta / steel blue,
##     ochre / violet, dusty rose / teal.
COLS <- c("CosMx up"  = "#CC8963", "CosMx down" = "#5C8CB8",
          "TCGA up"   = "#BFA25A", "TCGA down"  = "#8874A3",
          "CPTAC up"  = "#C6879B", "CPTAC down" = "#5FA3A0")

## cell-type row order, following the reference figure
CT_ORDER    <- c("ABS", "GOB", "TAC", "STM", "TUF", "EE", "CT", "ASC", "SSC")
BLOCK_LABEL <- c(NL = "NL models", ADSER = "Polyps models")

## 1. input -------------------------------------------------------------------
cc <- fread(GENE_CSV)
supcols <- paste0(RES$resource, "_DEsupport")
stopifnot(nrow(cc) == 347L, uniqueN(cc$gene) == 347L,
          all(c("rep_Dataset", "rep_CellType") %in% names(cc)),
          all(supcols %in% names(cc)))
stopifnot(all(unlist(cc[, ..supcols]) %chin%
                c("Up", "Down", "DE_discordant", "NS", "not_tested")))
stopifnot(setequal(unique(cc$rep_Dataset), c("NL", "ADSER")),
          all(cc$rep_CellType %chin% CT_ORDER))

## 2. counts: one gene in one row -------------------------------------------
m <- melt(cc[, c("gene", "gene_name", "rep_Dataset", "rep_CellType",
                 supcols), with = FALSE],
          id.vars = c("gene", "gene_name", "rep_Dataset", "rep_CellType"),
          variable.name = "resource", value.name = "call")
m[, resource := sub("_DEsupport$", "", as.character(resource))]

## complete grid of model x resource x direction, so a zero is a drawn zero
models <- unique(cc[, .(rep_Dataset, rep_CellType)])
grid <- CJ(idx = seq_len(nrow(models)), resource = RES$resource,
           direction = c("Up", "Down"), unique = TRUE)
grid <- cbind(models[grid$idx], grid[, .(resource, direction)])

obs <- m[call %chin% c("Up", "Down"),
         .(n_genes = .N), by = .(rep_Dataset, rep_CellType, resource,
                                 direction = call)]
k <- merge(grid, obs, by = c("rep_Dataset", "rep_CellType", "resource",
                             "direction"), all.x = TRUE)
k[is.na(n_genes), n_genes := 0L]

## denominator: genes whose representative model is this one
den <- cc[, .(model_genes = .N), by = .(rep_Dataset, rep_CellType)]
k <- merge(k, den, by = c("rep_Dataset", "rep_CellType"))
k <- merge(k, RES, by = "resource", sort = FALSE)

## the per-resource totals must reproduce the gene-level support counts
chk <- k[, .(total = sum(n_genes)), by = resource]
ref <- m[call %chin% c("Up", "Down"), .(ref = .N), by = resource]
stopifnot(nrow(merge(chk, ref, by = "resource")[total != ref]) == 0L)

k[, block := factor(BLOCK_LABEL[rep_Dataset], levels = BLOCK_LABEL)]
k[, CellType := factor(rep_CellType, levels = rev(CT_ORDER))]
k[, label := factor(label, levels = RES$label)]
k[, fill_key := factor(sprintf("%s %s", short, tolower(direction)),
                       levels = names(COLS))]
stopifnot(!anyNA(k$fill_key))
## down to the left of the baseline, up to the right
k[, x_signed := fifelse(direction == "Down", -n_genes, n_genes)]
setorder(k, resource, block, CellType, direction)

fwrite(k[, .(block, rep_Dataset, rep_CellType, model_genes, resource,
             resource_label = gsub("\n", " ", label),
             direction, n_genes)],
       file.path(DIR_KEY, "scTWAS_multisource_DEsupport_by_model.csv"))

## 3. figure -----------------------------------------------------------------
XMAX <- max(k$n_genes)

## how many of the 347 genes each resource tests at all, read from the data so
## the caption cannot drift from the table
tested <- m[, .(n_tested = sum(call != "not_tested")), by = resource]
tested <- merge(RES[, .(resource, label)], tested, by = "resource", sort = FALSE)
TESTED_TXT <- paste(sprintf("%s %d", sub("\n.*$", "", tested$label),
                            tested$n_tested), collapse = ", ")
## rows where neither direction is supported get a baseline stub, so an empty
## row reads as "tested, none supported" rather than as a missing row
## Rows tested with nothing supported. These carry a faint tick at zero, the
## convention already used in the companion TCGA / CPTAC figure, so an empty
## slot is not read as a missing category. Drawn as a point with shape 124
## rather than as a tile or a rect: a zero-row or constant-aes background layer
## silently took over the fill mapping in an earlier version of this figure and
## rendered every panel in the first palette.
stub <- k[, .(total = sum(n_genes)), by = .(block, CellType, label)][total == 0L]
stopifnot(!any(RES$one_sided))   # no panel needs a structurally-empty half

p <- ggplot(k, aes(x = x_signed, y = CellType, fill = fill_key)) +
  geom_vline(xintercept = 0, colour = "grey25", linewidth = 0.35) +
  geom_col(width = 0.68) +
  geom_point(data = stub, aes(x = 0, y = CellType), inherit.aes = FALSE,
             shape = 124, size = 1.6, colour = "grey60") +
  geom_text(data = k[n_genes > 0],
            aes(x = x_signed, y = CellType, label = n_genes, colour = fill_key,
                hjust = fifelse(direction == "Down", 1.35, -0.35)),
            inherit.aes = FALSE, size = 1.9, show.legend = FALSE) +
  facet_grid(block ~ label, scales = "free_y", space = "free_y") +
  scale_fill_manual(values = COLS, name = NULL, drop = FALSE) +
  scale_colour_manual(values = COLS, guide = "none") +
  scale_x_continuous(limits = c(-XMAX * 1.28, XMAX * 1.28),
                     breaks = seq(-10, 10, 5), labels = function(v) abs(v),
                     expand = expansion(mult = 0)) +
  labs(x = paste0("supported genes:   \u2190 down in carcinoma   |   ",
                  "up in carcinoma \u2192"),
       y = "cell-type model that nominated the gene",
       title = paste0("One normal-to-carcinoma contrast across three ",
                      "platforms: 79 genes supported in bulk RNA, 58 in ",
                      "spatial RNA, 15 in protein"),
       subtitle = paste0("347 scTWAS genes, each counted once in the row of ",
                         "its representative cell-type model"),
       caption = paste0(
         "Supported = padj < 0.05, |log2FC| > 0.5 and sign(log2FC) matching the gene's representative scTWAS Z, so up and down are both supported and\n",
         "record the direction of the expression change in carcinoma. A gene nominated in several cell-type models is placed only in the model with its\n",
         "smallest scTWAS P value, which means bars are counts of distinct genes but a row understates how many models implicate that cell type. Each\n",
         "platform has its own colour pair, warm for up and cool for down, and all three panels share one x scale. A faint grey tick at zero marks a\n",
         "model x platform combination that was tested with no supported gene, distinguishing it from a row that is absent. The platforms differ\n",
         "in how many of the 347 genes they measure at all: ", TESTED_TXT, ", which is why the protein counts are the\n",
         "smallest. The CosMx contrast is against the carcinoma core specifically, as distinct from the invasive front, which is not shown here.\n",
         "ASC and SSC exist only among the polyp models.")) +
  theme_classic(base_size = 8) +
  theme(plot.title    = element_text(size = 8.5, hjust = 0),
        plot.subtitle = element_text(size = 7, hjust = 0, colour = "grey30"),
        plot.caption  = element_text(size = 5.5, hjust = 0, colour = "grey30"),
        axis.text.x   = element_text(size = 6, colour = "black"),
        axis.text.y   = element_text(size = 6.5, colour = "black"),
        axis.title    = element_text(size = 7.5),
        axis.line.y   = element_blank(),
        axis.ticks.y  = element_blank(),
        strip.background = element_blank(),
        strip.text.x  = element_text(size = 7, lineheight = 1.1),
        strip.text.y  = element_text(size = 7.5, angle = -90),
        panel.grid.major.x = element_line(colour = "grey93", linewidth = 0.25),
        panel.spacing.x = unit(2.4, "mm"),
        panel.spacing.y = unit(2.6, "mm"),
        legend.position = "top",
        legend.text = element_text(size = 6.5),
        legend.key.size = unit(3, "mm"))

FIG_STEM <- file.path(DIR_FIG, "Fig_scTWAS_multisource_DEsupport_diverging")
fwrite(k[, .(block, rep_Dataset, rep_CellType, model_genes, resource,
             resource_label = gsub("\n", " ", label), direction,
             n_genes, x_plotted = x_signed, fill_colour = COLS[fill_key])],
       paste0(FIG_STEM, "_sourcedata.csv"))
ggsave(paste0(FIG_STEM, ".png"), plot = p, width = 7.6, height = 5.6,
       units = "in", dpi = 300)
ggsave(paste0(FIG_STEM, ".pdf"), plot = p, width = 7.6, height = 5.6,
       units = "in", useDingbats = FALSE)

## 4. console summary --------------------------------------------------------
cat("\n== supported genes per resource and direction ==\n")
print(dcast(k[, .(n = sum(n_genes)), by = .(resource, direction)],
            resource ~ direction, value.var = "n")[
              order(match(resource, RES$resource))])
cat("\n== models with any support, of 16 ==\n")
print(k[, .(models_with_support = uniqueN(paste(rep_Dataset, rep_CellType)[n_genes > 0]),
            largest_bar = max(n_genes)), by = resource][
              order(match(resource, RES$resource))])
cat("\n== largest bars ==\n")
print(head(k[order(-n_genes), .(resource, block, rep_CellType, direction,
                                n_genes)], 10))
cat("\nwrote ", FIG_STEM, ".{png,pdf}\n", sep = "")
