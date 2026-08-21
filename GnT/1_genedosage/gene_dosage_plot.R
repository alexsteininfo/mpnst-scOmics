#!/usr/bin/env Rscript
# Gene-dosage effect: does DNA copy number predict RNA expression in GnT cells?
#
# For each gene, compares two DNA copy-number states (n_high > n_low) observed
# across the 130 GnT cells with matched DNA+RNA. Under a strict linear-dosage
# model (each extra copy contributes equally, no compensation), the expected
# mean expression at n_high is just the mean actually observed at n_low,
# rescaled by (n_high / n_low). Plotting that expected value against what's
# actually observed at n_high lets genes that escape (or exceed) linear
# dosage stand out below/above the y=x line.
#
# Reproduces the GnT gene-dosage panel from the manuscript (published as
# Figure 2E; the mpnst_phase_1 code repo's own internal draft numbering calls
# the same analysis figure_3E.R).
#
# Input:
#   <GnT staging dir>/GnT_scRNA_MPNST_C.rds       Seurat object, 300 GnT RNA cells
#   <GnT staging dir>/MPNST_GnT_cell_seg_exp.rds  list of 130 per-cell (gene x CN, SCT_count) tables
# Output:
#   GnT/1_genedosage/results/gene_dosage_expected_vs_observed.png
#   GnT/1_genedosage/results/gene_dosage_pairs.tsv
#
# Requires: Seurat -- not installed in any micromamba env. Run via:
#   singularity exec /mnt/iribhm/software/singularity/single-cell.sif \
#     Rscript GnT/1_genedosage/gene_dosage_plot.R
#
# Must be run from the repository root (paths below are repo-relative).

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(scales)
})

# ── Config ────────────────────────────────────────────────────────────────

GNT_DIR          <- "/srv/home/aste0033/projects/MPNST/Haixi/GnT/zenodo"
CELL_SEG_EXP_RDS <- file.path(GNT_DIR, "MPNST_GnT_cell_seg_exp.rds")
SEURAT_RDS       <- file.path(GNT_DIR, "GnT_scRNA_MPNST_C.rds")

OUT_DIR    <- "GnT/1_genedosage/results"
MIN_CELLS  <- 10  # a (gene, CN) group must have >= this many cells to be used
EXP_CUTOFF <- 10  # a gene must have median SCT count >= this (across all 300
                   # GnT RNA cells, not just the 130 DNA-matched ones -- this
                   # matches the manuscript's own gene_dosage code) to be "highly expressed"
N_HIGHLIGHT <- 25 # number of largest expected-vs-observed discrepancies to flag

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ── 1. Match DNA cells (cell_seg_exp) to their RNA-side EGAF accession ──────
# GnT DNA and RNA libraries for the same cell were sequenced back-to-back and
# deposited at EGA with adjacent accessions: RNA accession = DNA accession - 1
# (see GnT/README.md). This is the only place that adjacency is used here --
# MPNST_GnT_cell_seg_exp.rds already carries per-cell SCT_count merged in via
# that mapping upstream.

cell_seg_exp <- readRDS(CELL_SEG_EXP_RDS)
gnt_scrna <- readRDS(SEURAT_RDS)

dna_barcodes <- unique(unlist(lapply(cell_seg_exp, function(df) df$barcode)))
rna_barcodes <- paste0("EGAF0000", as.numeric(gsub("EGAF", "", dna_barcodes)) - 1)
gnt_scrna_sub <- subset(gnt_scrna, cells = colnames(gnt_scrna)[colnames(gnt_scrna) %in% rna_barcodes])

# ── 2. Highly-expressed genes: median SCT count >= EXP_CUTOFF across all ───
# 300 GnT RNA cells (not restricted to the 130 DNA-matched ones)

sct_count_mtx <- as.matrix(gnt_scrna@assays[["SCT"]]@counts)
high_exp_genes <- rownames(sct_count_mtx)[apply(sct_count_mtx, 1, median) >= EXP_CUTOFF]

# ── 3. Per-(gene, CN) observed mean/median expression, >= MIN_CELLS cells ──

all_seg_exp <- do.call(rbind, cell_seg_exp) %>%
  filter(gene_name %in% high_exp_genes, CN > 0)

mean_exp_by_cn <- all_seg_exp %>%
  group_by(gene_name, CN) %>%
  filter(n() >= MIN_CELLS) %>%
  summarize(n_cells = n(),
            median_SCT_count = median(SCT_count, na.rm = TRUE),
            mean_SCT_count = mean(SCT_count, na.rm = TRUE),
            .groups = "drop")

# ── 4. Expected-vs-observed dosage pairs: every (CN_high, CN_low) pair of ──
# the same gene, expected_high = mean_observed_low * (CN_high / CN_low)

dosage_pairs <- mean_exp_by_cn %>%
  inner_join(mean_exp_by_cn, by = "gene_name", suffix = c("_high", "_low"),
             relationship = "many-to-many") %>%
  filter(CN_high > CN_low) %>%
  mutate(mean_expected_high = mean_SCT_count_low * (CN_high / CN_low),
         log10_abs_diff = log10(abs(mean_expected_high - mean_SCT_count_high))) %>%
  arrange(desc(log10_abs_diff)) %>%
  mutate(highlight = row_number() <= N_HIGHLIGHT)

write.table(dosage_pairs, file.path(OUT_DIR, "gene_dosage_pairs.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# ── 5. Plot: mean expected vs mean observed normalized UMI count (SCT) ─────

png(file.path(OUT_DIR, "gene_dosage_expected_vs_observed.png"),
    width = 1500, height = 1500, res = 200)
options(scipen = 999)
print(
  ggplot(dosage_pairs, aes(x = mean_expected_high, y = mean_SCT_count_high, colour = highlight)) +
    geom_point() +
    geom_abline(intercept = 0, slope = 1, colour = "red") +
    scale_x_log10(limits = c(1, 100000), labels = scales::comma) +
    scale_y_log10(limits = c(1, 100000), labels = scales::comma) +
    xlab("Mean expected normalized UMI count (SCT)") +
    ylab("Mean observed normalized UMI count (SCT)") +
    scale_color_manual(values = c(`TRUE` = "black", `FALSE` = "black")) +
    theme_classic(base_size = 16) +
    theme(aspect.ratio = 1, legend.position = "none")
)
dev.off()

cat(sprintf("Wrote %d dosage pairs (%d genes) to %s/\n",
            nrow(dosage_pairs), length(unique(dosage_pairs$gene_name)), OUT_DIR))
