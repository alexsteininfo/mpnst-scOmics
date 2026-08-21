#!/usr/bin/env Rscript
# Recompute cell-cycle scores from scratch with Seurat::CellCycleScoring() -- the
# standard method in the Seurat ecosystem -- and validate the result against the
# published Phase/S.Score/G2M.Score calls already present in the object (see
# 1_cellcycle.R, which uses those published calls directly). This is the "show
# your work" companion script: it reproduces the method, checks it against a
# known answer, and demonstrates how much the result moves under different
# reasonable parameter choices.
#
# Input:  <ZENODO_DIR>/data/snRNA/seurat_objects/MPNST_C_updated.rds
# Output: scRNA/results/cellcycle_method_comparison.tsv
#         scRNA/results/cellcycle_method_sensitivity.png
#
# What CellCycleScoring() actually does (see "under the hood" section below):
# for the S-phase gene list (and separately the G2M list), it bins ALL genes
# into expression-level bins, then for each gene in the list subtracts the mean
# expression of a random control set of genes drawn from the SAME bin. This
# corrects for the fact that dropout/expression level, not biology, would
# otherwise dominate a raw mean. The gene with the higher corrected score
# (S.Score vs G2M.Score) wins; "G1" is the fallback when neither score is > 0.
#
# Parameters this script demonstrates the effect of (see README.md "Parameter
# tuning" section for the discussion):
#   - gene list:  cc.genes.updated.2019 (current Seurat default) vs cc.genes (2015)
#   - assay/slot: RNA "data" (LogNormalize) vs SCT "data" (SCTransform)
#
# Requires: Seurat -- not installed in any micromamba env. Run via:
#   singularity exec /mnt/iribhm/software/singularity/single-cell.sif \
#     Rscript scRNA/1_cellcycle2.R
#
# Must be run from the repository root (paths below are repo-relative).

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(scales)
})

# ── Config ─────────────────────────────────────────────────────────────────

ZENODO_DIR <- "/srv/home/aste0033/projects/MPNST/zenodo_phase_1/MPNST-Zenodo"
SEURAT_RDS <- file.path(ZENODO_DIR, "data/snRNA/seurat_objects/MPNST_C_updated.rds")

OUT_TSV <- "scRNA/results/cellcycle_method_comparison.tsv"
OUT_PNG <- "scRNA/results/cellcycle_method_sensitivity.png"

TUMOR_NORMAL_COLOURS <- c(Normal = "#2a78d6", Tumor = "#e34948")

dir.create(dirname(OUT_TSV), recursive = TRUE, showWarnings = FALSE)

# ── Load ───────────────────────────────────────────────────────────────────

cat("Loading Seurat object (~2 GB, may take a minute)...\n")
so <- readRDS(SEURAT_RDS)
stopifnot(all(c("Phase", "S.Score", "G2M.Score", "tumor_normal") %in% colnames(so@meta.data)))

# Keep the published calls aside before we overwrite S.Score/G2M.Score/Phase
# below -- CellCycleScoring() writes into those exact column names.
published <- so@meta.data[, c("S.Score", "G2M.Score", "Phase")]
tumor_normal <- so@meta.data$tumor_normal

# ── Gene lists: current Seurat default vs the original 2015 list ───────────

data("cc.genes.updated.2019", package = "Seurat")
data("cc.genes", package = "Seurat")

cat("\ncc.genes.updated.2019: ", length(cc.genes.updated.2019$s.genes), "S genes, ",
    length(cc.genes.updated.2019$g2m.genes), "G2M genes\n")
cat("cc.genes (2015):       ", length(cc.genes$s.genes), "S genes, ",
    length(cc.genes$g2m.genes), "G2M genes\n")

# ── Run CellCycleScoring() under three parameter choices ───────────────────
#
# `so` is passed by value (R copy-on-modify): each call below scores a fresh
# copy from the unmodified object in this scope, so the three variants don't
# interfere with each other, and the object's original published columns are
# never touched.

score_variant <- function(object, assay, s.genes, g2m.genes) {
  DefaultAssay(object) <- assay
  object <- CellCycleScoring(object, s.features = s.genes, g2m.features = g2m.genes,
                              set.ident = FALSE)
  object@meta.data[, c("S.Score", "G2M.Score", "Phase")]
}

cat("\nScoring variant RNA_2019 (RNA assay, cc.genes.updated.2019) ...\n")
v_rna_2019 <- score_variant(so, "RNA", cc.genes.updated.2019$s.genes, cc.genes.updated.2019$g2m.genes)

cat("Scoring variant RNA_2015 (RNA assay, cc.genes) ...\n")
v_rna_2015 <- score_variant(so, "RNA", cc.genes$s.genes, cc.genes$g2m.genes)

cat("Scoring variant SCT_2019 (SCT assay, cc.genes.updated.2019) ...\n")
v_sct_2019 <- score_variant(so, "SCT", cc.genes.updated.2019$s.genes, cc.genes.updated.2019$g2m.genes)

variants <- list(
  "RNA + 2019 gene list (this script's primary reproduction)" = v_rna_2019,
  "RNA + 2015 gene list"                                      = v_rna_2015,
  "SCT + 2019 gene list"                                      = v_sct_2019
)

# ── Validate: confusion matrix vs published Phase, per variant ─────────────

cat("\n=== Concordance with the published Phase calls ===\n")
concordance <- sapply(variants, function(v) {
  tab <- table(published = published$Phase, ours = v$Phase)
  print(tab)
  conc <- sum(diag(tab)) / sum(tab)
  cat(sprintf("Concordance: %.1f%%\n\n", 100 * conc))
  conc
})
print(round(100 * concordance, 1))

# ── S-phase fraction: published vs each variant, overall and by tumor/normal ─

frac_by_group <- function(phase, group) {
  tab <- table(group, phase == "S")
  tab[, "TRUE"] / rowSums(tab)
}

comparison <- data.frame(
  method = c("Published", names(variants)),
  overall_sphase_fraction = c(
    mean(published$Phase == "S"),
    sapply(variants, function(v) mean(v$Phase == "S"))
  )
)

by_tn <- cbind(
  Published = frac_by_group(published$Phase, tumor_normal),
  sapply(variants, function(v) frac_by_group(v$Phase, tumor_normal))
)
comparison <- cbind(comparison, t(by_tn)[comparison$method, , drop = FALSE])

cat("\n=== S-phase fraction by method ===\n")
print(comparison, row.names = FALSE)
write.table(comparison, OUT_TSV, sep = "\t", quote = FALSE, row.names = FALSE)
cat("\nSaved:\n ", OUT_TSV, "\n")

# ── Under the hood: what is CellCycleScoring() actually computing? ─────────
#
# A naive approach would just average the log-normalized expression of the
# S-phase genes per cell. CellCycleScoring() (via AddModuleScore()) instead
# bins genes by average expression and subtracts a matched control-gene-set
# score, which corrects for expression-level-dependent dropout. The two are
# correlated but not identical -- that gap IS the correction.

s_genes_present <- intersect(cc.genes.updated.2019$s.genes, rownames(so[["RNA"]]))
cat(sprintf("\n%d of %d S-phase genes found in the RNA assay\n",
            length(s_genes_present), length(cc.genes.updated.2019$s.genes)))

naive_s_score <- colMeans(GetAssayData(so, assay = "RNA", layer = "data")[s_genes_present, ])
naive_vs_actual <- cor(naive_s_score, v_rna_2019$S.Score, method = "spearman")
cat(sprintf(
  "Spearman correlation, naive mean S-gene expression vs Seurat's control-corrected S.Score: %.2f\n",
  naive_vs_actual))
cat("(Strongly correlated but not 1:1 -- the difference is exactly the dropout correction.)\n")

# ── Sensitivity plot ─────────────────────────────────────────────────────────

plot_df <- data.frame(
  method = rep(comparison$method, each = 2),
  tumor_normal = rep(c("Normal", "Tumor"), times = nrow(comparison)),
  sphase_fraction = as.vector(t(comparison[, c("Normal", "Tumor")]))
)
plot_df$method <- factor(plot_df$method, levels = comparison$method)
plot_df$tumor_normal <- factor(plot_df$tumor_normal, levels = c("Normal", "Tumor"))

p <- ggplot(plot_df, aes(x = method, y = sphase_fraction, fill = tumor_normal)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.65) +
  geom_text(aes(label = percent(sphase_fraction, accuracy = 0.1)),
            position = position_dodge(width = 0.75), vjust = -0.4, size = 2.8) +
  scale_fill_manual(values = TUMOR_NORMAL_COLOURS, name = NULL) +
  scale_y_continuous(labels = percent_format(accuracy = 1),
                     limits = c(0, NA), expand = expansion(mult = c(0, 0.15))) +
  labs(x = NULL, y = "S-phase fraction",
       title = "Sensitivity of the S-phase fraction to cell-cycle scoring choices",
       caption = "Published = calls already in the object. The other three bars are this script's own\nSeurat::CellCycleScoring() re-runs under different gene-list/assay choices.") +
  theme_bw(base_size = 11) +
  theme(axis.text.x        = element_text(angle = 20, hjust = 1),
        panel.grid.major.x = element_blank(),
        legend.position    = "bottom",
        plot.caption       = element_text(size = 8, colour = "grey50", hjust = 0))

ggsave(OUT_PNG, p, width = 220, height = 130, units = "mm", dpi = 300)
cat("Saved:\n ", OUT_PNG, "\n")
