#!/usr/bin/env Rscript
# Plot S-phase fractions (healthy vs malignant, per region and per transcriptional
# cluster) using the published cell-cycle and tumour/normal annotation that ships
# with the snRNA-seq Seurat object.
#
# Input:  <ZENODO_DIR>/data/snRNA/seurat_objects/MPNST_C_updated.rds
#         (Zenodo-published, CNV/genotype-annotated snRNA object; see
#         scRNA/README.md for full provenance and column reference)
# Output: scRNA/results/sphase_summary.tsv
#         scRNA/results/sphase_fractions.png
#         scRNA/results/sphase_fractions_region_native.png
#
# This script only READS the Phase / tumor_normal / cell_type columns that ship
# with the object -- they were computed by the original authors (Seurat
# CellCycleScoring for Phase/S.Score/G2M.Score; CNV-based clustering for
# tumor_normal/cell_type). It answers "is there a tumour vs healthy difference
# in S-phase fraction, and does it vary by region/cluster?" using their calls.
#
# See 1_cellcycle2.R for a from-scratch recomputation of the cell-cycle scores
# with Seurat::CellCycleScoring() and a validation against these published calls.
#
# Requires: Seurat -- not installed in any micromamba env. Run via:
#   singularity exec /mnt/iribhm/software/singularity/single-cell.sif \
#     Rscript scRNA/1_cellcycle.R
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

OUT_TSV <- "scRNA/results/sphase_summary.tsv"
OUT_PNG <- "scRNA/results/sphase_fractions.png"

REGION_LEVELS <- c("P", "R1", "R2", "R3", "R4", "R5")
MIN_N_ALPHA   <- 15   # groups with fewer cells are shown at 40% opacity (matches
                      # the convention in scDNA/4_sprinter/persample/plot_sphase_persample.R)

# Categorical pair (blue/red, dataviz-skill slots 1 and 8 -- max separation,
# colourblind-safe): Healthy vs Malignant is an identity split, not a status.
CLASS_COLOURS <- c(Healthy = "#2a78d6", Malignant = "#e34948")

dir.create(dirname(OUT_TSV), recursive = TRUE, showWarnings = FALSE)

# ── Load ───────────────────────────────────────────────────────────────────

cat("Loading Seurat object (~2 GB, may take a minute)...\n")
so <- readRDS(SEURAT_RDS)
md <- so@meta.data

stopifnot(all(c("Phase", "tumor_normal", "cell_type") %in% colnames(md)))
cat("Cells:", nrow(md), "\n")

# ── Derive region / group / class ────────────────────────────────────────────
#
# Barcodes are "<region>_<10x-barcode>", e.g. "P_AAACCCAAGCGACTTT". Normal cells
# (tumor_normal == "Normal") are collapsed to a single "healthy" group per region
# -- cell_type would otherwise split them into Macrophage/T_cell/Endothelial/
# Skeletal_Muscle, which is finer-grained than the tumour/healthy question asked
# here. Malignant cells keep their published cell_type label, which mostly
# encodes a region-specific subclone (e.g. "Malignant_R1_1") but a few clusters
# ("Malignant_G2MS_R1R5", "Malignant_Ribo", ...) span multiple regions or reflect
# a shared transcriptional state rather than a single region of origin -- these
# are shown, faithfully, in every region panel where they actually have cells.

md$region <- sub("_.*", "", rownames(md))
stopifnot(all(md$region %in% REGION_LEVELS))
md$region <- factor(md$region, levels = REGION_LEVELS)

md$class <- ifelse(md$tumor_normal == "Normal", "Healthy", "Malignant")
md$group <- ifelse(md$class == "Healthy", "healthy", as.character(md$cell_type))
md$is_s  <- as.integer(md$Phase == "S")

# ── Aggregate: n_cells / n_sphase / sphase_fraction per region x group ──────

n_cells_df  <- aggregate(is_s ~ region + group + class, data = md, FUN = length)
names(n_cells_df)[names(n_cells_df) == "is_s"] <- "n_cells"

n_sphase_df <- aggregate(is_s ~ region + group + class, data = md, FUN = sum)
names(n_sphase_df)[names(n_sphase_df) == "is_s"] <- "n_sphase"

summary_df <- merge(n_cells_df, n_sphase_df, by = c("region", "group", "class"))
summary_df$sphase_fraction <- summary_df$n_sphase / summary_df$n_cells
summary_df <- summary_df[order(summary_df$region, -summary_df$n_cells), ]

write.table(summary_df, OUT_TSV, sep = "\t", quote = FALSE, row.names = FALSE)
cat("Saved:\n ", OUT_TSV, "\n")

# ── Order groups for plotting: healthy first, then by descending total n ────

totals <- aggregate(n_cells ~ group, data = summary_df, sum)
totals <- totals[order(-totals$n_cells), ]
group_order <- c("healthy", setdiff(totals$group, "healthy"))

summary_df$group <- factor(summary_df$group, levels = group_order)
summary_df$class <- factor(summary_df$class, levels = c("Healthy", "Malignant"))
summary_df$alpha_val <- ifelse(summary_df$n_cells >= MIN_N_ALPHA, 1.0, 0.4)

# ── Plot ─────────────────────────────────────────────────────────────────────

p <- ggplot(summary_df, aes(x = group, y = sphase_fraction,
                             fill = class, alpha = alpha_val)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = n_cells),
            angle = 90, hjust = -0.1, vjust = 0.5,
            size = 2, colour = "grey30", alpha = 1) +
  scale_fill_manual(values = CLASS_COLOURS, name = NULL) +
  scale_alpha_identity(guide = "none") +
  scale_y_continuous(labels = percent_format(accuracy = 1),
                     limits = c(0, NA),
                     expand = expansion(mult = c(0, 0.15))) +
  facet_wrap(~ region, nrow = 1, scales = "free_x") +
  labs(x = NULL, y = "S-phase fraction",
       caption = paste(
         "Bar labels = n cells. Bars with n_cells < 15 shown at 40% opacity.",
         "Clusters that span multiple regions (e.g. Malignant_G2MS_R1R5) appear",
         "in every region panel where they have cells.",
         sep = "\n")) +
  theme_bw(base_size = 11) +
  theme(axis.text.x        = element_text(angle = 45, hjust = 1, size = 7),
        strip.background   = element_rect(fill = "grey92"),
        strip.text         = element_text(face = "bold"),
        panel.grid.major.x = element_blank(),
        legend.position    = "bottom",
        plot.caption       = element_text(size = 8, colour = "grey50", hjust = 0))

ggsave(OUT_PNG, p, width = 300, height = 130, units = "mm", dpi = 300)
cat("Saved:\n ", OUT_PNG, "\n")

# ── Plot (region-native clusters only) ──────────────────────────────────────
#
# Same data, but each region panel keeps only "healthy" plus clusters whose
# published cell_type label actually names that region (e.g. the "R1" panel
# keeps Malignant_R1_1/R1_2 but drops the 162 Malignant_R5 cells that are
# physically located in R1). Clusters with no single native region at all
# (Malignant_Ribo, Malignant_G2MS_*) are dropped everywhere. This mirrors the
# region-native filter in scDNA/4_sprinter/persample/plot_sphase_persample.R
# ("Keep only region-native clones and healthy; exclude cross-region and
# unassigned").

native_region <- sub("_.*", "", sub("^Malignant_", "", as.character(summary_df$group)))
is_native <- summary_df$group == "healthy" | native_region == as.character(summary_df$region)

fully_dropped <- setdiff(unique(as.character(summary_df$group[!is_native])),
                          unique(as.character(summary_df$group[is_native])))
leakage_rows  <- summary_df[!is_native & !(as.character(summary_df$group) %in% fully_dropped),
                             c("region", "group", "n_cells")]

cat("\nRegion-native plot -- clusters with no single native region (excluded everywhere):\n ",
    paste(fully_dropped, collapse = ", "), "\n")
if (nrow(leakage_rows) > 0) {
  cat("Region-native plot -- cross-region rows dropped (cluster kept only in its own region):\n")
  print(leakage_rows, row.names = FALSE)
}

summary_df_native <- summary_df[is_native, ]
summary_df_native$group <- factor(as.character(summary_df_native$group), levels = group_order)

OUT_PNG_NATIVE <- "scRNA/results/sphase_fractions_region_native.png"

p_native <- ggplot(summary_df_native, aes(x = group, y = sphase_fraction,
                                           fill = class, alpha = alpha_val)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = n_cells),
            angle = 90, hjust = -0.1, vjust = 0.5,
            size = 2, colour = "grey30", alpha = 1) +
  scale_fill_manual(values = CLASS_COLOURS, name = NULL) +
  scale_alpha_identity(guide = "none") +
  scale_y_continuous(labels = percent_format(accuracy = 1),
                     limits = c(0, NA),
                     expand = expansion(mult = c(0, 0.15))) +
  facet_wrap(~ region, nrow = 1, scales = "free_x") +
  labs(x = NULL, y = "S-phase fraction",
       caption = paste(
         "Bar labels = n cells. Bars with n_cells < 15 shown at 40% opacity.",
         "Only region-native clusters kept per panel (e.g. P keeps only Malignant_P);",
         "cross-region rows and non-region-specific clusters (Ribo, G2MS_*) excluded.",
         sep = "\n")) +
  theme_bw(base_size = 11) +
  theme(axis.text.x        = element_text(angle = 45, hjust = 1, size = 8),
        strip.background   = element_rect(fill = "grey92"),
        strip.text         = element_text(face = "bold"),
        panel.grid.major.x = element_blank(),
        legend.position    = "bottom",
        plot.caption       = element_text(size = 8, colour = "grey50", hjust = 0))

ggsave(OUT_PNG_NATIVE, p_native, width = 300, height = 130, units = "mm", dpi = 300)
cat("Saved:\n ", OUT_PNG_NATIVE, "\n")

# ── Headline numbers ─────────────────────────────────────────────────────────

overall <- aggregate(cbind(n_cells, n_sphase) ~ class, data = summary_df, sum)
overall$sphase_fraction <- overall$n_sphase / overall$n_cells
cat("\n=== Overall S-phase fraction, Healthy vs Malignant (all regions pooled) ===\n")
print(overall)
