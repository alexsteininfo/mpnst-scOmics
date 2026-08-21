#!/usr/bin/env Rscript
# Compare S-phase fractions across all three technologies (10X, DLP+, scRNA)
# side by side, per region and per clone -- the three-technology extension of
# scDNA/4_sprinter/persample/plot_sphase_persample.R (which only had 10X/DLP+).
#
# Input:  scDNA/4_sprinter/results/sphase_summary_persample.tsv  (10X, DLP+)
#         scRNA/results/sphase_summary.tsv                       (scRNA)
# Output: scMultiOmics/results/sphase_comparison.tsv
#         scMultiOmics/results/sphase_comparison.png
#
# Requires: ggplot2, scales (present in the "R" micromamba env; no Seurat/
# singularity needed -- this script only reads the two summary TSVs above).
#   micromamba run -n R Rscript scMultiOmics/sphase_comparison.R
#
# Must be run from the repository root (paths below are repo-relative), after
# both scDNA/4_sprinter/persample/summarize_sphase_persample.R and
# scRNA/1_cellcycle.R have produced their summary TSVs.
#
# IMPORTANT CAVEAT -- read before trusting any single clone-level bar here:
# scDNA clone_label ("R1_1", "R4_2", ...) comes from k-means on CNA profiles;
# scRNA's clusters come from an independent transcriptome/inferCNV-based
# clustering (see scRNA/README.md). These are two unrelated partitions of the
# malignant cells. Where their labels happen to share a region+index (this
# only occurs for R4, by coincidence of both having 3 subclusters), that does
# NOT mean "R4_1" refers to the same set of cells in both technologies -- no
# cross-technology clone matching was performed. Where scRNA does not
# subdivide a region into multiple clusters at all (P, R2, R3, R5 -- see
# scRNA/README.md), its one malignant cluster is kept as its own
# "<region>_all" bar rather than folded into any specific scDNA clone.

library(ggplot2)
library(scales)

# ── Paths ──────────────────────────────────────────────────────────────────

SCDNA_TSV <- "scDNA/4_sprinter/results/sphase_summary_persample.tsv"
SCRNA_TSV <- "scRNA/results/sphase_summary.tsv"

OUT_TSV <- "scMultiOmics/results/sphase_comparison.tsv"
OUT_PNG <- "scMultiOmics/results/sphase_comparison.png"

TECH_LEVELS <- c("10X", "DLP", "scRNA")
MIN_N_ALPHA <- 15   # matches the convention used throughout scDNA/4_sprinter and scRNA

dir.create(dirname(OUT_TSV), recursive = TRUE, showWarnings = FALSE)

# ── Load scDNA (10X + DLP+), keep region-native clones + healthy ───────────
# Same filter as scDNA/4_sprinter/persample/plot_sphase_persample.R.

scdna <- read.table(SCDNA_TSV, header = TRUE, sep = "\t",
                    colClasses = c(region = "character", technology = "character",
                                   clone_label = "character", cell_type = "character",
                                   n_cells = "integer", n_sphase = "integer",
                                   sphase_fraction = "numeric"))

scdna <- scdna[scdna$clone_label == "healthy" |
               startsWith(scdna$clone_label, paste0(scdna$region, "_")), ]
scdna <- scdna[, c("region", "technology", "clone_label", "cell_type",
                    "n_cells", "n_sphase", "sphase_fraction")]

# ── Load scRNA, derive a comparable clone_label, keep region-native only ───

scrna_raw <- read.table(SCRNA_TSV, header = TRUE, sep = "\t",
                        colClasses = c(region = "character", group = "character",
                                       class = "character", n_cells = "integer",
                                       n_sphase = "integer", sphase_fraction = "numeric"))

tag <- sub("^Malignant_", "", scrna_raw$group)
native_region <- sub("_.*", "", tag)
is_native <- scrna_raw$group == "healthy" | native_region == scrna_raw$region
scrna <- scrna_raw[is_native, ]
tag <- tag[is_native]

# A tag equal to its own region (e.g. "P", "R2") means scRNA did not subdivide
# that region into multiple clusters -- label it "<region>_all" so it reads
# as a whole-region estimate, not a specific numbered clone.
undivided <- tag == scrna$region
clone_label <- ifelse(scrna$group == "healthy", "healthy",
                ifelse(undivided, paste0(scrna$region, "_all"), tag))

scrna <- data.frame(
  region          = scrna$region,
  technology      = "scRNA",
  clone_label     = clone_label,
  cell_type       = ifelse(clone_label == "healthy", "healthy", "tumour"),
  n_cells         = scrna$n_cells,
  n_sphase        = scrna$n_sphase,
  sphase_fraction = scrna$sphase_fraction,
  stringsAsFactors = FALSE
)

cat("scDNA rows (10X + DLP, region-native + healthy):", nrow(scdna), "\n")
cat("scRNA rows (region-native + healthy):", nrow(scrna), "\n")
cat("scRNA whole-region ('_all') labels:",
    paste(sort(unique(clone_label[undivided & scrna_raw$group[is_native] != "healthy"])), collapse = ", "), "\n")

# ── Combine ──────────────────────────────────────────────────────────────────

dat <- rbind(scdna, scrna)
dat$region <- factor(dat$region, levels = c("P", "R1", "R2", "R3", "R4", "R5"))

# x_order: healthy = 0, "<region>_all" = 999 (always last), else numeric suffix
dat$x_order <- 0L
tumour_rows <- dat$clone_label != "healthy"
suffix <- sub(".*_", "", dat$clone_label[tumour_rows])
suffix_int <- suppressWarnings(as.integer(suffix))   # "all" -> NA by design, replaced below
dat$x_order[tumour_rows] <- ifelse(suffix == "all", 999L, suffix_int)

x_levels <- unique(dat[order(dat$region, dat$x_order), "clone_label"])
dat$clone_label <- factor(dat$clone_label, levels = x_levels)

# ── Complete to strictly 3 bars (10X / DLP / scRNA) per clone per region ───

group_keys <- unique(data.frame(clone_label = as.character(dat$clone_label),
                                 region      = as.character(dat$region),
                                 stringsAsFactors = FALSE))

all_combos <- do.call(rbind, lapply(seq_len(nrow(group_keys)), function(i) {
  data.frame(clone_label = group_keys$clone_label[i],
             region      = group_keys$region[i],
             technology  = TECH_LEVELS,
             stringsAsFactors = FALSE)
}))

dat_vals <- data.frame(clone_label     = as.character(dat$clone_label),
                        region          = as.character(dat$region),
                        technology      = as.character(dat$technology),
                        sphase_fraction = dat$sphase_fraction,
                        n_cells         = dat$n_cells,
                        n_sphase        = dat$n_sphase,
                        stringsAsFactors = FALSE)

dat <- merge(all_combos, dat_vals, by = c("clone_label", "region", "technology"), all.x = TRUE)

dat$sphase_fraction[is.na(dat$sphase_fraction)] <- 0
dat$n_cells[is.na(dat$n_cells)]                 <- 0L
dat$n_sphase[is.na(dat$n_sphase)]                <- 0L

dat$clone_label <- factor(dat$clone_label, levels = x_levels)
dat$technology  <- factor(dat$technology,  levels = TECH_LEVELS)
dat$region      <- factor(dat$region,      levels = c("P", "R1", "R2", "R3", "R4", "R5"))
dat$alpha_val   <- ifelse(dat$n_cells >= MIN_N_ALPHA, 1.0, 0.4)
dat$n_label     <- ifelse(dat$n_cells == 0, "", as.character(dat$n_cells))

dat <- dat[order(dat$region, dat$clone_label, dat$technology), ]
write.table(dat[, c("region", "technology", "clone_label", "n_cells", "n_sphase", "sphase_fraction")],
            OUT_TSV, sep = "\t", quote = FALSE, row.names = FALSE)
cat("Saved:\n ", OUT_TSV, "\n")

# ── Colour palette ─────────────────────────────────────────────────────────
# 10X/DLP colours match the established scDNA/4_sprinter figures; scRNA adds
# a third, well-separated categorical hue (dataviz-skill slot 1).

tech_colours <- c(`10X` = "#E67E22", DLP = "#27AE60", scRNA = "#2a78d6")

# ── Plot ───────────────────────────────────────────────────────────────────

p <- ggplot(dat, aes(x = clone_label, y = sphase_fraction,
                     fill = technology, alpha = alpha_val)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(aes(label = n_label),
            position  = position_dodge(width = 0.8),
            angle = 90, hjust = -0.1, vjust = 0.5,
            size = 2, colour = "grey30", alpha = 1) +
  scale_fill_manual(values = tech_colours, name = "Technology") +
  scale_alpha_identity(guide = "none") +
  scale_y_continuous(labels = percent_format(accuracy = 1),
                     limits = c(0, NA),
                     expand = expansion(mult = c(0, 0.15))) +
  facet_wrap(~ region, nrow = 1, scales = "free_x") +
  labs(x = NULL, y = "S-phase fraction",
       caption = paste(
         "Bar labels = n cells. Bars with n_cells < 15 shown at 40% opacity.",
         "Missing technology modes shown as zero-height bars.",
         "'<region>_all' = scRNA's one undivided malignant cluster for that region.",
         "scDNA and scRNA clone/cluster labels come from independent clusterings --",
         "matching label text (e.g. R4_1) does not imply the same cells (see script header).",
         sep = "\n")) +
  theme_bw(base_size = 11) +
  theme(axis.text.x        = element_text(angle = 45, hjust = 1, size = 8),
        strip.background   = element_rect(fill = "grey92"),
        strip.text         = element_text(face = "bold"),
        panel.grid.major.x = element_blank(),
        legend.position    = "bottom",
        plot.caption       = element_text(size = 7.5, colour = "grey50", hjust = 0))

ggsave(OUT_PNG, p, width = 320, height = 135, units = "mm", dpi = 300)
cat("Saved:\n ", OUT_PNG, "\n")
