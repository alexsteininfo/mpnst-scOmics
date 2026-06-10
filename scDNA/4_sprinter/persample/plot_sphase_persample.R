#!/usr/bin/env Rscript
# Plot per-clone S-phase fractions from per-sample SPRINTER runs.
#
# Input:  scDNA/4_sprinter/results/sphase_summary_persample.tsv
# Output: scDNA/4_sprinter/results/sphase_fractions_persample.png
#
# Layout: one facet per region (P | R1–R5). Within each facet, two dodged bars
# per clone label (10X and DLP+). Healthy group is shown first, tumour clones
# follow in clone-number order, unassigned last. Missing technologies appear as
# zero-height bars. Cell counts shown as vertical labels above each bar.
# Groups with n_cells < 15 shown at 40% opacity.

library(ggplot2)
library(scales)

# ── Paths ──────────────────────────────────────────────────────────────────────

IN_FILE <- "scDNA/4_sprinter/results/sphase_summary_persample.tsv"
OUT_PNG <- "scDNA/4_sprinter/results/sphase_fractions_persample.png"

# ── Load and wrangle ───────────────────────────────────────────────────────────

dat <- read.table(IN_FILE, header = TRUE, sep = "\t",
                  colClasses = c(region          = "character",
                                 technology      = "character",
                                 clone_label     = "character",
                                 cell_type       = "character",
                                 n_cells         = "integer",
                                 n_sphase        = "integer",
                                 sphase_fraction = "numeric"))

dat$region <- factor(dat$region, levels = c("P", "R1", "R2", "R3", "R4", "R5"))

# x_order: healthy = 0, tumour clones by numeric suffix, unassigned = Inf
dat$x_order <- ifelse(
  dat$clone_label == "healthy",    0L,
  ifelse(dat$clone_label == "unassigned", 99999L,
         as.integer(sub(".*_", "", dat$clone_label)))
)

dat$alpha_val <- ifelse(dat$n_cells >= 15, 1.0, 0.4)

# Ordered x-axis levels: across all regions, sort by region then x_order
x_levels <- unique(
  dat[order(dat$region, dat$x_order), "clone_label"]
)
dat$clone_label <- factor(dat$clone_label, levels = x_levels)

# ── Complete to strictly 2 bars per clone per region ──────────────────────────

group_keys <- unique(data.frame(
  clone_label = as.character(dat$clone_label),
  region      = as.character(dat$region),
  stringsAsFactors = FALSE
))

all_combos <- do.call(rbind, lapply(seq_len(nrow(group_keys)), function(i) {
  data.frame(
    clone_label = group_keys$clone_label[i],
    region      = group_keys$region[i],
    technology  = c("10X", "DLP"),
    stringsAsFactors = FALSE
  )
}))

dat_vals <- data.frame(
  clone_label     = as.character(dat$clone_label),
  region          = as.character(dat$region),
  technology      = as.character(dat$technology),
  sphase_fraction = dat$sphase_fraction,
  n_cells         = dat$n_cells,
  alpha_val       = dat$alpha_val,
  stringsAsFactors = FALSE
)

dat <- merge(all_combos, dat_vals,
             by    = c("clone_label", "region", "technology"),
             all.x = TRUE)

dat$sphase_fraction[is.na(dat$sphase_fraction)] <- 0
dat$n_cells[is.na(dat$n_cells)]                 <- 0L
dat$alpha_val[is.na(dat$alpha_val)]             <- 0.4

dat$clone_label <- factor(dat$clone_label, levels = x_levels)
dat$technology  <- factor(dat$technology,  levels = c("10X", "DLP"))
dat$region      <- factor(dat$region,      levels = c("P", "R1", "R2", "R3", "R4", "R5"))

# n_cells label: blank for absent modes (n = 0)
dat$n_label <- ifelse(dat$n_cells == 0, "", as.character(dat$n_cells))

# ── Colour palette ─────────────────────────────────────────────────────────────

tech_colours <- c(`10X` = "#E67E22", DLP = "#27AE60")

# ── Plot ───────────────────────────────────────────────────────────────────────

p <- ggplot(dat, aes(x     = clone_label,
                     y     = sphase_fraction,
                     fill  = technology,
                     alpha = alpha_val)) +
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
  labs(x       = NULL,
       y       = "S-phase fraction",
       caption = "Bar labels = n cells. Bars with n_cells < 15 shown at 40% opacity. Missing technology modes shown as zero-height bars.") +
  theme_bw(base_size = 11) +
  theme(axis.text.x        = element_text(angle = 45, hjust = 1, size = 9),
        strip.background   = element_rect(fill = "grey92"),
        strip.text         = element_text(face = "bold"),
        panel.grid.major.x = element_blank(),
        legend.position    = "bottom",
        plot.caption       = element_text(size = 8, colour = "grey50"))

# ── Save ───────────────────────────────────────────────────────────────────────

ggsave(OUT_PNG, p, width = 240, height = 120, units = "mm", dpi = 300)
cat("Saved:\n ", OUT_PNG, "\n")
