#!/usr/bin/env Rscript
# Plot per-clone S-phase fractions from SPRINTER summary.
#
# Input:  scDNA/4_sprinter/results/sphase_summary.tsv
# Output: scDNA/4_sprinter/results/sphase_fractions.png
#
# Layout: one facet per region of origin (P | R1–R5), healthy group first in
# each facet, tumour clones following in clone-number order. Strictly 3 dodged
# bars per group (combined / 10X / DLP+); modes with no data appear as
# zero-height bars. Cell counts shown as vertical labels above each bar.
# Groups with n_cells < 15 shown at 40% opacity.

library(ggplot2)
library(scales)

# ── Paths ──────────────────────────────────────────────────────────────────────

IN_FILE <- "scDNA/4_sprinter/results/sphase_summary.tsv"
OUT_PNG <- "scDNA/4_sprinter/results/sphase_fractions.png"

# ── Load and wrangle ───────────────────────────────────────────────────────────

dat <- read.table(IN_FILE, header = TRUE, sep = "\t",
                  colClasses = c(group_name      = "character",
                                 cell_type       = "character",
                                 clone_label     = "character",
                                 region          = "character",
                                 technology_mode = "character",
                                 n_cells         = "integer",
                                 n_sphase        = "integer",
                                 sphase_fraction = "numeric"))

dat <- dat[dat$clone_label != "unassigned", ]

dat$region_of_origin <- ifelse(
  dat$cell_type == "healthy",
  dat$region,
  sub("_.*", "", dat$clone_label)   # "R1_4" → "R1"
)

dat$group_short <- ifelse(dat$cell_type == "healthy", "healthy", dat$clone_label)

dat$x_order <- ifelse(
  dat$cell_type == "healthy",
  0L,
  as.integer(sub(".*_", "", dat$clone_label))  # "R1_4" → 4
)

dat$alpha_val <- ifelse(dat$n_cells >= 15, 1.0, 0.4)

dat$region_of_origin <- factor(dat$region_of_origin,
                                levels = c("P", "R1", "R2", "R3", "R4", "R5"))

# Ordered x-axis levels: healthy first within each region, clones by number
x_levels    <- unique(dat[order(dat$region_of_origin, dat$x_order), "group_short"])
dat$group_short <- factor(dat$group_short, levels = x_levels)

# ── Complete to strictly 3 bars per group ──────────────────────────────────────

# "healthy" is shared by 6 regional groups — use (group_short, region_of_origin)
# as the composite key so each healthy region is treated independently.
group_keys <- unique(data.frame(
  group_short      = as.character(dat$group_short),
  region_of_origin = as.character(dat$region_of_origin),
  stringsAsFactors = FALSE
))

all_combos <- do.call(rbind, lapply(seq_len(nrow(group_keys)), function(i) {
  data.frame(
    group_short      = group_keys$group_short[i],
    region_of_origin = group_keys$region_of_origin[i],
    technology_mode  = c("combined", "10X", "DLP"),
    stringsAsFactors = FALSE
  )
}))

dat_vals <- data.frame(
  group_short      = as.character(dat$group_short),
  region_of_origin = as.character(dat$region_of_origin),
  technology_mode  = as.character(dat$technology_mode),
  sphase_fraction  = dat$sphase_fraction,
  n_cells          = dat$n_cells,
  alpha_val        = dat$alpha_val,
  stringsAsFactors = FALSE
)

dat <- merge(all_combos, dat_vals,
             by    = c("group_short", "region_of_origin", "technology_mode"),
             all.x = TRUE)

dat$sphase_fraction[is.na(dat$sphase_fraction)] <- 0
dat$n_cells[is.na(dat$n_cells)]                 <- 0L
dat$alpha_val[is.na(dat$alpha_val)]             <- 0.4

dat$group_short      <- factor(dat$group_short,      levels = x_levels)
dat$technology_mode  <- factor(dat$technology_mode,  levels = c("combined", "10X", "DLP"))
dat$region_of_origin <- factor(dat$region_of_origin, levels = c("P", "R1", "R2", "R3", "R4", "R5"))

# n_cells label: blank for absent modes (n = 0 = filled-in missing)
dat$n_label <- ifelse(dat$n_cells == 0, "", as.character(dat$n_cells))

# ── Colour palette ─────────────────────────────────────────────────────────────

tech_colours <- c(combined = "#2C3E50", `10X` = "#E67E22", DLP = "#27AE60")

# ── Plot ───────────────────────────────────────────────────────────────────────

p <- ggplot(dat, aes(x     = group_short,
                     y     = sphase_fraction,
                     fill  = technology_mode,
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
  facet_wrap(~ region_of_origin, nrow = 1, scales = "free_x") +
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
