# Validate knee-point cell calling against Cell Ranger DNA ground truth.
#
# For each sample, plots the barcode rank curve with both thresholds overlaid:
#   red  = knee-point threshold
#   blue = Cell Ranger DNA threshold (from original analysis RDS)
#
# Also produces a read-count density plot per sample showing the four barcode
# categories: shared (both agree), knee-only, CR-only, and background.
#
# Input:  <OUTDIR>/<sample>/barcode_counts[COUNT_SUFFIX].txt
#         <OUTDIR>/<sample>/valid_barcodes_cellranger.txt
#         <OUTDIR>/<sample>/valid_barcodes_knee[COUNT_SUFFIX].txt   (optional)
# Output: <OUTDIR>/validation_knee_vs_cr[COUNT_SUFFIX].pdf
#
# Set COUNT_SUFFIX <- "_wl" for whitelist-corrected files, or pass before source():
#   Rscript -e 'COUNT_SUFFIX <- "_wl"; source("validate_cell_calling.R")'

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

OUTDIR  <- "/srv/home/aste0033/projects/MPNST/scDNA/10X/bwa_output"
SAMPLES <- c("FIT208A3", "FIT208A4", "FIT208A5", "FIT208A6", "FIT208A7", "FIT208A8")
if (!exists("COUNT_SUFFIX")) COUNT_SUFFIX <- ""

plots <- list()

for (samp in SAMPLES) {
  counts_file <- file.path(OUTDIR, samp, paste0("barcode_counts", COUNT_SUFFIX, ".txt"))
  cr_file     <- file.path(OUTDIR, samp, "valid_barcodes_cellranger.txt")
  knee_file   <- file.path(OUTDIR, samp, paste0("valid_barcodes_knee", COUNT_SUFFIX, ".txt"))

  if (!file.exists(counts_file) || !file.exists(cr_file)) {
    message(samp, ": missing ", basename(counts_file), " or valid_barcodes_cellranger.txt — skipping")
    next
  }

  bc <- read.table(counts_file, col.names = c("count", "barcode"))
  bc <- bc[order(bc$count, decreasing = TRUE), ]
  bc$rank <- seq_len(nrow(bc))

  cr_bc   <- readLines(cr_file)
  cr_rank <- length(cr_bc)        # Cell Ranger DNA threshold = rank of last called cell
  cr_thr  <- bc$count[cr_rank]

  has_knee <- file.exists(knee_file)
  if (has_knee) {
    knee_bc   <- readLines(knee_file)
    knee_rank <- length(knee_bc)
    knee_thr  <- bc$count[knee_rank]
  }

  # ── Category label per barcode ──────────────────────────────────────────────
  bc$category <- "background"
  bc$category[bc$barcode %in% cr_bc]                              <- "CR only"
  if (has_knee) {
    bc$category[bc$barcode %in% knee_bc]                          <- "knee only"
    bc$category[bc$barcode %in% cr_bc & bc$barcode %in% knee_bc] <- "shared"
  }
  bc$category <- factor(bc$category,
                        levels = c("shared", "knee only", "CR only", "background"))

  cat_colours <- c(shared     = "#2ca02c",
                   "knee only" = "#d62728",
                   "CR only"   = "#1f77b4",
                   background  = "grey80")

  # ── Rank curve with both thresholds ─────────────────────────────────────────
  p_rank <- ggplot(bc, aes(rank, count)) +
    geom_line(linewidth = 0.3, colour = "grey40") +
    geom_vline(xintercept = cr_rank, colour = "#1f77b4", linewidth = 0.7,
               linetype = "dashed") +
    annotate("text", x = cr_rank * 1.3, y = max(bc$count) * 0.7,
             label = paste0("CR: ", cr_rank, " cells\n≥", cr_thr, " reads"),
             colour = "#1f77b4", hjust = 0, size = 2.8) +
    {if (has_knee) list(
      geom_vline(xintercept = knee_rank, colour = "#d62728", linewidth = 0.7,
                 linetype = "dashed"),
      annotate("text", x = knee_rank * 1.3, y = max(bc$count) * 0.45,
               label = paste0("knee: ", knee_rank, " cells\n≥", knee_thr, " reads"),
               colour = "#d62728", hjust = 0, size = 2.8)
    )} +
    scale_x_log10(labels = scales::comma) +
    scale_y_log10(labels = scales::comma) +
    labs(title = samp, x = "Barcode rank", y = "Read count") +
    theme_bw(base_size = 9)

  # ── Density of read counts by category ──────────────────────────────────────
  # Limit to barcodes with ≥ cr_thr/10 reads to keep the plot readable
  bc_sub <- bc[bc$count >= max(2, cr_thr / 10), ]
  p_dens <- ggplot(bc_sub, aes(x = count, fill = category, colour = category)) +
    geom_density(alpha = 0.35, linewidth = 0.4) +
    scale_x_log10(labels = scales::comma) +
    scale_fill_manual(values = cat_colours, drop = FALSE) +
    scale_colour_manual(values = cat_colours, drop = FALSE) +
    labs(x = "Read count", y = "Density", fill = NULL, colour = NULL) +
    theme_bw(base_size = 9) +
    theme(legend.position = "bottom")

  plots[[paste0(samp, "_rank")]] <- p_rank
  plots[[paste0(samp, "_dens")]] <- p_dens
}

# ── Save PDF: pairs of plots (rank | density) per sample ─────────────────────
n_samples <- length(SAMPLES)
pdf_path  <- file.path(OUTDIR, paste0("validation_knee_vs_cr", COUNT_SUFFIX, ".pdf"))
pdf(pdf_path, width = 14, height = 5 * n_samples)

for (samp in SAMPLES) {
  rk <- plots[[paste0(samp, "_rank")]]
  dn <- plots[[paste0(samp, "_dens")]]
  if (!is.null(rk) && !is.null(dn)) print(rk + dn)
}

invisible(dev.off())
message("Validation plots → ", pdf_path)
message("")
message("How to read the plots:")
message("  Rank curve  — if the blue (CR) line is far to the right of the red (knee)")
message("                line, the knee algorithm is cutting off too early.")
message("  Density     — 'CR only' cells should have a read-count distribution")
message("                clearly above 'background'. If they do, they are real cells")
message("                the knee algorithm missed, and you should lower the threshold.")
