# Knee-point cell calling for FIT208 10x scDNA samples.
#
# Input:  <OUTDIR>/<sample>/barcode_counts.txt  (from extract_barcode_counts.sh)
# Output: <OUTDIR>/<sample>/valid_barcodes.txt  (one barcode per line)
#         <OUTDIR>/knee_plots.pdf               (6-panel knee plot)
#         <OUTDIR>/cell_calling_summary.tsv     (per-sample summary table)
#
# After inspecting the plots, adjust MANUAL_THRESHOLDS if the automatic knee
# misses the inflection point, then re-run. NA = use automatic detection.

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

OUTDIR <- "/srv/home/aste0033/projects/MPNST/scDNA/10X/bwa_output"

# Set COUNT_SUFFIX to "_wl" to use whitelist-corrected barcode_counts_wl.txt
# (produced by apply_whitelist_to_counts.py). Leave "" for raw counts.
# Can be overridden by the caller: Rscript -e 'COUNT_SUFFIX <- "_wl"; source(...)'
if (!exists("COUNT_SUFFIX")) COUNT_SUFFIX <- ""

SAMPLES <- c("FIT208A3", "FIT208A4", "FIT208A5", "FIT208A6", "FIT208A7", "FIT208A8")

EXPECTED_CELLS <- c(
  FIT208A3 = 2220,
  FIT208A4 = 1440,
  FIT208A5 = 106,
  FIT208A6 = 363,
  FIT208A7 = 1513,
  FIT208A8 = 1154
)

# Thresholds derived from the original Cell Ranger DNA cell calls (read from the
# analysis RDS). The automatic Kneedle algorithm finds a knee within the cell
# population rather than at the cell/empty-droplet boundary for this data type,
# because scDNA cells have a wide read-count spread. These manual thresholds match
# the last Cell Ranger DNA-called cell in each sample's sorted barcode list.
MANUAL_THRESHOLDS <- c(
  FIT208A3 = 6812,
  FIT208A4 = 5658,
  FIT208A5 = 10319,
  FIT208A6 = 10026,
  FIT208A7 = 8361,
  FIT208A8 = 13863
)

# Kneedle algorithm: returns the rank index of the knee.
# Normalises the log10(rank) vs log10(count) curve to [0,1] and finds the
# point with maximum perpendicular distance above the first-to-last diagonal.
find_knee <- function(counts) {
  x <- log10(seq_along(counts))
  y <- log10(counts)
  x_n <- (x - min(x)) / (max(x) - min(x))
  y_n <- (y - min(y)) / (max(y) - min(y))
  d <- (x_n + y_n - 1) / sqrt(2)   # positive = above diagonal
  which.max(d)
}

plots       <- list()
summary_rows <- list()

for (samp in SAMPLES) {
  bc_file <- file.path(OUTDIR, samp, paste0("barcode_counts", COUNT_SUFFIX, ".txt"))
  if (!file.exists(bc_file)) {
    message(samp, ": ", basename(bc_file), " not found — skipping")
    next
  }

  bc       <- read.table(bc_file, col.names = c("count", "barcode"))
  bc       <- bc[order(bc$count, decreasing = TRUE), ]
  bc$rank  <- seq_len(nrow(bc))

  manual_thr <- MANUAL_THRESHOLDS[[samp]]
  if (!is.na(manual_thr)) {
    knee_rank <- max(1L, which(bc$count < manual_thr)[1] - 1L)
    method    <- "manual"
  } else {
    knee_rank <- find_knee(bc$count)
    method    <- "knee"
  }

  threshold <- bc$count[knee_rank]
  cells     <- bc$barcode[bc$count >= threshold]

  writeLines(cells, file.path(OUTDIR, samp, paste0("valid_barcodes_knee", COUNT_SUFFIX, ".txt")))

  summary_rows[[samp]] <- data.frame(
    sample          = samp,
    total_barcodes  = nrow(bc),
    cells_called    = length(cells),
    expected_cells  = EXPECTED_CELLS[[samp]],
    threshold_reads = threshold,
    method          = method,
    stringsAsFactors = FALSE
  )

  # Knee plot (log-log)
  plots[[samp]] <- ggplot(bc, aes(rank, count)) +
    geom_line(linewidth = 0.4) +
    geom_vline(xintercept = knee_rank, color = "firebrick", linetype = "dashed") +
    annotate(
      "text",
      x     = knee_rank * 2,
      y     = max(bc$count) * 0.5,
      label = paste0(scales::comma(length(cells)), " cells\n≥", scales::comma(threshold), " reads"),
      hjust = 0,
      size  = 3
    ) +
    scale_x_log10(labels = scales::comma) +
    scale_y_log10(labels = scales::comma) +
    labs(
      title    = paste0(samp, "  [", method, "]"),
      subtitle = paste0("expected ~", scales::comma(EXPECTED_CELLS[[samp]]), " cells"),
      x        = "Barcode rank",
      y        = "Read count"
    ) +
    theme_bw(base_size = 10)

  message(samp, ": ", length(cells), " cells called (threshold = ",
          threshold, " reads, method = ", method, ")")
}

# ── Save combined knee plot ────────────────────────────────────────────────────
pdf_path <- file.path(OUTDIR, paste0("knee_plots", COUNT_SUFFIX, ".pdf"))
pdf(pdf_path, width = 14, height = 8)
print(wrap_plots(plots, ncol = 3))
invisible(dev.off())
message("Knee plots → ", pdf_path)

# ── Save summary table ─────────────────────────────────────────────────────────
summary_df    <- do.call(rbind, summary_rows)
rownames(summary_df) <- NULL
summary_path  <- file.path(OUTDIR, paste0("cell_calling_summary", COUNT_SUFFIX, ".tsv"))
write.table(summary_df, summary_path, sep = "\t", quote = FALSE, row.names = FALSE)

cat("\n── Cell calling summary ─────────────────────────────────────────────────\n")
print(summary_df)
message("\nSummary → ", summary_path)
message("Valid barcode lists → ", OUTDIR, "/<sample>/valid_barcodes_knee", COUNT_SUFFIX, ".txt")
message("\nIf cell counts differ substantially from expected, set MANUAL_THRESHOLDS")
message("for affected samples and re-run this script.")
