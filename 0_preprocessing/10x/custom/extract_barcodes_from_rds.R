# Extract Cell Ranger DNA-called barcodes from the original analysis RDS and
# compare them against the knee-point calls produced by knee_point_cell_calling.R.
#
# The RDS matrix (cells × bins) was produced by the original Haixi analysis
# using Cell Ranger DNA. Its rownames encode the called cells in the format
# "REGION_BARCODE" (e.g. "R1_AAACCTGCAACCTAGT"). These are the exact cells
# used in every downstream script in the original analysis.
#
# Input:  results/pcf/MPNST_scDNA_nopcf_raw_mtx.rds          (read-only)
#         <OUTDIR>/<sample>/valid_barcodes_knee.txt           (optional, for comparison)
# Output: <OUTDIR>/<sample>/valid_barcodes_cellranger.txt     (one 16-mer per line)
#         <OUTDIR>/cellranger_barcodes_summary.tsv

RDS_PATH <- "/mnt/iribhm/people/mtarabichi/MPNST/pvanloo_20250618/Haixi/10x_DNA/results/pcf/MPNST_scDNA_nopcf_raw_mtx.rds"
OUTDIR   <- "/srv/home/aste0033/projects/MPNST/scDNA/10X/bwa_output"

# Set KNEE_SUFFIX to "_wl" to compare against whitelist-corrected knee calls
# (produced by knee_point_cell_calling.R with COUNT_SUFFIX <- "_wl"). Leave "" for raw.
# Can be overridden by the caller: Rscript -e 'KNEE_SUFFIX <- "_wl"; source(...)'
if (!exists("KNEE_SUFFIX")) KNEE_SUFFIX <- ""

REGION_TO_SAMPLE <- c(
  R1 = "FIT208A3",
  R2 = "FIT208A4",
  R3 = "FIT208A5",
  R4 = "FIT208A6",
  R5 = "FIT208A7",
  P  = "FIT208A8"
)

# ── Extract barcodes from RDS ─────────────────────────────────────────────────
cat("Loading RDS...\n")
mtx          <- readRDS(RDS_PATH)
all_barcodes <- rownames(mtx)   # "R1_AAACCTGCAACCTAGT", "P_TGACGGCAGTGCAAGC", ...
cat("Matrix:", nrow(mtx), "cells ×", ncol(mtx), "bins\n\n")

summary_rows <- list()

for (region in names(REGION_TO_SAMPLE)) {
  sample <- REGION_TO_SAMPLE[[region]]
  mask   <- startsWith(all_barcodes, paste0(region, "_"))

  # Strip region prefix to recover bare 16-mer matching CB tags in our BAMs
  barcodes_16mer <- sub(paste0("^", region, "_"), "", all_barcodes[mask])

  outfile <- file.path(OUTDIR, sample, "valid_barcodes_cellranger.txt")
  writeLines(barcodes_16mer, outfile)

  summary_rows[[sample]] <- data.frame(
    sample             = sample,
    region             = region,
    cellranger_n_cells = length(barcodes_16mer),
    stringsAsFactors   = FALSE
  )
  message(sample, " (", region, "): ", length(barcodes_16mer),
          " cells → ", outfile)
}

# ── Save summary ──────────────────────────────────────────────────────────────
summary_df <- do.call(rbind, summary_rows)
rownames(summary_df) <- NULL

summary_path <- file.path(OUTDIR, "cellranger_barcodes_summary.tsv")
write.table(summary_df, summary_path, sep = "\t", quote = FALSE, row.names = FALSE)

# ── Compare with knee-point calls (if available) ──────────────────────────────
knee_files <- setNames(
  file.path(OUTDIR, summary_df$sample, paste0("valid_barcodes_knee", KNEE_SUFFIX, ".txt")),
  summary_df$sample
)
cr_files <- setNames(
  file.path(OUTDIR, summary_df$sample, "valid_barcodes_cellranger.txt"),
  summary_df$sample
)

if (any(file.exists(knee_files))) {
  cat("\n── Comparison: Cell Ranger DNA vs knee-point ────────────────────────────\n")
  comp_rows <- list()
  for (samp in summary_df$sample) {
    cr_bc <- readLines(cr_files[[samp]])
    if (!file.exists(knee_files[[samp]])) {
      comp_rows[[samp]] <- data.frame(
        sample        = samp,
        cellranger    = length(cr_bc),
        knee          = NA_integer_,
        overlap       = NA_integer_,
        jaccard       = NA_real_,
        only_cr       = NA_integer_,
        only_knee     = NA_integer_,
        stringsAsFactors = FALSE
      )
      next
    }
    knee_bc <- readLines(knee_files[[samp]])
    both    <- intersect(cr_bc, knee_bc)
    union_  <- union(cr_bc, knee_bc)
    comp_rows[[samp]] <- data.frame(
      sample     = samp,
      cellranger = length(cr_bc),
      knee       = length(knee_bc),
      overlap    = length(both),
      jaccard    = round(length(both) / length(union_), 3),
      only_cr    = length(setdiff(cr_bc, knee_bc)),
      only_knee  = length(setdiff(knee_bc, cr_bc)),
      stringsAsFactors = FALSE
    )
  }
  comp_df <- do.call(rbind, comp_rows)
  rownames(comp_df) <- NULL
  print(comp_df)

  comp_path <- file.path(OUTDIR, paste0("cell_calling_comparison", KNEE_SUFFIX, ".tsv"))
  write.table(comp_df, comp_path, sep = "\t", quote = FALSE, row.names = FALSE)
  message("\nComparison → ", comp_path)
  message("jaccard = overlap / union; 1.0 = perfect agreement")
} else {
  cat("\nNo valid_barcodes_knee", KNEE_SUFFIX, ".txt files found — skipping comparison.\n", sep = "")
  cat("Run knee_point_cell_calling.R first if you want a comparison.\n")
}

cat("\n── Cell Ranger DNA summary ──────────────────────────────────────────────\n")
print(summary_df)
message("\nSummary → ", summary_path)
