#!/usr/bin/env Rscript
# Summarize S-phase fractions from per-sample SPRINTER outputs.
#
# Each SPRINTER run covers all cells from one region × technology. The output
# file contains one row per cell (CELL, IS-S-PHASE). This script joins those
# results with the metadata to recover per-clone breakdowns within each run.
#
# Input:   /srv/home/aste0033/projects/MPNST/scDNA/SPRINTER/persample/
#          scDNA/2_clustering/metadata_cells.tsv
# Output:  scDNA/4_sprinter/results/sphase_summary_persample.tsv

# ── Paths ─────────────────────────────────────────────────────────────────────

BASE_OUTDIR <- "/srv/home/aste0033/projects/MPNST/scDNA/SPRINTER/persample"
METADATA    <- "/srv/home/aste0033/GitHub/mpnst-scOmics/scDNA/2_clustering/metadata_cells.tsv"
OUT_FILE    <- "scDNA/4_sprinter/results/sphase_summary_persample.tsv"

# ── Find all SPRINTER output files ────────────────────────────────────────────

out_files <- list.files(BASE_OUTDIR, pattern = "sprinter\\.output\\.tsv\\.gz$",
                        recursive = TRUE, full.names = TRUE)

cat("Found", length(out_files), "SPRINTER output files\n")

if (length(out_files) == 0) {
  stop("No SPRINTER outputs found in ", BASE_OUTDIR,
       "\nRun persample/run_sprinter_persample.sh first.")
}

# ── Load metadata ─────────────────────────────────────────────────────────────

meta <- read.table(METADATA, header = TRUE, sep = "\t",
                   colClasses = c(barcode         = "character",
                                  region          = "character",
                                  technology      = "character",
                                  kmeans_cluster  = "integer",
                                  clone_label     = "character",
                                  cell_type       = "character",
                                  chisel_rdr_path = "character"))

# Removed cells are excluded from SPRINTER runs; drop them from the lookup too
meta <- meta[meta$cell_type != "removed", ]

# ── Process each output file ──────────────────────────────────────────────────

results <- lapply(out_files, function(f) {
  # Parse region and technology from path:
  # .../persample/<region>/<technology>/sprinter.output.tsv.gz
  parts      <- strsplit(f, .Platform$file.sep)[[1]]
  technology <- parts[length(parts) - 1]
  region     <- parts[length(parts) - 2]

  tryCatch({
    df <- read.table(gzfile(f), header = TRUE, sep = "\t", check.names = FALSE)

    if (!"IS-S-PHASE" %in% names(df)) {
      warning("No IS-S-PHASE column in ", f, " — skipping")
      return(NULL)
    }
    if (!"CELL" %in% names(df)) {
      warning("No CELL column in ", f, " — skipping")
      return(NULL)
    }

    # One row per cell
    cell_df <- unique(df[, c("CELL", "IS-S-PHASE")])

    # Join with metadata to get clone_label and cell_type per cell
    meta_sub <- meta[meta$region == region & meta$technology == technology,
                     c("barcode", "clone_label", "cell_type")]

    cell_df <- merge(cell_df, meta_sub,
                     by.x = "CELL", by.y = "barcode",
                     all.x = TRUE)

    if (any(is.na(cell_df$clone_label))) {
      n_missing <- sum(is.na(cell_df$clone_label))
      warning(sprintf("%s/%s: %d cells in SPRINTER output not found in metadata — they will be excluded",
                      region, technology, n_missing))
      cell_df <- cell_df[!is.na(cell_df$clone_label), ]
    }

    if (nrow(cell_df) == 0) return(NULL)

    # Group by clone_label
    clones <- unique(cell_df[, c("clone_label", "cell_type")])

    do.call(rbind, lapply(seq_len(nrow(clones)), function(i) {
      cl    <- clones$clone_label[i]
      ctype <- clones$cell_type[i]
      rows  <- cell_df[cell_df$clone_label == cl, ]
      n     <- nrow(rows)
      ns    <- sum(rows[["IS-S-PHASE"]] == "True", na.rm = TRUE)
      data.frame(
        region          = region,
        technology      = technology,
        clone_label     = cl,
        cell_type       = ctype,
        n_cells         = n,
        n_sphase        = ns,
        sphase_fraction = if (n > 0) ns / n else NA_real_,
        stringsAsFactors = FALSE
      )
    }))

  }, error = function(e) {
    warning("Error reading ", f, ": ", conditionMessage(e))
    NULL
  })
})

results <- do.call(rbind, Filter(Negate(is.null), results))

if (is.null(results) || nrow(results) == 0) stop("All output files failed to parse.")

# ── Order and write ───────────────────────────────────────────────────────────

region_levels <- c("P", "R1", "R2", "R3", "R4", "R5")
results$region <- factor(results$region, levels = region_levels)
results <- results[order(results$region, results$technology, results$cell_type,
                         results$clone_label), ]
results$region <- as.character(results$region)

write.table(results, OUT_FILE, sep = "\t", quote = FALSE, row.names = FALSE)

cat("\nS-phase summary (per-sample) written to:", OUT_FILE, "\n")
cat("Rows:", nrow(results), "\n\n")
print(results, row.names = FALSE)
