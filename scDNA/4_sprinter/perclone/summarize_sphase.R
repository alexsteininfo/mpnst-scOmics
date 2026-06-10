#!/usr/bin/env Rscript
# Summarize S-phase fractions from per-clone SPRINTER outputs.
#
# Reads all sprinter.output.tsv.gz files under the perclone output directory,
# computes per-group S-phase fractions, and writes a tidy summary table.
#
# Input:   /srv/home/aste0033/projects/MPNST/scDNA/SPRINTER/perclone/
# Output:  scDNA/4_sprinter/results/sphase_summary.tsv

# ── Paths ─────────────────────────────────────────────────────────────────────

BASE_OUTDIR <- "/srv/home/aste0033/projects/MPNST/scDNA/SPRINTER/perclone"
METADATA    <- "/srv/home/aste0033/GitHub/mpnst-scOmics/scDNA/2_clustering/metadata_cells.tsv"
OUT_FILE    <- "scDNA/4_sprinter/results/sphase_summary.tsv"

# ── Find all SPRINTER output files ────────────────────────────────────────────

out_files <- list.files(BASE_OUTDIR, pattern = "sprinter\\.output\\.tsv\\.gz$",
                        recursive = TRUE, full.names = TRUE)

cat("Found", length(out_files), "SPRINTER output files\n")

if (length(out_files) == 0) {
  stop("No SPRINTER outputs found in ", BASE_OUTDIR,
       "\nRun perclone/run_sprinter_perclone.sh first.")
}

# ── Load metadata for group annotations ───────────────────────────────────────

meta <- read.table(METADATA, header = TRUE, sep = "\t",
                   colClasses = c(barcode         = "character",
                                  region          = "character",
                                  technology      = "character",
                                  kmeans_cluster  = "integer",
                                  clone_label     = "character",
                                  cell_type       = "character",
                                  chisel_rdr_path = "character"))

# Clone label → cell_type lookup (one row per unique clone_label)
clone_meta <- unique(meta[, c("clone_label", "cell_type")])

# ── Process each output file ──────────────────────────────────────────────────

results <- lapply(out_files, function(f) {
  # Parse group_name and mode from the path:
  # .../perclone/<group_name>/<mode>/sprinter.output.tsv.gz
  parts     <- strsplit(f, .Platform$file.sep)[[1]]
  mode      <- parts[length(parts) - 1]
  group_name <- parts[length(parts) - 2]

  # Parse clone_label from group_name:
  # "healthy_R1"  → clone_label = "healthy",  region = "R1"
  # "clone_R1_1"  → clone_label = "R1_1",     region = NA (spans regions)
  if (startsWith(group_name, "healthy_")) {
    clone_label <- "healthy"
    region      <- sub("^healthy_", "", group_name)
  } else {
    clone_label <- sub("^clone_", "", group_name)
    region      <- NA_character_
  }

  tryCatch({
    df <- read.table(gzfile(f), header = TRUE, sep = "\t",
                     check.names = FALSE)

    # The S-phase column is named "IS-S-PHASE" (with hyphen)
    if (!"IS-S-PHASE" %in% names(df)) {
      warning("No IS-S-PHASE column in ", f, " — skipping")
      return(NULL)
    }
    if (!"CELL" %in% names(df)) {
      warning("No CELL column in ", f, " — skipping")
      return(NULL)
    }

    # One row per cell (SPRINTER outputs one row per cell in the summary)
    # IS-S-PHASE is a Python boolean string ("True"/"False"), not R logical
    cell_df  <- unique(df[, c("CELL", "IS-S-PHASE")])
    n_cells  <- nrow(cell_df)
    n_sphase <- sum(cell_df[["IS-S-PHASE"]] == "True", na.rm = TRUE)

    data.frame(
      group_name      = group_name,
      clone_label     = clone_label,
      region          = region,
      technology_mode = mode,
      n_cells         = n_cells,
      n_sphase        = n_sphase,
      sphase_fraction = if (n_cells > 0) n_sphase / n_cells else NA_real_,
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    warning("Error reading ", f, ": ", conditionMessage(e))
    NULL
  })
})

results <- do.call(rbind, Filter(Negate(is.null), results))

if (is.null(results) || nrow(results) == 0) stop("All output files failed to parse.")

# ── Annotate with cell_type ───────────────────────────────────────────────────

results <- merge(results, clone_meta, by = "clone_label", all.x = TRUE)

# ── Order and write ───────────────────────────────────────────────────────────

results <- results[, c("group_name", "cell_type", "clone_label", "region",
                        "technology_mode", "n_cells", "n_sphase", "sphase_fraction")]
results <- results[order(results$cell_type, results$clone_label, results$technology_mode), ]

write.table(results, OUT_FILE, sep = "\t", quote = FALSE, row.names = FALSE)

cat("\nS-phase summary written to:", OUT_FILE, "\n")
cat("Rows:", nrow(results), "\n\n")
print(results, row.names = FALSE)
