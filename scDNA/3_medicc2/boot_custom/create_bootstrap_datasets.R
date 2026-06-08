#!/usr/bin/env Rscript
# create_bootstrap_datasets.R
#
# Generate bootstrap TSV datasets from a MEDICC2 input TSV for downstream
# full-pipeline CNA burden estimation with proper integer branch lengths.
#
# Two resampling methods:
#   chr-wise  — sample chromosomes with replacement, preserving within-chromosome
#               correlation of CNA events. Sequential integer labels (1..N) are
#               assigned to each draw so MEDICC2 always sees a valid genome.
#   cell-wise — sample cells with replacement; duplicate draws are renamed with
#               _copy2, _copy3 suffixes. Assesses sensitivity to cell sampling.
#
# Usage: micromamba run -n R Rscript create_bootstrap_datasets.R
#
# Input TSV columns (with header): sample_id, chrom, start, end, cn_a, cn_b
#
# Output:
#   <out_dir>/<method>/<sample_name>/boot_001.tsv ... boot_N.tsv
#   <out_dir>/<method>/<sample_name>/metadata.tsv
#     Chr-wise  columns: replicate, draw_index, new_chrom, original_chrom
#     Cell-wise columns: replicate, draw_index, new_cell_id, original_cell_id

suppressPackageStartupMessages(library(data.table))

# ---------------------------------------------------------------------------
# Configuration — adjust before running
# ---------------------------------------------------------------------------

minseg    <- "1.5"        # "1.5" or "2.5"
method    <- "chr-wise"  # "chr-wise" or "cell-wise"
n_boot    <- 100L
seed      <- 42L

input_dir <- sprintf("/srv/home/aste0033/projects/MPNST/Haixi/scDNA/MEDICC2_input/per_cluster_10X_DLP/minseg_%s", minseg)
out_dir   <- "/srv/home/aste0033/projects/MPNST/scDNA/MEDICC2/bootstrap_custom/datasets"

# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------

if (!method %in% c("chr-wise", "cell-wise"))
  stop("method must be 'chr-wise' or 'cell-wise'")
if (!dir.exists(input_dir))
  stop("Input directory not found: ", input_dir)

tsv_files <- list.files(input_dir, pattern = "\\.tsv$", full.names = TRUE)
if (length(tsv_files) == 0)
  stop("No TSV files found in: ", input_dir)

cat(sprintf("minseg: %s | method: %s | n-boot: %d | seed: %d\n", minseg, method, n_boot, seed))
cat(sprintf("Processing %d TSV files from %s\n\n", length(tsv_files), input_dir))

# ---------------------------------------------------------------------------
# Loop over all cluster TSVs
# ---------------------------------------------------------------------------

for (input in tsv_files) {

  dt          <- fread(input)
  sample_name <- tools::file_path_sans_ext(basename(input))
  out_sub     <- file.path(out_dir, method, sample_name)
  dir.create(out_sub, recursive = TRUE, showWarnings = FALSE)

  set.seed(seed)

  cat(sprintf("[%s] %d rows, %d cells -> %s\n",
              basename(input), nrow(dt), length(unique(dt$sample_id)), out_sub))

  # -------------------------------------------------------------------------
  # Chr-wise bootstrap
  # -------------------------------------------------------------------------

  if (method == "chr-wise") {

    chroms    <- sort(unique(dt$chrom))
    n_chr     <- length(chroms)
    chr_parts <- split(dt, by = "chrom", keep.by = TRUE)
    names(chr_parts) <- as.character(sort(unique(dt$chrom)))

    meta_list <- vector("list", n_boot)

    for (i in seq_len(n_boot)) {
      draw <- sample(chroms, n_chr, replace = TRUE)

      pieces <- vector("list", n_chr)
      for (j in seq_len(n_chr)) {
        piece        <- copy(chr_parts[[as.character(draw[j])]])
        piece[, chrom := j]
        pieces[[j]]  <- piece
      }
      result <- rbindlist(pieces)
      fwrite(result, file.path(out_sub, sprintf("boot_%03d.tsv", i)), sep = "\t")

      meta_list[[i]] <- data.table(
        replicate      = i,
        draw_index     = seq_len(n_chr),
        new_chrom      = seq_len(n_chr),
        original_chrom = draw
      )

      if (i %% 10 == 0) cat(sprintf("  %d / %d replicates done\n", i, n_boot))
    }

    meta <- rbindlist(meta_list)
    fwrite(meta, file.path(out_sub, "metadata.tsv"), sep = "\t")
    cat(sprintf("  Done. %d chr-wise bootstrap TSVs written.\n\n", n_boot))

  # -------------------------------------------------------------------------
  # Cell-wise bootstrap
  # -------------------------------------------------------------------------

  } else {

    cells      <- unique(dt$sample_id)
    n_cell     <- length(cells)
    cell_parts <- split(dt, by = "sample_id", keep.by = TRUE)

    meta_list <- vector("list", n_boot)

    for (i in seq_len(n_boot)) {
      draw <- sample(cells, n_cell, replace = TRUE)

      current_count        <- setNames(integer(length(cells)), cells)
      new_ids              <- character(n_cell)
      for (j in seq_len(n_cell)) {
        orig               <- draw[j]
        current_count[orig] <- current_count[orig] + 1L
        k                  <- current_count[orig]
        new_ids[j]         <- if (k == 1L) orig else paste0(orig, "_copy", k)
      }

      pieces <- vector("list", n_cell)
      for (j in seq_len(n_cell)) {
        piece              <- copy(cell_parts[[draw[j]]])
        piece[, sample_id := new_ids[j]]
        pieces[[j]]        <- piece
      }
      result <- rbindlist(pieces)
      fwrite(result, file.path(out_sub, sprintf("boot_%03d.tsv", i)), sep = "\t")

      meta_list[[i]] <- data.table(
        replicate        = i,
        draw_index       = seq_len(n_cell),
        new_cell_id      = new_ids,
        original_cell_id = draw
      )

      if (i %% 10 == 0) cat(sprintf("  %d / %d replicates done\n", i, n_boot))
    }

    meta <- rbindlist(meta_list)
    fwrite(meta, file.path(out_sub, "metadata.tsv"), sep = "\t")
    cat(sprintf("  Done. %d cell-wise bootstrap TSVs written.\n\n", n_boot))
  }

}
