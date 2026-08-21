#!/usr/bin/env Rscript
# Export raw snRNA-seq counts (Tumor cells only) for TreeAlign, restricted to
# the gene set that survived 0_prepare_scdna_cn.R's segment->gene mapping.
#
# Input:
#   - /srv/home/aste0033/projects/MPNST/zenodo_phase_1/MPNST-Zenodo/data/snRNA/seurat_objects/MPNST_C_updated.rds
#   - <TreeAlign inputs dir>/gene_by_cell_cnv.csv.gz (from 0_prepare_scdna_cn.R; run that first)
#
# Output (under OUT_DIR, see config below), 10x-style sparse triplet:
#   - expr_counts.mtx.gz    gene x cell raw RNA counts (Matrix::writeMM)
#   - expr_genes.tsv.gz     row names, one gene per line, same order as the matrix
#   - expr_cells.tsv.gz     column names, one cell barcode per line
#   - scrna_cell_metadata.csv  cell_id, region, cell_type — for 4_compare_results.R,
#     so it doesn't need to reload the 1.9GB Seurat object a third time
#
# Requires: Seurat via the single-cell.sif container (nothing scRNA-related is
# installed in any micromamba env — see scRNA/README.md for why).
# Usage:
#   singularity exec /mnt/iribhm/software/singularity/single-cell.sif \
#     Rscript scMultiOmics/2_TreeAlign/1_prepare_scrna_expr.R

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(data.table)
})

# ── Configuration ─────────────────────────────────────────────────────────────

MODE <- Sys.getenv("TREEALIGN_MODE", unset = "test")  # "test" or "full"
stopifnot(MODE %in% c("test", "full"))

TEST_REGION <- "R1"  # test mode: restrict scRNA cells to this physical region
                      # (barcode prefix) — the only axis available on the scRNA
                      # side, since scRNA cells have no scDNA clone_label yet;
                      # this is exactly what TreeAlign is being run to determine.

SEURAT_RDS <- "/srv/home/aste0033/projects/MPNST/zenodo_phase_1/MPNST-Zenodo/data/snRNA/seurat_objects/MPNST_C_updated.rds"
CN_DIR  <- file.path("/srv/home/aste0033/projects/MPNST/scMultiOmics/TreeAlign/inputs", MODE)
OUT_DIR <- file.path("/srv/home/aste0033/projects/MPNST/scMultiOmics/TreeAlign/inputs", MODE)
RESULTS_DIR <- "scMultiOmics/2_TreeAlign/results"

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)

cat(sprintf("[1_prepare_scrna_expr.R] MODE = %s, output -> %s\n", MODE, OUT_DIR))

# ── 1. Gene set surviving step 0 (defines which genes we bother exporting) ──

cnv_genes <- fread(file.path(CN_DIR, "gene_by_cell_cnv.csv.gz"), select = "gene")$gene
cat(sprintf("Gene set from 0_prepare_scdna_cn.R: %d genes\n", length(cnv_genes)))

# ── 2. Load the Seurat object (large — ~1.9GB RDS, loaded once) ─────────────

cat("Loading Seurat object (this is a ~1.9GB RDS, may take a couple of minutes)...\n")
obj <- readRDS(SEURAT_RDS)
cat(sprintf("Loaded: %d genes x %d cells, assays: %s\n",
            nrow(obj), ncol(obj), paste(Assays(obj), collapse = ", ")))

# ── 3. Cell selection: Tumor only, optionally region-restricted for test ────

keep_cells <- colnames(obj)[obj$tumor_normal == "Tumor"]
if (MODE == "test") {
  keep_cells <- keep_cells[startsWith(keep_cells, paste0(TEST_REGION, "_"))]
}
cat(sprintf("Cells in scope: %d (Tumor%s)\n", length(keep_cells),
            if (MODE == "test") sprintf(", region %s", TEST_REGION) else ""))

# ── 4. Raw RNA counts, NOT the default SCT assay — see README for why ───────

counts <- GetAssayData(obj, assay = "RNA", layer = "counts")
counts <- counts[, keep_cells, drop = FALSE]

keep_genes <- intersect(rownames(counts), cnv_genes)
cat(sprintf("Genes overlapping the scDNA CN gene set: %d / %d\n", length(keep_genes), length(cnv_genes)))
counts <- counts[keep_genes, , drop = FALSE]

# ── 5. Write as a 10x-style sparse triplet (raw counts are very sparse) ─────

invisible(writeMM(as(counts, "CsparseMatrix"), file.path(OUT_DIR, "expr_counts.mtx")))
invisible(system(sprintf("gzip -f %s", file.path(OUT_DIR, "expr_counts.mtx"))))
fwrite(data.table(gene = rownames(counts)), file.path(OUT_DIR, "expr_genes.tsv.gz"), col.names = FALSE)
fwrite(data.table(cell_id = colnames(counts)), file.path(OUT_DIR, "expr_cells.tsv.gz"), col.names = FALSE)

cell_meta <- data.table(
  cell_id = keep_cells,
  region = sub("_.*$", "", keep_cells),
  cell_type = obj$cell_type[keep_cells]
)
fwrite(cell_meta, file.path(OUT_DIR, "scrna_cell_metadata.csv"))

# ── 6. Small QC summary kept in the repo ────────────────────────────────────

qc <- data.table(
  mode = MODE,
  n_cells = ncol(counts),
  n_genes = nrow(counts),
  pct_nonzero = round(100 * Matrix::nnzero(counts) / (nrow(counts) * ncol(counts)), 3)
)
fwrite(qc, file.path(RESULTS_DIR, sprintf("1_prepare_scrna_expr_qc_%s.tsv", MODE)), sep = "\t")

cat("Done.\n")
print(qc)
