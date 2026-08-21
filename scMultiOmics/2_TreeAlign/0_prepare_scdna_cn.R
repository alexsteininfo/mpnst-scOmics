#!/usr/bin/env Rscript
# Build TreeAlign's scDNA-side inputs: a gene x cell copy-number matrix, a
# cell -> clone_label table, and a matching (copied/pruned) phylogenetic tree.
#
# Input:
#   - MPNST_all_single_cells_2.5_ds_6039_seed100_final_cn_profiles.tsv
#     (all-cohort MEDICC2 run, minseg=2.5Mb, 10X+DLP combined; per-cell
#     allele-specific segmented CN, plus 6,038 ancestral "internal_N" node
#     profiles and 1 "diploid" root reference which are filtered out below)
#   - MPNST_all_single_cells_2.5_ds_6039_seed100_final_tree.new (same run)
#   - scDNA/2_clustering/metadata_cells.tsv (barcode -> clone_label)
#   - GENCODE v27 GTF (hg38-gatk/gencode.v27.primary_assembly.annotation.gtf)
#
# Output (under OUT_DIR, see config below):
#   - gene_by_cell_cnv.csv.gz   gene x cell total copy number (cn_a + cn_b)
#   - cell_clone_labels.csv     cell_id, clone_id (for CloneAlignClone)
#   - tree.newick               tree pruned/copied to match the cell set above
#   - gene_positions.tsv        gene, chrom, start, end (small QC reference)
#
# Requires: R micromamba env (GenomicRanges, rtracklayer, data.table, ape —
# all already installed, no new packages needed).
# Usage: micromamba run -n R Rscript scMultiOmics/2_TreeAlign/0_prepare_scdna_cn.R

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(rtracklayer)
  library(ape)
})

# ── Configuration ─────────────────────────────────────────────────────────────

MODE <- Sys.getenv("TREEALIGN_MODE", unset = "test")  # "test" or "full"
stopifnot(MODE %in% c("test", "full"))

TEST_REGION <- "R1"          # test mode: restrict scDNA cells (by clone_label prefix) to this region

# NOTE: test mode does NOT restrict to a single chromosome (an earlier version
# tried chr21-only, but TreeAlign's CloneAlignClone drops any gene whose
# per-clone consensus copy number is identical across all clones being
# compared — chr21 turned out to be completely CN-invariant across R1's own
# four clones, i.e. every gene failed that filter and the run errored with
# "No valid genes or snps exist in the matrix after filtering". Genome-wide
# coverage is needed so real inter-clone CNAs (wherever they are) are
# actually included; the region-based cell/clone subsetting below is what
# keeps test mode fast, not a gene-count restriction.

CN_PROFILES_TSV <- "/srv/home/aste0033/projects/MPNST/Haixi/scDNA/MEDICC2_output/all_10X_DLP/minseg_2.5/MPNST_all_single_cells_2.5_ds_6039_wgd_seed100/MPNST_all_single_cells_2.5_ds_6039_seed100_final_cn_profiles.tsv"
TREE_NEWICK      <- "/srv/home/aste0033/projects/MPNST/Haixi/scDNA/MEDICC2_output/all_10X_DLP/minseg_2.5/MPNST_all_single_cells_2.5_ds_6039_wgd_seed100/MPNST_all_single_cells_2.5_ds_6039_seed100_final_tree.new"
METADATA_CELLS   <- "scDNA/2_clustering/metadata_cells.tsv"
GTF_FILE         <- "/mnt/iribhm/genomes/hg38-gatk/gencode.v27.primary_assembly.annotation.gtf"

OUT_DIR <- file.path("/srv/home/aste0033/projects/MPNST/scMultiOmics/TreeAlign/inputs", MODE)
RESULTS_DIR <- "scMultiOmics/2_TreeAlign/results"

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)

cat(sprintf("[0_prepare_scdna_cn.R] MODE = %s, output -> %s\n", MODE, OUT_DIR))

# ── 1. Gene model: GENCODE v27, protein-coding, autosomes only ───────────────

cat("Loading GENCODE v27 GTF...\n")
gtf <- rtracklayer::import(GTF_FILE, feature.type = "gene")
genes_gr <- gtf[gtf$gene_type == "protein_coding" & seqnames(gtf) %in% paste0("chr", 1:22)]
# Same genome-wide gene set in both modes — see NOTE above on why test mode
# does not also restrict to a subset of chromosomes.

# GENCODE occasionally reuses a gene_name across distinct gene_ids (paralogs /
# readthrough transcripts); make.unique keeps every entry but disambiguates
# the column names we'll use downstream.
mcols(genes_gr)$gene_name <- make.unique(genes_gr$gene_name)
names(genes_gr) <- genes_gr$gene_name

cat(sprintf("Gene model: %d protein-coding genes on %s\n",
            length(genes_gr), paste(unique(as.character(seqnames(genes_gr))), collapse = ", ")))

# ── 2. Cell -> clone_label lookup, restricted to the 18 tumour clones ────────

metadata <- fread(METADATA_CELLS)
metadata <- metadata[!clone_label %in% c("healthy", "removed", "unassigned")]

# sample_id in the MEDICC2 output is "<region>_<technology>_<raw_id>", where
# raw_id keeps the 10x "-1" GEM-group suffix or the DLP+ ".bam" suffix. This
# is the SAME raw format used as tree tip labels, so we deliberately do not
# strip it any further than needed to join against metadata_cells.tsv (whose
# `barcode` column is already fully stripped of both).
parse_sample_id <- function(sample_id) {
  m <- regmatches(sample_id, regexec("^([^_]+)_(10X|DLP)_(.*)$", sample_id))
  region     <- vapply(m, `[`, character(1), 2)
  technology <- vapply(m, `[`, character(1), 3)
  raw        <- vapply(m, `[`, character(1), 4)
  barcode <- ifelse(technology == "10X", sub("-1$", "", raw), sub("\\.bam$", "", raw))
  data.table(sample_id = sample_id, region = region, technology = technology, barcode = barcode)
}

# ── 3. Load the all-cohort CN profile, keep only real single cells ──────────

cat("Loading all-cohort MEDICC2 CN profile (this is a large file, may take a minute)...\n")
cn <- fread(CN_PROFILES_TSV)
cn <- cn[grepl("_10X_|_DLP_", sample_id)]  # drop ancestral "internal_N" + "diploid" rows
cat(sprintf("Real single-cell segment rows: %d (%d distinct cells)\n",
            nrow(cn), uniqueN(cn$sample_id)))

cell_ids <- unique(cn$sample_id)
cell_lookup <- parse_sample_id(cell_ids)
cell_lookup <- merge(cell_lookup, metadata[, .(barcode, region, technology, clone_label)],
                      by = c("barcode", "region", "technology"), all.x = TRUE)

n_unmatched <- sum(is.na(cell_lookup$clone_label))
if (n_unmatched > 0) {
  stop(sprintf("%d cells in the CN profile could not be matched to a clone_label in %s — investigate before proceeding.",
               n_unmatched, METADATA_CELLS))
}

if (MODE == "test") {
  # Deliberately filtered by clone_label prefix, NOT the `region` column:
  # region (physical sampling site) and clone_label (genomic k-means cluster)
  # are different axes, and cross-region clone "leakage" is real and
  # documented (scMultiOmics/README.md) — filtering on `region` would pull in
  # a long tail of cells whose clone_label belongs to other regions entirely.
  cell_lookup <- cell_lookup[grepl(paste0("^", TEST_REGION, "_"), clone_label)]
}

cat(sprintf("Cells in scope: %d, spanning %d clones (%s)\n",
            nrow(cell_lookup), uniqueN(cell_lookup$clone_label),
            paste(sort(unique(cell_lookup$clone_label)), collapse = ", ")))

cn <- cn[sample_id %in% cell_lookup$sample_id]

# ── 4. Segment -> gene overlap (single vectorised findOverlaps, not a loop) ──

cat("Mapping segments to genes (GenomicRanges::findOverlaps)...\n")
seg_gr <- GRanges(seqnames = cn$chrom, ranges = IRanges(start = cn$start, end = cn$end))
mcols(seg_gr)$sample_id <- cn$sample_id
mcols(seg_gr)$total_cn  <- cn$cn_a + cn$cn_b

hits <- findOverlaps(seg_gr, genes_gr)
overlap_width <- width(pintersect(seg_gr[queryHits(hits)], genes_gr[subjectHits(hits)]))

overlap_dt <- data.table(
  sample_id = mcols(seg_gr)$sample_id[queryHits(hits)],
  gene      = names(genes_gr)[subjectHits(hits)],
  total_cn  = mcols(seg_gr)$total_cn[queryHits(hits)],
  overlap_width = overlap_width
)

# A gene rarely spans >1 segment (segments are Mb-scale, genes are not); when
# it does, keep the segment with the larger overlap.
setorder(overlap_dt, gene, sample_id, -overlap_width)
overlap_dt <- unique(overlap_dt, by = c("gene", "sample_id"))

# ── 5. Long -> wide (gene x cell), drop genes with any missing cell ─────────

cnv_wide <- dcast(overlap_dt, gene ~ sample_id, value.var = "total_cn")
n_genes_before <- nrow(cnv_wide)
complete_rows <- complete.cases(cnv_wide)
cnv_wide <- cnv_wide[complete_rows]
cat(sprintf("Genes with complete coverage across all %d cells: %d / %d (dropped %d with a coverage gap)\n",
            length(cell_ids <- unique(cn$sample_id)), nrow(cnv_wide), n_genes_before,
            n_genes_before - nrow(cnv_wide)))

# ── 6. Write outputs ──────────────────────────────────────────────────────────

fwrite(cnv_wide, file.path(OUT_DIR, "gene_by_cell_cnv.csv.gz"))

clone_table <- cell_lookup[, .(cell_id = sample_id, clone_id = clone_label)]
fwrite(clone_table, file.path(OUT_DIR, "cell_clone_labels.csv"))

gene_positions <- as.data.table(genes_gr)[, .(gene = gene_name, chrom = seqnames, start, end)]
gene_positions <- gene_positions[gene %in% cnv_wide$gene]
fwrite(gene_positions, file.path(OUT_DIR, "gene_positions.tsv"), sep = "\t")

cat("Preparing matching tree...\n")
tree <- ape::read.tree(TREE_NEWICK)
# Always prune to exactly the cell set in scope, in both modes: the tree has
# 6,040 tips, not 6,039 — MEDICC2 includes a "diploid" reference genome as an
# actual outgroup tip (root), which correctly never appears in the CN profile
# (it's not a real cell). Assuming "full mode needs no pruning" was wrong;
# test mode's drop.tip happened to remove "diploid" for free as part of
# dropping everything outside the region, which is why this only surfaced
# once full mode's stricter setequal() check ran.
drop_tips <- setdiff(tree$tip.label, cell_lookup$sample_id)
tree <- ape::drop.tip(tree, drop_tips)
cat(sprintf("Pruned tree: %d tips retained (dropped %d: %s)\n",
            length(tree$tip.label), length(drop_tips), paste(head(drop_tips, 5), collapse = ", ")))
stopifnot(setequal(tree$tip.label, cell_lookup$sample_id))
ape::write.tree(tree, file.path(OUT_DIR, "tree.newick"))

# ── 7. Small QC summary kept in the repo (not the big matrices) ─────────────

qc <- data.table(
  mode = MODE,
  n_cells = nrow(cell_lookup),
  n_clones = uniqueN(cell_lookup$clone_label),
  n_genes_before_coverage_filter = n_genes_before,
  n_genes_final = nrow(cnv_wide),
  n_tree_tips = length(tree$tip.label)
)
fwrite(qc, file.path(RESULTS_DIR, sprintf("0_prepare_scdna_cn_qc_%s.tsv", MODE)), sep = "\t")

cat("Done.\n")
print(qc)
