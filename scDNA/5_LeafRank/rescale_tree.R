# Purpose: Convert the tumor-clade tree from prepare_tree.R (branch lengths =
#          MEDICC2 event counts, not ultrametric) into an ultrametric
#          chronogram LeafRank requires (branch lengths proportional to time).
# Input:   <group_output_dir>/<group_name>_tumor_tree.rds (from prepare_tree.R)
# Output:  <group_output_dir>/<group_name>_ultrametric_tree.rds (ape "chronos"/
#          "phylo" object, root-to-tip distance normalized to 1)
# Requires: micromamba run -n R Rscript scDNA/5_LeafRank/rescale_tree.R
#           (invoke from the repo root)
#
# NOTE ON METHOD CHOICE (decided 2026-08-03, per user direction): this cancer's
# WGD is treated as fully truncal throughout -- any MEDICC2 "subclonal WGD"
# calls (e.g. the per-cluster 0-2-flagged-node pattern found during initial
# research) are considered likely inference artifacts, not real biology. This
# means every tip in every tree is WGD+, so LeafRank's two-rate WGD-aware
# rescaling (get_ultrametric_prepared() -> MATLAB modified_chronos ->
# get_ultrametric(), requiring tips_WGD to have both 0s and 1s to estimate two
# separate rates) is unnecessary here -- it would be degenerate/uninformative
# on a uniformly-WGD+ vector. Instead we use the paper's own simpler path for
# truncal-WGD tumors: plain constant-rate `ape::chronos(tree, model = "clock")`
# (a strict molecular clock, no MATLAB required).
#
# KEEPING THE WGD-AWARE PATH FOR LATER: if a future tumor/cohort genuinely has
# subclonal WGD, do NOT reuse this script as-is -- implement the two-rate path
# documented in scDNA/5_LeafRank/todo.md step 3 instead (get_ultrametric_prepared()
# + MATLAB's modified_chronos.m + get_ultrametric()). MATLAB module availability
# was tested 2026-08-03 (module loads, but license checkout currently fails --
# see todo.md) and is NOT needed for the path implemented here.

source("scDNA/5_LeafRank/config.R")

suppressPackageStartupMessages(library(ape))

log_msg <- function(fmt, ...) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), sprintf(fmt, ...)))
}

tumor_tree_file <- file.path(group_output_dir, paste0(group_name, "_tumor_tree.rds"))
if (!file.exists(tumor_tree_file)) {
  stop("Tumor tree not found: ", tumor_tree_file, " -- run prepare_tree.R first.")
}

tree_tumor <- readRDS(tumor_tree_file)
log_msg("Loaded tumor tree: %d tips, %d internal nodes, %d zero-length branches",
        ape::Ntip(tree_tumor), tree_tumor$Nnode, sum(tree_tumor$edge.length == 0))

log_msg("Fitting strict-clock chronogram via ape::chronos(model = 'clock') -- this took ~4.4 min for the 1571-tip R1 tree")
t0 <- Sys.time()
tree_ultrametric <- ape::chronos(tree_tumor, model = "clock")
elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
log_msg("chronos() finished in %.1f sec", elapsed)

# chronos() commonly emits a "false convergence" warning from its underlying
# nlminb optimizer on large trees -- this is usually benign (the fit still
# converges close enough to be usable) but worth confirming: check that the
# tree is genuinely ultrametric and root-to-tip distance is the expected 1.
stopifnot(ape::is.ultrametric(tree_ultrametric))
root_to_tip <- ape::node.depth.edgelength(tree_ultrametric)[1:ape::Ntip(tree_ultrametric)]
log_msg("Root-to-tip distance range: [%.6f, %.6f] (expect ~1, ~1)",
        min(root_to_tip), max(root_to_tip))

# ape::chronos() returns class c("chronos","phylo") with extra attributes;
# LeafRank's message-passing code only touches standard phylo slots
# (edge, edge.length, tip.label, Nnode), so this is safe to use directly.
class(tree_ultrametric) <- "phylo"

dir.create(group_output_dir, recursive = TRUE, showWarnings = FALSE)
out_file <- file.path(group_output_dir, paste0(group_name, "_ultrametric_tree.rds"))
saveRDS(tree_ultrametric, out_file)
log_msg("Saved ultrametric tree to: %s", out_file)
