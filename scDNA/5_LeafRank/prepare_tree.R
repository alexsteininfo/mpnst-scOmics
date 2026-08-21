# Purpose: Load a MEDICC2 tree for one group (region or cluster), reduce it to
#          a strictly-bifurcating tumor-only clade, and derive the per-tip
#          whole-genome-doubling (WGD) status vector LeafRank needs.
# Input:   <group_name>_final_tree.new, <group_name>_copynumber_events_df.tsv
#          under medicc2_output_base (see config.R).
# Output:  <group_output_dir>/<group_name>_tumor_tree.rds   (ape phylo, tumor
#          clade only, strictly bifurcating)
#          <group_output_dir>/<group_name>_tips_wgd.rds     (integer 0/1
#          vector, positionally aligned to tumor_tree$tip.label order)
#          (group_output_dir is under the project folder, not this repo --
#          see config.R for why.)
# Requires: micromamba run -n R Rscript scDNA/5_LeafRank/prepare_tree.R
#           (invoke from the repo root; sources config.R via a repo-relative path)

source("scDNA/5_LeafRank/config.R")

suppressPackageStartupMessages(library(ape))

log_msg <- function(fmt, ...) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), sprintf(fmt, ...)))
}

if (!file.exists(tree_file))      stop("Tree file not found: ", tree_file)
if (!file.exists(cn_events_file)) stop("CN events file not found: ", cn_events_file)

log_msg("Loading tree for group_mode=%s minseg=%s group_name=%s", group_mode, minseg, group_name)
log_msg("Tree file: %s", tree_file)
tree_raw <- ape::read.tree(tree_file)
log_msg("Loaded %d tips, %d internal nodes", ape::Ntip(tree_raw), tree_raw$Nnode)

### Extract tumor clade (drop diploid pseudo-tip / root trunk) ---------------
# Adapted from CNA-PhyloAnalysis/code/2_phylogenies/base_trees.R:44-61
# (tree_to_plot_data()) -- ported locally rather than sourced across repos.
# MEDICC2 tree structure: root has two children -- the tumor MRCA (long trunk
# branch) and a zero-length "diploid" pseudo-tip.
extract_tumor_clade <- function(tree) {
  if (!("diploid" %in% tree$tip.label)) return(tree)

  root_node   <- setdiff(tree$edge[, 1], tree$edge[, 2])
  root_kids   <- tree$edge[tree$edge[, 1] == root_node, 2]
  diploid_idx <- which(tree$tip.label == "diploid")
  mrca_idx    <- root_kids[root_kids != diploid_idx]

  if (length(mrca_idx) == 1L && mrca_idx > ape::Ntip(tree)) {
    ape::extract.clade(tree, mrca_idx)
  } else {
    ape::drop.tip(tree, "diploid")
  }
}

tree_tumor <- extract_tumor_clade(tree_raw)
log_msg("Tumor clade: %d tips, %d internal nodes (dropped diploid pseudo-tip/trunk)",
        ape::Ntip(tree_tumor), tree_tumor$Nnode)

### Ensure strictly bifurcating -----------------------------------------------
# LeafRank's message-passing code (R/message_passing.R: up_m_des, up_m_sib,
# get_sibling) hard-assumes exactly two children per internal node and will
# silently misbehave (not error) on a polytomy.
if (!ape::is.binary(tree_tumor)) {
  n_poly <- sum(table(tree_tumor$edge[, 1]) > 2)
  log_msg("Tree has %d polytomies -- resolving via ape::multi2di()", n_poly)
  tree_tumor <- ape::multi2di(tree_tumor)
} else {
  log_msg("Tree is already strictly bifurcating")
}

### Derive per-tip WGD status vector ------------------------------------------
# Confirmed mechanism (see todo.md step 2): <group_name>_copynumber_events_df.tsv
# has one row per actual WGD event (type == "wgd"), naming the node where it
# occurred (e.g. "internal_591"). All descendant tips of that node are WGD+;
# every other tip is WGD-. Reuses the same logic pattern as CNA-PhyloAnalysis's
# near_clonal_cnas.R has_wgd derivation.
cn_events <- read.delim(cn_events_file, stringsAsFactors = FALSE)
wgd_nodes <- unique(cn_events$sample_id[cn_events$type == "wgd"])
log_msg("Found %d WGD event(s) recorded on node(s): %s",
        length(wgd_nodes), paste(wgd_nodes, collapse = ", "))

tips_wgd <- setNames(integer(ape::Ntip(tree_tumor)), tree_tumor$tip.label)

for (node_label in wgd_nodes) {
  tip_match      <- which(tree_tumor$tip.label == node_label)
  internal_match <- which(tree_tumor$node.label == node_label)

  if (length(tip_match) == 1L) {
    tips_wgd[tip_match] <- 1L
  } else if (length(internal_match) == 1L) {
    node_id <- internal_match + ape::Ntip(tree_tumor)
    descendant_tips <- ape::extract.clade(tree_tumor, node_id)$tip.label
    tips_wgd[descendant_tips] <- 1L
  } else {
    warning("WGD event node '", node_label, "' not found in the tumor clade ",
            "(possibly pruned during clade extraction/multi2di, or on the ",
            "trunk/diploid side already removed) -- skipping this event.")
  }
}

log_msg("WGD+ tips: %d / %d (%.1f%%)",
        sum(tips_wgd), length(tips_wgd), 100 * mean(tips_wgd))

# LeafRank's get_ultrametric_prepared() indexes tips_WGD POSITIONALLY, not by
# name (confirmed from LeafRank's R/tree_utils.R source) -- strip names right
# before saving; order here is guaranteed to match tree_tumor$tip.label.
tips_wgd_vec <- unname(tips_wgd[tree_tumor$tip.label])

### Save -----------------------------------------------------------------
dir.create(group_output_dir, recursive = TRUE, showWarnings = FALSE)
tree_out <- file.path(group_output_dir, paste0(group_name, "_tumor_tree.rds"))
wgd_out  <- file.path(group_output_dir, paste0(group_name, "_tips_wgd.rds"))

saveRDS(tree_tumor, tree_out)
saveRDS(tips_wgd_vec, wgd_out)

log_msg("Saved tumor tree to: %s", tree_out)
log_msg("Saved tips_WGD vector to: %s", wgd_out)
