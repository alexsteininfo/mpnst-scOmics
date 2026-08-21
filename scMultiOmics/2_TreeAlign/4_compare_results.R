#!/usr/bin/env Rscript
# Cross-tabulate TreeAlign's two output modes against each other and against
# the published scRNA cell_type clusters, to actually close the "independent
# clusterings, not cross-matched" gap flagged in scMultiOmics/README.md.
#
# Input:
#   - <OUT_DIR>/clonelabel/clone_assign_df.csv  (2_run_treealign_clonelabel.py)
#   - <OUT_DIR>/tree/clone_assign_df.csv        (3_run_treealign_tree.py)
#   - <IN_DIR>/scrna_cell_metadata.csv          (1_prepare_scrna_expr.R)
#
# Output (results/):
#   - clonelabel_vs_celltype.tsv, treemode_vs_celltype.tsv
#   - clonelabel_vs_treemode.tsv, clonelabel_vs_region.tsv
#   - clonelabel_vs_celltype.png (heatmap)
#
# Requires: R micromamba env (data.table, ggplot2).
# Usage: micromamba run -n R Rscript scMultiOmics/2_TreeAlign/4_compare_results.R

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# ── Configuration ─────────────────────────────────────────────────────────────

MODE <- Sys.getenv("TREEALIGN_MODE", unset = "test")
stopifnot(MODE %in% c("test", "full"))

IN_DIR  <- file.path("/srv/home/aste0033/projects/MPNST/scMultiOmics/TreeAlign/inputs", MODE)
OUT_DIR <- file.path("/srv/home/aste0033/projects/MPNST/scMultiOmics/TreeAlign/outputs", MODE)
RESULTS_DIR <- "scMultiOmics/2_TreeAlign/results"
dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)

cat(sprintf("[4_compare_results.R] MODE = %s\n", MODE))

# ── Load ──────────────────────────────────────────────────────────────────────

clonelabel <- fread(file.path(OUT_DIR, "clonelabel", "clone_assign_df.csv"))
setnames(clonelabel, "clone_id", "clone_id_labelmode")

treemode <- fread(file.path(OUT_DIR, "tree", "clone_assign_df.csv"))
setnames(treemode, "clone_id", "clone_id_treemode")

meta <- fread(file.path(IN_DIR, "scrna_cell_metadata.csv"))

merged <- merge(meta, clonelabel, by = "cell_id", all.x = TRUE)
merged <- merge(merged, treemode, by = "cell_id", all.x = TRUE)

cat(sprintf("scRNA cells: %d total, %d assigned by clone-label mode, %d assigned by tree mode\n",
            nrow(merged), sum(!is.na(merged$clone_id_labelmode)), sum(!is.na(merged$clone_id_treemode))))

write_crosstab <- function(dt, col1, col2, fname) {
  tab <- dt[, .N, by = c(col1, col2)]
  setorderv(tab, c(col1, col2))
  fwrite(tab, file.path(RESULTS_DIR, fname), sep = "\t")
  tab
}

ct_label_celltype <- write_crosstab(merged, "clone_id_labelmode", "cell_type",
                                     sprintf("clonelabel_vs_celltype_%s.tsv", MODE))
invisible(write_crosstab(merged, "clone_id_treemode", "cell_type",
                          sprintf("treemode_vs_celltype_%s.tsv", MODE)))
invisible(write_crosstab(merged, "clone_id_labelmode", "clone_id_treemode",
                          sprintf("clonelabel_vs_treemode_%s.tsv", MODE)))
invisible(write_crosstab(merged, "clone_id_labelmode", "region",
                          sprintf("clonelabel_vs_region_%s.tsv", MODE)))

# ── Headline figure: assigned scDNA clone (clone-label mode) x published scRNA cell_type ──

p <- ggplot(ct_label_celltype, aes(x = clone_id_labelmode, y = cell_type, fill = N)) +
  geom_tile() +
  scale_fill_viridis_c(name = "n cells") +
  labs(title = sprintf("TreeAlign clone-label assignment vs. published scRNA cell_type (%s)", MODE),
       x = "TreeAlign-assigned scDNA clone", y = "Published scRNA cell_type") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(RESULTS_DIR, sprintf("clonelabel_vs_celltype_%s.png", MODE)), p, width = 9, height = 6, dpi = 150)

cat("Done. Wrote cross-tabs and heatmap to", RESULTS_DIR, "\n")
