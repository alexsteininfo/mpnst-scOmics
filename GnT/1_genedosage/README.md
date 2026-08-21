# GnT/1_genedosage — gene-dosage effect (DNA copy number vs RNA expression)

Does a gene's DNA copy number in a cell predict its RNA expression in that
same cell? GnT is the one technology in this repo where that question can be
asked directly, cell-by-cell, because DNA and RNA come from the same lysed
cell (see [`../README.md`](../README.md)). This stage reproduces the
gene-dosage panel from the manuscript (published as **Figure 2E**; the
`mpnst_phase_1` code repo's own internal draft numbering calls the same
analysis `figure_3/figure_3E.R` — content matches down to the `min_cells=10`
/ `exp_cutoff=10` thresholds, it's the same figure under a different draft
number).

## The two axes, with formulas

Both quantities are Seurat **SCTransform-corrected UMI counts**
(`SCT` assay, `counts` slot) — SCTransform fits a regularized
negative-binomial regression of each gene's raw UMI count against
sequencing depth, then re-expresses every cell's count as what it would be
at a common reference depth (the median depth across cells). That makes it a
UMI count **normalized for sequencing depth**, not log-transformed — hence
"normalized UMI count" here, rather than the "raw UMI count (SCT)" label
used in the original script (a shorthand for "count-scale, as opposed to the
log1p `data` slot", not "unprocessed").

**Mean observed** — for gene *g* at DNA copy-number state *n*, the plain
average SCT count across the GnT cells actually carrying *n* copies of *g*
(kept only if `n_cells >= MIN_CELLS`):

```
mean_observed(g, n) = mean( SCT_count(g, c) for c in cells where CN(g, c) == n )
```

**Mean expected** — a *prediction*, not a measurement. Take a reference
copy-number state `n_low` of the same gene with its own observed mean. Under
strict linear dosage (each extra copy contributes equally, no compensation),
the expected mean at a higher state `n_high` is that reference rescaled by
the copy-number ratio:

```
mean_expected(g, n_high | n_low) = mean_observed(g, n_low) * (n_high / n_low)
```

Every row in `gene_dosage_pairs.tsv` is one `(gene, n_high, n_low)` triple
with `n_high > n_low`. The figure plots `mean_expected_high` (x) against the
truly observed `mean_SCT_count_high` (y), both log10. Points on the red
y = x line have expression that scales exactly with copy number; below it
= dosage compensation/attenuation; above it = super-linear amplification.

## Script

| Script | Purpose |
|--------|---------|
| `gene_dosage_plot.R` | Matches the 130 GnT cells with usable DNA+RNA (DNA barcode ↔ RNA barcode via the EGAF ±1 accession adjacency, see `../README.md`), restricts to genes with median SCT count ≥10 across all 300 GnT RNA cells, computes `mean_observed`/`mean_expected` per copy-number pair (≥10 cells per group), and plots. |

Needs Seurat, not installed in any micromamba env:

```bash
singularity exec /mnt/iribhm/software/singularity/single-cell.sif \
  Rscript GnT/1_genedosage/gene_dosage_plot.R
```

Reads from the staged GnT copy at
`/srv/home/aste0033/projects/MPNST/Haixi/GnT/zenodo/` (see the top-level
[`GnT/README.md`](../README.md) for what's there and how it was produced).
Must be run from the repository root.

### Outputs (`results/`)

| File | Description |
|------|-------------|
| `gene_dosage_expected_vs_observed.png` | The log-log scatter described above. |
| `gene_dosage_pairs.tsv` | Every `(gene, CN_high, CN_low)` pair behind it: `n_cells_*`, `median_SCT_count_*`, `mean_SCT_count_*`, `mean_expected_high`, `log10_abs_diff`, and `highlight` (the 25 largest `\|expected - observed\|` discrepancies on log scale). |

**Note on `highlight`:** the original manuscript script computes this same
top-25-outlier flag but maps both `TRUE` and `FALSE` to black in the actual
plot — i.e. the published figure does not visually distinguish the
outliers, only the underlying analysis code identifies them. That's
reproduced faithfully here; use the `highlight` column in the TSV (or the
`gene_name` values `arrange(desc(log10_abs_diff))`-sorted at the top) if you
want to actually see which genes they are.

### Current run

753 genes pass the expression filter, yielding 1,473 dosage pairs (≥10
cells per copy-number group).
