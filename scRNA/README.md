# scRNA — cell-cycle analysis of the MPNST single-nucleus RNA-seq data

This stage estimates S-phase fractions from single-nucleus RNA-seq (snRNA-seq) data
for the same patient and six spatial regions (P, R1–R5) covered by `scDNA/`, and asks
whether tumour and healthy cells differ in how much they proliferate.

**Headline result:** yes. Pooling all six regions, malignant cells are in S-phase far
more often than the normal (immune/stromal) cells embedded in the same samples:

| | n cells | S-phase fraction |
|---|---|---|
| Healthy | 3,450 | **12.9%** |
| Malignant | 34,266 | **27.2%** |

---

## Where the data actually lives (read this before writing new scripts here)

Two different sets of Seurat objects exist for this patient's scRNA data, at very
different stages of processing. **Only the second one is used by the scripts in this
directory.**

### 1. Bare per-region objects (not used here)

`/srv/home/aste0033/projects/MPNST/Haixi/scRNA/seurat/MPNST_{P,R1,R2,R3,R4,R5}.rds`

Six separate Seurat v3.1.2 objects (need `UpdateSeuratObject()` before use with
current Seurat), one per region, ~700 MB–1.9 GB each. QC'd, log-normalized, scaled,
and clustered (`seurat_clusters`, resolution 0.5) — but with **no cell-cycle score
and no tumour/healthy or cell-type annotation**. Cell counts are close to but not
identical to the Zenodo object below (e.g. region P: 3,984 cells here vs 4,019 in the
Zenodo object) — a slightly different/earlier QC pass on the same underlying cells.
Gene symbols (not Ensembl IDs) are used as feature names.

### 2. Published, fully annotated object (used by `1_cellcycle.R` / `1_cellcycle2.R`)

`/srv/home/aste0033/projects/MPNST/zenodo_phase_1/MPNST-Zenodo/data/snRNA/seurat_objects/MPNST_C_updated.rds`

This is the actual object behind the published analysis — see **Provenance** below.
All six regions merged into a single Seurat object (37,716 cells, `RNA` + `SCT`
assays), barcodes prefixed by region (`"P_AAACCCAAGCGACTTT"`, `"R1_..."`, ...).
Unlike the bare objects, it already carries:

| `meta.data` column | What it is |
|---|---|
| `S.Score`, `G2M.Score`, `Phase` | Seurat `CellCycleScoring()` output, already computed by the original authors |
| `tumor_normal` | Per-cell **Tumor** / **Normal** call, from CNV-based clustering (not a marker-gene heuristic) |
| `cell_type` | Finer label: malignant subclones (`Malignant_R1_1`, `Malignant_R5`, ...) for Tumor cells; `Macrophage` / `T_cell` / `Endothelial` / `Skeletal_Muscle` for Normal cells |
| `inferCNV_clusters` | The underlying CNV-based clustering (k=24) that `tumor_normal`/`cell_type` were derived from |
| `scenic_clusters`, `scRNA_haplo_group` | Not used here; SCENIC regulon clusters and SNP-haplotype-based doublet/genotype groups respectively |

We chose to build on this object rather than the bare per-region ones because it has
an actual CNV-validated tumour/normal ground truth. The alternative — inventing a
marker-gene heuristic to split tumour vs healthy ourselves on the bare objects — was
considered and rejected once this object was found, since it would have been strictly
weaker evidence for the same question.

### Provenance

This is snRNA-seq (single-**nucleus**, not single-cell) data from:

> Van Loo lab. *Chromosomal instability shapes spatial and temporal phenotypic
> diversity in a malignant peripheral nerve sheath tumour.*
> Data: [Zenodo record 19653314](https://zenodo.org/records/19653314) (`MPNST-Zenodo.zip`, 6.9 GB).
> Figure code: [github.com/VanLoo-lab/mpnst_phase_1](https://github.com/VanLoo-lab/mpnst_phase_1).

The Zenodo archive was downloaded and unzipped to
`/srv/home/aste0033/projects/MPNST/zenodo_phase_1/MPNST-Zenodo/`. It contains far more
than the snRNA cell-cycle data used here — `data/scDNA/`, `data/spRNA/` (spatial
Visium + CARD deconvolution), `data/GnT/`, `data/LCM/` (laser-capture microdissection
genotyping), `data/bulk/` (Battenberg/DPClust/CCF) — none of which is used by this
stage. The `results/figure_3/figure_3C.R` and `figure_3D.R` scripts in the GitHub repo
are what pointed us at `inferCNV_clusters` and `MPNST_C_updated.rds` in the first
place.

**snRNA vs scRNA matters here:** nuclear RNA under-samples some cytoplasmic and
late-M-phase transcripts relative to whole-cell RNA, which is a plausible source of
disagreement between any cell-cycle score computed here and results from a
true single-cell (not single-nucleus) protocol.

---

## Scripts

| Script | Purpose |
|--------|---------|
| `1_cellcycle.R` | Reports S-phase fractions per region and per cluster using the **published** `Phase` / `tumor_normal` / `cell_type` columns. Fast, no new computation. |
| `1_cellcycle2.R` | Recomputes cell-cycle scores **from scratch** with `Seurat::CellCycleScoring()`, validates against the published `Phase` calls, and quantifies how much the result shifts under different gene-list/assay choices. |

Both must be run from the repository root, via the container that has Seurat
(nothing scRNA-related is installed in any micromamba env):

```bash
singularity exec /mnt/iribhm/software/singularity/single-cell.sif \
  Rscript scRNA/1_cellcycle.R
singularity exec /mnt/iribhm/software/singularity/single-cell.sif \
  Rscript scRNA/1_cellcycle2.R
```

### Outputs (`results/`)

| File | From | Description |
|------|------|--------------|
| `sphase_summary.tsv` | `1_cellcycle.R` | `region \| group \| class \| n_cells \| n_sphase \| sphase_fraction`, one row per region × (healthy or malignant cluster) |
| `sphase_fractions.png` | `1_cellcycle.R` | One panel per region; bars = healthy + every malignant cluster with cells present in that region (including cross-region leakage), coloured by class |
| `sphase_fractions_region_native.png` | `1_cellcycle.R` | Same, but each panel keeps only "healthy" + clusters *native* to that region (e.g. `P` keeps only `Malignant_P`); cross-region rows and clusters with no single native region (`Malignant_Ribo`, `Malignant_G2MS_*`) are dropped, mirroring the region-native filter in `scDNA/4_sprinter/persample/plot_sphase_persample.R`. The cleaner headline view. |
| `cellcycle_method_comparison.tsv` | `1_cellcycle2.R` | S-phase fraction (overall, Normal, Tumor) for the published calls and three of our own re-runs |
| `cellcycle_method_sensitivity.png` | `1_cellcycle2.R` | Bar chart of the same comparison |

### A caveat when reading `sphase_fractions.png`

A handful of `cell_type` clusters (`Malignant_G2MS_R1R5`, `Malignant_G2MS_R2R3`,
`Malignant_G2MS_R4`) are **defined by having high G2M/S scores** — checking, their
G2M fraction is 54–74% by construction. Their S-phase-fraction bars are not an
independent proliferation measurement; they are close to tautological. This is a
general risk in scRNA/snRNA clustering: cell-cycle genes carry a strong transcriptional
signal and can dominate unsupervised clustering unless explicitly regressed out (see
`CellCycleScoring()` + `ScaleData(vars.to.regress = ...)` in the Seurat vignette). It's
unclear whether the original clustering regressed cell-cycle out before defining these
particular clusters; treat their bars as descriptive, not as evidence of a real
proliferation difference.

Also, most `cell_type` cluster labels carry a region suffix (e.g. `Malignant_R1_1`)
but are not perfectly confined to that region — e.g. `Malignant_R5` has 6,341 cells in
R5 but also 162 in R1. `sphase_fractions.png` plots each cluster in every region where
it actually has cells, faithfully, rather than filtering this out — the paper's own
title is about *spatial* phenotypic diversity, so this cross-region mixing may be a
real finding rather than noise. Bars from n < 15 cells are shown at 40% opacity as a
reminder that they're not reliable regardless.

---

## Cell-cycle scoring: methods, and why Seurat's is used here

### The standard method: `Seurat::CellCycleScoring()`

This is the default approach in the Seurat ecosystem, and what the original authors
used (`1_cellcycle2.R` reproduces their exact `Phase` calls at 100% concordance using
the `RNA` assay's log-normalized data and the `cc.genes.updated.2019` gene list — see
below). Mechanically: it is `AddModuleScore()` run twice, once for a curated list of
S-phase marker genes and once for G2M-phase markers. `AddModuleScore()` bins all genes
by average expression, and for each gene in the list subtracts the mean expression of
a random control gene set drawn from the *same expression bin* — this corrects for the
fact that raw expression/dropout rate, not cell-cycle biology, would otherwise dominate
a naive mean. Each cell gets an `S.Score` and a `G2M.Score`; `Phase` is `"G2M"` if
G2M.Score > S.Score and > 0, `"S"` if S.Score wins and is > 0, else `"G1"`.

### Alternatives considered

| Method | Pros | Cons | Used here? |
|---|---|---|---|
| **Seurat `CellCycleScoring`** | Native, fast, zero extra deps, matches the published pipeline | Discrete G1/S/G2M only; can't separate true quiescence (G0) from G1; gene list curated from cell lines, may transfer imperfectly | **Yes — primary method** |
| `scran::cyclone` (Bioconductor) | Pairs-based, different normalization assumptions | Much slower; needs `SingleCellExperiment`; still discrete phases | No — `scran`/`scater`/`SingleCellExperiment` are installed in `single-cell.sif` if ever needed |
| `tricycle` (Bioconductor) | Continuous cell-cycle "angle" instead of 3 bins; can reveal subtle proliferation gradients | Newer, less battle-tested; awkward circular-coordinate interpretation; not installed | No |
| `Revelio` / `ccAFv2` | `ccAFv2` distinguishes true G0 from G1 — directly relevant to a tumour-vs-normal quiescence question | Niche; `ccAFv2` needs a Python/keras model; more setup for marginal gain given the result already separates cleanly | No |

### Parameter tuning — what `1_cellcycle2.R` demonstrates

| Parameter | Choice used | Effect of changing it (measured in `1_cellcycle2.R`) |
|---|---|---|
| Gene list | `cc.genes.updated.2019` (current Seurat recommendation) | Reverting to the original 2015 `cc.genes` drops concordance with the published calls to 92.4%. Partly mechanical: `MLF1IP`, `FAM64A`, `HN1` (old symbols in the 2015 list) aren't in this object at all — they were renamed to `CENPU`, `PIMREG`, `JPT1` and only the updated list uses the current symbols. **Gene-symbol/nomenclature drift is a real, silent failure mode** — genes that "aren't found" are dropped with only a warning, not an error. |
| Assay / slot | `RNA` assay, `data` slot (log-normalized counts) | Must be log-normalized data, never raw counts or `scale.data`. Switching to the `SCT` assay's `data` slot drops concordance to 93.3% and compresses the Normal/Tumor gap (9.3% vs 24.8%, down from 12.9% vs 27.2%). SCTransform's variance-stabilizing correction interacts with `AddModuleScore()`'s own binning in a way plain `RNA` data doesn't. |
| `nbin`/`ctrl` in `AddModuleScore` | Seurat defaults | Not an issue at this scale (thousands of cells per region); can misbehave with very small groups (tens of cells) — not encountered here. |

The core biological result — malignant cells cycle roughly twice as often as
normal cells — holds across every variant tested, even though the exact percentages
move by several points. That robustness-despite-parameter-sensitivity is itself worth
knowing: don't over-trust the second decimal place of any single S-phase fraction, but
do trust a >2-fold, cross-method-consistent difference.

### General things to watch for when doing cell-cycle analysis on scRNA/snRNA data

- **Normalization first.** Cell-cycle scoring needs log-normalized expression, not
  raw counts — the scores are differences of means on a log scale.
- **G0 vs G1 is invisible to this method.** Seurat calls everything that isn't
  confidently S or G2M "G1", which conflates truly quiescent cells with cycling cells
  that just happen to be caught between S and G2M. This matters most for exactly the
  question asked here (do healthy cells rest more than tumour cells?) — a `ccAFv2`-style
  method that separates G0 would give a sharper answer.
- **Cell-cycle signal can dominate clustering.** If clusters are defined without
  regressing out cell-cycle scores first, "clusters" can end up being proliferation
  states rather than cell identities (see the `Malignant_G2MS_*` caveat above).
- **Small groups are noisy.** A handful of cells gives a S-phase fraction that swings
  wildly (0% or 100%) on one or two cells. Always carry `n_cells` alongside any
  fraction and discount/flag anything under a few dozen cells.
- **Gene symbol drift is silent.** `CellCycleScoring()` warns but doesn't error when
  list genes are missing from the object — always check how many of the 43+54 genes
  were actually found (see `1_cellcycle2.R`).
- **snRNA ≠ scRNA.** Nuclear preparations under-sample some cytoplasmic and late-M
  transcripts; don't assume a score computed here transfers unchanged to whole-cell
  data from the same tumour.
