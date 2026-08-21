# scMultiOmics — cross-technology comparison

Combines S-phase fraction estimates from all three technologies covering this
patient's six regions (P, R1–R5) into one figure: **10X** and **DLP+** (scDNA,
from `scDNA/4_sprinter/persample/`) and **scRNA** (from `scRNA/1_cellcycle.R`).

## Script

| Script | Purpose |
|--------|---------|
| `sphase_comparison.R` | Reads both stages' summary TSVs, aligns them onto a common per-region clone axis, and plots three dodged bars (10X / DLP / scRNA) per clone, analogous to `scDNA/4_sprinter/results/sphase_fractions_persample.png` but with the third technology added. |

No Seurat/singularity needed — it only reads two already-computed TSVs:

```bash
micromamba run -n R Rscript scMultiOmics/sphase_comparison.R
```

### Outputs (`results/`)

| File | Description |
|------|--------------|
| `sphase_comparison.tsv` | `region \| technology \| clone_label \| n_cells \| n_sphase \| sphase_fraction`, completed to all 3 technologies × clone × region combinations (missing = 0 cells) |
| `sphase_comparison.png` | The three-technology bar chart |

## The clone-label caveat (read this before drawing conclusions from a single bar)

**scDNA and scRNA clones are two independent clusterings, not a joint one.**
`clone_label` values like `R1_1`, `R4_2` in the scDNA summary come from k-means
on CNA profiles (`scDNA/2_clustering/`). scRNA's clusters (`Malignant_R1_1`,
`Malignant_R4_2`, ... in `scRNA/results/sphase_summary.tsv`) come from a
transcriptome/inferCNV-based clustering the original authors ran (see
`scRNA/README.md`). **No cross-technology cell-level matching was performed.**
Where a label happens to coincide between technologies (this only happens for
R4, which coincidentally has 3 subclusters in both), the shared label does
**not** mean the same underlying cells are being compared — it's the same
region and the same numbered slot in two unrelated partitions.

scRNA does not subdivide malignant cells into multiple clusters at all for P,
R2, R3, and R5 (one `Malignant_<region>` cluster covers the whole region's
tumour compartment there). Those are plotted as their own `<region>_all` bar
rather than being folded into any specific scDNA clone, so at least those bars
are honestly labelled as "whole malignant compartment" rather than implying a
false match to one particular scDNA clone.

**What you *can* trust from this figure:**
- The **healthy** bars are a clean three-way comparison — "healthy" means the
  same thing (non-malignant cells) in all three technologies.
- The **overall shape** — malignant S-phase fraction consistently higher than
  healthy, across nearly every region and technology — is a robust,
  cross-technology-consistent signal.

**What you should not overinterpret:**
- Any claim that "clone R4_1 is more proliferative by DLP+ than by scRNA" —
  that compares two different partitions of the same region's cells, not the
  same clone measured twice.
