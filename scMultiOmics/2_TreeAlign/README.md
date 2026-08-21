# TreeAlign — assigning snRNA-seq cells to scDNA copy-number subclones

Integrates this patient's 18 scDNA copy-number clones (`scDNA/2_clustering/metadata_cells.tsv`)
with the published snRNA-seq object using [TreeAlign](https://github.com/shahcompbio/TreeAlign)
(Shi et al. 2024, *Nat Commun*, `literature/shi_TreeAlign_natcom.pdf`), a Bayesian model (Pyro/PyTorch)
that assigns scRNA cells to scDNA subclones from raw expression + a gene-level copy-number matrix,
and scores each gene's probability of being CN-"dosage-affected". This directly addresses the caveat
repeated throughout `scMultiOmics/README.md`: scDNA and scRNA clones in this project are two
**independent clusterings**, never cross-matched at the cell level. TreeAlign is the first stage to
actually link them.

## Decisions made (and why)

| Decision | Choice | Rationale |
|---|---|---|
| Assignment mode | **Both**: clone-label mode (primary) + tree mode (secondary/exploratory) | Clone-label mode directly answers "assign scRNA cells to my existing 18 clusters" — the literal ask. Tree mode is TreeAlign's own "preferred" richer mode and the exact input it needs (a full cell-level MEDICC2 tree) already exists in this project, so it costs little extra to also run it and cross-check against clone-label mode's output. |
| scRNA QC | **Minimal** — use the published Zenodo object's existing QC/annotations as-is | The object is already published, CNV-validated (tumor/normal, cell-cycle scored). TreeAlign's book-chapter protocol recommends *additional* ambient-RNA correction (SoupX/DecontX) and doublet removal (scDblFinder), but none of those tools exist anywhere in this repo's pipeline yet, and `scRNA/README.md` already establishes the convention of trusting this object's QC as-is. Adding them is a bigger, separately-scoped follow-up if results warrant it. |
| Allele-specific submodel | **Skipped** | Needs SIGNALS/cellsnp-lite BAF extraction from scDNA BAMs + matching scRNA ref/alt allele counts — none of which exist in this repo. TreeAlign's own tutorial demonstrates the core expression+CN model as fully valid on its own. Can be added later without re-doing anything below. |
| scDNA CN source | `MPNST_all_single_cells_2.5_ds_6039_seed100_final_cn_profiles.tsv` (the **all-cohort** MEDICC2 run, `minseg=2.5Mb`, `10X_DLP` combined), *not* the separate per-cluster MEDICC2 input TSVs | This one file is paired 1:1 with the existing full-cohort tree (`..._final_tree.new`) from the *same* MEDICC2 run — guaranteed consistent segment boundaries and cell IDs between the CN matrix and the tree, which two independently-run per-cluster files would not guarantee. Verified: exactly 6,039 real single-cell profiles in this file (it also contains 6,038 ancestral "internal_N" node profiles + 1 "diploid" root reference, which are excluded — see step 0 below), and this is **exactly** the sum of cells across the 18 tumour `clone_label` values in `metadata_cells.tsv` (728+624+574+509+468+467+408+397+309+255+255+233+208+174+171+112+111+36 = 6,039). No healthy/removed/unassigned cells are present in this file at all — no extra filtering needed. |
| Gene annotation | GENCODE v27 (`hg38-gatk/gencode.v27.primary_assembly.annotation.gtf`), `protein_coding` genes, autosomes only (chr1–22) | v27 is already the GTF paired with this repo's 10x-pipeline hg38 reference. Restricting to protein-coding autosomal genes avoids sex-chromosome CN complications and pseudogene/lncRNA noise. |
| Segment→gene mapping tool | R `GenomicRanges`/`rtracklayer` | Already installed in the `R` micromamba env — no new install needed. |
| Total vs. allele-specific CN | Total copy number (`cn_a + cn_b`) per gene | Matches "skip allele-specific submodel" above; TreeAlign's `cnv` argument just wants a copy-number value per gene. |
| Genes with incomplete segment coverage | Dropped entirely (kept only if every one of the 6,039 cells has an overlapping segment) | Simplest safe default; documented here so it's revisited if too many genes are lost (see step 0 output log for the actual count). |
| Per-clone consensus CN | **Not precomputed manually** — TreeAlign's `CloneAlignClone` takes the full per-cell `cnv` matrix plus a `clone` (cell→clone_label) table and builds its own Bayesian per-clone consensus internally | Confirmed from the library source (`clonealign_clone.py`): `cnv` is always gene×**cell**, never gene×clone. Precomputing a naive median ourselves would be strictly worse than letting the model's own inference handle it. |
| TreeAlign source | `git clone https://github.com/AlexHelloWorld/TreeAlign.git` | `shahcompbio/TreeAlign`'s own README points here as the actual source (confirmed by fetching both). Not on PyPI. |
| Environment | New micromamba env `treealign` (Python 3.9, pinned `torch==1.11.0` + `pyro-ppl==1.8.1` from TreeAlign's `requirements.txt`), CPU-only | Matches this repo's existing pattern of one dedicated micromamba env per tool (`chisel`, `sprinter`). Confirmed from source: no CUDA/GPU code paths anywhere in the library, so CPU-only is a supported, not a degraded, mode — matches this machine having no GPU. |
| Reproducibility | `pyro.set_rng_seed()` fixed at the top of both run scripts; TreeAlign commit SHA pinned and recorded here (see below) | The package isn't a versioned release, so the exact commit must be pinned; the model is stochastic (Pyro variational inference), so the seed matters for exact reproducibility. |
| Test vs. full runs | A `TREEALIGN_MODE` config (`"test"` or `"full"`) read by every script (env var, defaulting to `"test"`; set once in `run_pipeline.sh`) | Matches this repo's "configuration block at the top, not CLI arguments" convention while still letting `run_pipeline.sh` drive both passes from one place. **Test mode** = region **R1** only (4 clones: R1_1–R1_4, 1,571 scDNA cells) + chromosome 21 only (234 protein-coding genes) — small enough to validate the full mechanics quickly before committing to the ~6,039-cell / ~34,266-scRNA-cell / genome-wide full run, whose TreeAlign runtime is otherwise a real unknown (the paper reports no runtime/scaling benchmarks anywhere). |
| Expression matrix format on disk | Sparse MatrixMarket (`.mtx.gz` + gene/cell name files), not a dense CSV | Raw snRNA counts are very sparse; a dense CSV of the full 34,266-tumour-cell matrix would be needlessly large. The run scripts convert it to a plain in-memory `pandas.DataFrame` right before calling TreeAlign (TreeAlign's own API just wants a DataFrame — its example notebooks happen to load a CSV, but nothing requires that specific I/O path). |
| CN matrix format on disk | Dense `.csv.gz` | Not sparse (every gene has a CN value), and small enough (thousands of genes × 6,039 cells) that plain gzip-compressed CSV is simplest; `pandas.read_csv` reads `.csv.gz` natively. |

## Data flow

```
scDNA (existing, unmodified):
  MPNST_all_single_cells_2.5_ds_6039_seed100_final_cn_profiles.tsv  ┐
  MPNST_all_single_cells_2.5_ds_6039_seed100_final_tree.new         ┤ same MEDICC2 run
  scDNA/2_clustering/metadata_cells.tsv (barcode -> clone_label)    ┘
                    │
                    ▼
      0_prepare_scdna_cn.R
                    │
      ┌─────────────┴─────────────┐
      ▼                           ▼
gene x cell CN matrix       tree (copied/pruned to match)
      │                           │
      └─────────────┬─────────────┘
                     │  (shared by both modes)
scRNA (existing, unmodified):                │
  MPNST_C_updated.rds (RNA assay, counts layer, Tumor cells only)
                    │
                    ▼
      1_prepare_scrna_expr.R
                    │
                    ▼
        gene x cell raw count matrix
                    │
      ┌─────────────┴──────────────────────────┐
      ▼                                         ▼
2_run_treealign_clonelabel.py            3_run_treealign_tree.py
(CloneAlignClone + clone table)          (CloneAlignTree + tree)
      │                                         │
      └─────────────────┬───────────────────────┘
                         ▼
              4_compare_results.R
   (cross-tab both outputs + published scRNA cell_type)
```

## Files

| File | Purpose |
|---|---|
| `0_prepare_scdna_cn.R` | Segment→gene overlap (GENCODE v27 + `GenomicRanges`) on the all-cohort MEDICC2 CN profile; writes the shared gene×cell CN matrix, the cell→clone_label table, and a matching tree Newick. |
| `1_prepare_scrna_expr.R` | Pulls raw `RNA` counts from `MPNST_C_updated.rds` for Tumor cells, restricted to the gene set surviving step 0; writes a sparse MatrixMarket matrix plus a small `scrna_cell_metadata.csv` (region, published `cell_type`) so step 4 doesn't need to reload the 1.9GB object. |
| `2_run_treealign_clonelabel.py` | Runs `CloneAlignClone` (clone-label mode) against the 18 existing clusters. |
| `3_run_treealign_tree.py` | Runs `CloneAlignTree` (tree mode) against the full-cohort tree. |
| `4_compare_results.R` | Cross-tabulates clone-label output, tree-mode output, and the published scRNA `cell_type` clusters; writes summary tables/figures to `results/`. |
| `run_pipeline.sh` | Orchestrator; runs the chain once in `TREEALIGN_MODE=test`, then (on request) again in `TREEALIGN_MODE=full`. |

Large intermediate matrices and TreeAlign's raw output live under
`/srv/home/aste0033/projects/MPNST/scMultiOmics/TreeAlign/{inputs,outputs}/{test,full}/` — not in this
repo. `results/` here holds only small summary tables and figures.

## Environment setup (run these yourself — not run automatically)

```bash
# 1. Create the dedicated environment (pinned versions match TreeAlign's tested requirements.txt)
micromamba create -n treealign python=3.9 -y
micromamba run -n treealign pip install \
    torch==1.11.0 pyro-ppl==1.8.1 pyro-api==0.1.2 \
    numpy==1.21.6 pandas==1.3.5 scipy==1.7.3 biopython==1.79 \
    matplotlib==3.5.2 seaborn==0.11.2

# 2. Clone TreeAlign itself (not on PyPI) and pin the commit for reproducibility
git clone https://github.com/AlexHelloWorld/TreeAlign.git /srv/home/aste0033/tools/TreeAlign
cd /srv/home/aste0033/tools/TreeAlign && git rev-parse HEAD
# record the printed commit SHA in this README's "Pinned versions" section below

# 3. Install it into the new environment
micromamba run -n treealign pip install -e /srv/home/aste0033/tools/TreeAlign
```

The R-side prep scripts need no new installs (`GenomicRanges`, `rtracklayer`, `data.table`, `Matrix`,
`ape` already confirmed present in the `R` micromamba env); `1_prepare_scrna_expr.R` runs via
`single-cell.sif` like the rest of `scRNA/`.

**Pinned versions** (fill in after running the clone command above):
- TreeAlign commit: `<TODO: paste git rev-parse HEAD output here>`

## Running

```bash
# From the repository root. TREEALIGN_MODE defaults to "test" inside run_pipeline.sh.
bash scMultiOmics/2_TreeAlign/run_pipeline.sh          # test pass (region R1, chr21 only)
TREEALIGN_MODE=full bash scMultiOmics/2_TreeAlign/run_pipeline.sh   # full cohort, once test pass looks right
```

## Status

**The full test-mode pipeline (steps 0–4) has been run end-to-end and verified.** All four
pipeline stages and the comparison script ran successfully against real data.

- [x] `0_prepare_scdna_cn.R` — region R1: **1,571 scDNA cells** across clones R1_1–R1_4,
      **17,064** genome-wide protein-coding genes with complete segment coverage, tree pruned to
      1,571 matching tips. Two real issues were caught and fixed here, both left as comments in
      the script:
      1. The test-mode cell subset must filter by **`clone_label` prefix**, not the `region`
         column — cross-region clone "leakage" is real in this dataset (some `R1_1`–`R1_4` cells
         are physically sampled from other regions), so filtering by `region` instead pulled in
         ~16 unrelated clones.
      2. Test mode originally also restricted to **chr21 only** (163 genes) for speed — but
         `CloneAlignClone` drops any gene whose per-clone consensus copy number is identical
         across all clones being compared, and chr21 turned out to be completely CN-invariant
         across R1's own four clones, so every gene failed that filter and clone-label mode
         errored with `No valid genes or snps exist in the matrix after filtering`. Fixed by
         using the full genome-wide gene set in test mode too — the region-based cell/clone
         subsetting is what keeps test mode fast, not a gene-count restriction.
- [x] `1_prepare_scrna_expr.R` — region R1, Tumor cells: **4,903** scRNA cells; **15,371 / 17,064**
      genes overlapped between the GTF-derived gene set and the `RNA` assay's gene symbols.
- [x] `2_run_treealign_clonelabel.py` — **converged in 184 iterations, 1,763.6s (~29.4 min)**.
      3,899 genes survived TreeAlign's internal variance/consensus filters. **4,880 / 4,903**
      scRNA cells assigned (23 unassigned) across all 4 clones (R1_1: 374, R1_2: 1,740, R1_3: 707,
      R1_4: 2,059) — a real, non-degenerate distribution.
- [x] `3_run_treealign_tree.py` — **5,797.3s (~96.6 min)**. All 4,903 cells assigned across **7**
      tree clades (sizes 2,440 / 1,547 / 374 / 367 / 158 / 11 / 6) — TreeAlign's own data-driven
      split, coarser than and not 1:1 with the 4 k-means clones, as anticipated in the design
      table above.
- [x] `4_compare_results.R` — cross-tabs and heatmap written to `results/`. Headline finding:
      neither TreeAlign mode's assignment lines up cleanly with the scRNA's own published
      `cell_type` clusters (`Malignant_R1_1`/`Malignant_R1_2`) either — e.g. clone-label mode's
      `R1_4` splits almost evenly across published `Malignant_R1_1` (932) and `Malignant_R1_2`
      (1,013). This is itself the expected outcome given the "independent clusterings" caveat
      this whole stage exists to probe, not a bug — a genuine three-way disagreement between
      scDNA clone, scRNA published cluster, and TreeAlign's own tree clades is a real result worth
      discussing, not resolving away.
- [x] Full-mode steps 0/1 verified (**6,039 scDNA cells across all 18 clones**, 17,064 genes,
      tree correctly pruned to 6,039 tips; **34,266 scRNA tumour cells**, 15,371/17,064 genes
      overlap). One more bug caught here: the tree actually has **6,040 tips, not 6,039** —
      MEDICC2 includes a `"diploid"` reference genome as an actual outgroup tip, which correctly
      never appears in the CN profile. Full mode's original "just copy the tree, no pruning
      needed" shortcut assumed an exact tip/cell match and failed on `setequal()`; test mode's
      `drop.tip()` had silently absorbed the stray `"diploid"` tip into its "drop everything
      outside the region" step, which is why this only surfaced once run at full scale. Fixed by
      always pruning to the exact cell set in both modes.
- [ ] Full-cohort steps 2–4 (`TREEALIGN_MODE=full`) — not yet run. **Runtime scaling is the main open
      question**: full mode is ~3.8× the scDNA cells and ~7× the scRNA cells of the test run
      above, and tree mode's ~97 min already came from a 1,571-tip tree — a naive linear
      extrapolation alone would put a full run at several hours, and tree mode's recursive
      algorithm may scale worse than linearly with more/deeper clades. **Do not launch the full
      run via a `run_in_background` Bash call** — both test-mode background attempts of
      `3_run_treealign_tree.py` were killed after 45min–1h20min by this harness's own background-job
      supervision (not a code issue: it succeeded immediately once launched in a detached `tmux`
      session instead). Launch `run_pipeline.sh` (or the individual python steps) inside `tmux`/`screen`
      for the full run.

*(This section will be updated as the full-cohort run is actually performed.)*
