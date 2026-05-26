# MEDICC2 preprocessing

Scripts for running [MEDICC2](https://bitbucket.org/schwarzlab/medicc2/src/master/) on 10x scDNA copy-number profiles and unpacking bootstrap output for downstream analysis.

**Reference:** Kaufmann et al. (2022). *MEDICC2: whole-genome doubling aware copy-number phylogenies for cancer evolution.* Genome Biology, 23, 241. <https://doi.org/10.1186/s13059-022-02794-9>

---

## Files

| Script | Grouping | Bootstrap | `minseg` |
|--------|----------|-----------|----------|
| `run_medicc2_percluster.sh` | Per cluster | No | 1.5 Mb |
| `run_medicc2_perregion.sh` | Per region | No | 2.5 Mb |
| `run_medicc2_all.sh` | All cells | No | 1.5 Mb |
| `run_medicc2_percluster_bootstrapping.sh` | Per cluster | Yes (chr-wise, n=100) | 1.5 Mb |

All scripts use `--input-type t` (transposed TSV) and `--events`. The `--filter-segment-length` flag filters short noisy segments before tree inference. Output goes to `/mnt/iribhm/homes/aste0033/projects/CNA-PhyloAnalysis/MEDICC2/10X_DLP/`.

**`unpack_pickle.py`** — MEDICC2 stores bootstrap replicate trees as a pandas DataFrame in `*_bootstrap_trees_df.pickle` (unique topologies + occurrence counts). This script expands by count and writes one Newick string per replicate to a multi-tree `.new` file readable by `ape::read.tree()`. Requires `biopython` and `pandas` (e.g. `micromamba run -n python python unpack_pickle.py`).
- Input: `per_cluster_bootstrap_n20/minseg_1.5/<sample>/*_bootstrap_trees_df.pickle`
- Output: `per_cluster_bootstrap_n20/unpacked_minseg_1.5/<sample>_bootstrap_trees.new`


---

## Chromosome-wise bootstrapping

### How it works

Standard bootstrapping resamples individual segments, but adjacent segments on the same chromosome are highly correlated (arm-level events affect many bins). Resampling segments independently breaks this structure and can overstate confidence. Chr-wise bootstrapping ([`--bootstrap-method chr-wise`](https://medicc.readthedocs.io/en/latest/)) resamples whole chromosomes as units instead, preserving within-chromosome correlation. Each replicate draws *N* chromosomes with replacement (*N* = number of chromosomes in input, typically 22 autosomes), implemented in `chr_wise_bootstrap_df()` in [`medicc/bootstrap.py`](https://bitbucket.org/schwarzlab/medicc2/src/master/medicc/bootstrap.py).

### Why bootstrap branch lengths are non-integer

The final tree and bootstrap trees use different code paths in [`medicc/core.py`](https://bitbucket.org/schwarzlab/medicc2/src/master/medicc/core.py):

- **Final tree** (`ancestral_reconstruction=True`): NJ is run on integer pairwise MED distances → ancestral copy-number states are reconstructed → `update_branch_lengths()` overwrites NJ values with FST MED scores between reconstructed ancestor profiles (`float(fstlib.score(...))`). These are integer discrete event counts.
- **Bootstrap replicates** (`ancestral_reconstruction=False`): NJ is run on the resampled distance matrix and the pipeline stops there. Branch lengths are raw NJ outputs from [`medicc/nj.py`](https://bitbucket.org/schwarzlab/medicc2/src/master/medicc/nj.py): `d_left = 0.5 * D[left, right] + 1/(2*(r-2)) * (sum_row_left - sum_row_right)` — division by 2 and `2*(r-2)` produces fractional values from integer inputs.

Each chromosome is sampled with equal probability regardless of size, but larger chromosomes have more bins and contribute more events to the pairwise distances `D`. A replicate where chr1 is drawn twice has a systematically elevated `D` and correspondingly higher NJ branch lengths — this size effect is not corrected by NJ's row-sum term, which accounts for unequal rates between lineages, not for genome size differences between replicates. The per-replicate distance matrices are computed fresh per resample and **not saved to disk**; `_pairwise_distances.tsv` reflects the full original genome only.

### Using bootstrap edge distances to estimate branch length variance

Bootstrap edge distances can be used for this, but two obstacles must be understood:

1. **Different quantity.** NJ branch lengths and `update_branch_lengths()` branch lengths are not the same thing. Their means need not agree and their variances are not on the same scale, so a bootstrap NJ branch length cannot be directly compared to the corresponding integer in the final tree as a measure of "uncertainty".

2. **Variance conflates two effects.** High variance on a bootstrap edge reflects both (a) statistical uncertainty and (b) chromosome dependency — a branch driven by a single chromosomal event collapses whenever that chromosome is absent from the resample. These two effects are completely confounded.

**What bootstrap edge distances can be used for:**
- **Topology support** (primary use, reported by MEDICC2): frequency of clade recovery across replicates — unaffected by the branch-length issues above.
- **Relative stability**: edges with low CV are more robustly estimated than high-CV edges.
- **Chromosome dependency screening**: very high variance suggests the branch depends on few chromosomes.

To estimate variance in the same quantity as the final tree, `update_branch_lengths()` would need to be run on each replicate after ancestral reconstruction, with values normalised by the total segment count of each bootstrap genome. This is not currently implemented and would be computationally expensive.

### Known limitation

Because chr-wise resampling treats each chromosome as an equal-probability unit, bootstrap support for a branch driven predominantly by large chromosomes (e.g. chr1, chr3) may be inflated relative to one driven by small chromosomes (e.g. chr21, chr22).
