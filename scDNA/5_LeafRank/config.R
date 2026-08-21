# Purpose: Shared configuration for the LeafRank pipeline (scDNA/5_LeafRank/).
#          Sourced by every stage script (prepare_tree.R; rescale_wgd_tree.R
#          and diagnostics.R, both pending) so paths and parameters aren't
#          duplicated across files. Each variable can be pre-set before
#          sourcing (e.g. `group_name <- "R2"; source("scDNA/5_LeafRank/config.R")`)
#          or left to its default here, following this repo's existing
#          scripting convention (see e.g. CNA-PhyloAnalysis/code/2_phylogenies/base_trees.R).
# Input:   none (defines paths/parameters only).
# Output:  none (variables only; each stage script does its own I/O).
# Requires: invoke/source from the repository root (repo-relative `results_dir`).

### Which tree to run on ----------------------------------------------------

# "per_region" | "per_cluster" -- per_region first (per user decision), per_cluster later.
if (!exists("group_mode")) group_mode <- "per_region"

# "1.5" | "2.5" -- MEDICC2 minimum-segment-size filter subfolder.
if (!exists("minseg")) minseg <- "1.5"

# Which group within group_mode to run.
#   per_region groups:  "P", "R1", "R2", "R3", "R4", "R5"
#   per_cluster groups (for later): e.g. "scDNA_22_kmeans_10X_DLP_1.5_R1_1", "..._P_2"
#
# CAVEAT: despite the name, per_region groups are NOT strict anatomical
# partitions. Inspecting R1_final_cn_profiles.tsv directly shows cells
# prefixed P_, R1_, R2_, and R5_ all present in the "R1" tree (1571 unique
# cells). Treat group_name as an opaque tree identifier, not a region filter,
# until this is clarified.
if (!exists("group_name")) group_name <- "R1"

### Input paths (this project's current MEDICC2 output) ---------------------
# See this repo's CLAUDE.md for why an older path under
# /mnt/iribhm/homes/.../CNA-PhyloAnalysis/MEDICC2/... is stale/broken on this
# host; this one is the one actually populated and current (verified 2026-08-03).
medicc2_output_base <- "/srv/home/aste0033/projects/MPNST/scDNA/MEDICC2/default"
group_dir <- file.path(medicc2_output_base, group_mode, paste0("minseg_", minseg), group_name)

tree_file        <- file.path(group_dir, paste0(group_name, "_final_tree.new"))
cn_profiles_file <- file.path(group_dir, paste0(group_name, "_final_cn_profiles.tsv"))
cn_events_file   <- file.path(group_dir, paste0(group_name, "_copynumber_events_df.tsv"))
summary_file     <- file.path(group_dir, paste0(group_name, "_summary.tsv"))

### Output paths -------------------------------------------------------------
# Intermediate/processed pipeline outputs (processed trees, tips_WGD vectors,
# MATLAB rescaling files, LeafRank's own raw per-cell posterior RDS) go under
# the project folder, mirroring this repo's MEDICC2/SPRINTER convention
# (CLAUDE.md: "Code lives here; data never does, except small tables/results").
leafrank_output_base <- "/srv/home/aste0033/projects/MPNST/scDNA/LeafRank"
group_output_dir <- file.path(leafrank_output_base, group_mode, paste0("minseg_", minseg), group_name)

# Small, git-appropriate results (summary tables, diagnostic figures) go here,
# matching every other stage's own results/ subfolder (e.g. scDNA/4_sprinter/results/).
results_dir <- "scDNA/5_LeafRank/results"

### Sampling probability rho -------------------------------------------------
# rho = (# cells in this tree) / (estimated total tumor cell population).
# n_sampled_cells here is only a placeholder (used for early sanity-check
# printouts before any tree is loaded) -- run_leafrank.R OVERRIDES it with the
# real tip count of the loaded tree right before computing rho for real, since
# that's what actually matters to get_full_pars()/LeafRank() and what the
# calibration below was validated against.
#
# total_population = 1e10 was the user's original rough estimate ("~1e3 cells
# from 10e9"), but a full LeafRank() run on R1 with rho=1e-7 (from 1e3/1e10)
# showed severe "boundary saturation" (74.4% of tips pinned near the maximum
# fitness value -- see todo.md step 6). Progression of real full-tree runs
# since then (each also paired with the V change below):
#   V=16, pop=1e10 -> 74.4% saturated
#   V=16, pop=1e9  -> 70.3% saturated
#   V=24, pop=1e8  -> 37.6% saturated (real improvement, not fully resolved)
# Per user direction (2026-08-20), pushed to 1e7 (real full run: 32.9%
# saturated, down from 37.6% at 1e8 -- smaller but real continued
# improvement). Now pushing one more order of magnitude down to **1e6**
# (2026-08-21) -- the earlier 50-tip subsample directional check (2026-08-19)
# already flagged pop=1e6 (type 2 of 24) as nudging toward the LOW boundary,
# so this is expected to be close to the practical floor for this lever
# before risking the opposite (low-end) saturation problem -- not yet
# empirically confirmed at full-tree scale.
if (!exists("n_sampled_cells"))  n_sampled_cells  <- 1e3
if (!exists("total_population")) total_population <- 1e6
rho <- n_sampled_cells / total_population

### Multi-type branching-process parameter seed (see todo.md steps 4 and 6) ---
# V=16 (the paper's own truncal-WGD HGSOC preset) caused severe boundary
# saturation on the real R1 run: 70.3% of tips pinned within 1e-4 of the max
# birth rate (0.835). Widened to V=24 (2026-08-14/15, per user direction) to
# give the model room above the old ceiling rather than changing rho/tau
# further -- confirmed on a 50-tip subsample: 0% saturation, and only ~40%
# slower than V=16 (not the ~2.25x a naive quadratic-in-V scaling would
# predict). Real full runs at V=24 cut saturation from 70.3% -> 37.6%
# (pop=1e8) -> 32.9% (pop=1e7) -- real, if diminishing, progress. Per user
# direction (2026-08-21), widened further to **V=30** (max b_rate 3.173, see
# the V-scaling table in todo.md step 6) combined with total_population=1e6
# as a second lever. NOTE: unlike V=24, V=30 has only been checked for
# get_full_pars() numerical stability (confirmed fine) and its subsample
# TIMING (aborted mid-run once V=24 was judged sufficient at the time) --
# its actual effect on saturation/relative-type-position has NOT been
# subsample-tested the way V=24 and growth_factor=1.14 were. This run is the
# first real test of V=30. V=40 was tested and found to break get_full_pars()
# outright (root-finder failure) -- do not go that high.
if (!exists("V")) V <- 30
#
# growth_factor: the geometric step between consecutive b_rates. Default 1.1
# is the paper's own preset. Increasing it widens the birth-rate SPREAD at
# FIXED V, raising the ceiling without adding types -- unlike raising V, this
# doesn't add per-node computation cost (cost scales with V, not with the
# b_rate values themselves), so it's a cheaper lever than V=30+ for the same
# goal. Tried widening to 1.14 (2026-08-20) but REVERTED back to the default
# 1.10: a 50-tip subsample check showed factor=1.14 pushed cells to a HIGHER
# relative position within the type range (types 7-12 of 24, i.e. topping out
# ~50% of the way up) than factor=1.10 did (types 3-8, ~33% of the way up) --
# despite the higher absolute ceiling. Reason: widening growth_factor also
# raises `tau` substantially (far more than total_population does -- factor
# 1.10->1.14 moved tau from 7.09e-3 to 1.56e-2, a 2.2x jump), which shrinks
# the implied time budget T=1/tau and requires relatively higher growth rates
# to explain the same branch lengths -- directly counteracting the extra
# headroom from the raised ceiling. factor=1.14's subsample position (33-50%)
# sits higher than the best validated full-tree config's implied position
# (V=24/pop=1e8/factor=1.10, ~17-42%, which gave the real 37.6% improvement)
# -- a red flag that widening the spread this way could be a net step
# backward, not forward. Per user direction (2026-08-20), kept factor=1.10
# and used total_population (below) as the next lever instead. Checked
# numerically stable up to factor=1.20 (max b_rate 13.2) via get_full_pars()
# if this is revisited later -- no failures observed in that range, unlike
# V=40's outright break -- but any future attempt should re-check relative
# type-position on the subsample, not just saturation-at-boundary, given
# what was learned here.
if (!exists("growth_factor")) growth_factor <- 1.10
b_rates <- 0.2 * growth_factor ^ (0:(V - 1))
d_rates <- rep(0.18, V)
nu      <- 1e-4
d_t     <- 0.01

### Compute settings ----------------------------------------------------------
# Conservative default: this HPC has no queuing system (no SLURM, no
# login/compute-node separation -- see global CLAUDE.md), so every session
# shares physical cores directly with other concurrent users. Raise only with
# awareness of what else is running.
if (!exists("n_threads")) n_threads <- 8
