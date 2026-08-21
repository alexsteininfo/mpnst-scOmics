# LeafRank — execution todo

Goal: run LeafRank (Wu et al., 2026 — see `../../literature/README.md#leafrank`) on this
project's MEDICC2 copy-number trees to infer per-cell relative fitness rankings. Local tool
clone: `/srv/home/aste0033/GitHub/LeafRank`. This connects to `CNA-PhyloAnalysis/next_steps.md`
items 4, 7, 8 (WGD-aware molecular-clock correction, selection-coefficient estimator
prototype, pruning-rationale check) — once this stage produces output, that repo's planned
work can consume it.

**WGD handling decided (2026-08-03, per user direction):** this cancer is treated as having
**fully truncal WGD throughout** — the user is confident of this biologically, and considers any
MEDICC2 "subclonal WGD" call (e.g. the per-cluster 0-2-flagged-node pattern found during initial
research, or the internal-node-specific calls in per-region trees) likely an inference artifact
rather than real biology. Consequently **the full WGD-aware two-rate rescaling (MATLAB
`modified_chronos`) is skipped for now** — step 3 uses plain constant-rate `ape::chronos(tree,
model = "clock")` instead (no MATLAB needed). This is consistent with what was actually observed
on the `R1` per-region tree: the single WGD event (`internal_591`) sits so close to the tree's
root that all 1571/1571 tumor tips came back WGD+ (see step 2) — i.e. genuinely truncal within
that tree too. **Keep the WGD-aware path for later**: if a future tumor/patient genuinely needs
it, implement it as a separate script (don't reuse `rescale_tree.R`) following the two-rate path
still documented in step 3 below; MATLAB itself was tested and works (module loads) modulo the
license-server issue also noted below.

**File/output organization decided (2026-08-03):** code is split by pipeline stage — a shared
`config.R` (paths/parameters, sourced by every stage script) plus one script per stage
(`prepare_tree.R`, `rescale_tree.R`, `run_leafrank.R` — a `diagnostics.R` to follow later),
matching this repo's existing per-stage-script convention, e.g. `4_sprinter/persample/`.
Processed/intermediate outputs (tumor-clade trees, `tips_WGD` vectors, the ultrametric tree,
LeafRank's raw per-cell posteriors) go to the project folder at
`/srv/home/aste0033/projects/MPNST/scDNA/LeafRank/<group_mode>/minseg_<minseg>/<group_name>/`
(mirroring the MEDICC2/SPRINTER convention) — **not** committed to git. Only small,
git-appropriate outputs (summary tables, diagnostic figures) go to `scDNA/5_LeafRank/results/`,
matching every other stage's own `results/` subfolder. There is no shared top-level
`scDNA/results/` (none exists elsewhere in this repo; each stage keeps its own).

**Status (2026-08-03): steps 0-4 are done; step 5 (the actual LeafRank run) is implemented but
NOT executed by Claude** — an empirical timing probe (100-tip and 250-tip subsamples of `R1`,
same config, `n_threads=8`) showed even the 100-tip case hadn't finished after >10 minutes, so
the full 1571-tip run is expected to take hours, well past the user's 60-minute
run-it-yourself-in-tmux threshold. See "Running the real job" below for the exact command.

## Running the real job (step 5) — run this yourself in tmux

`rescale_tree.R` (step 3) has already been run for `R1` (~4.4 min, done). `run_leafrank.R`
(steps 4+5) is ready but the actual `LeafRank()` call is expected to take hours for the
1571-tip `R1` tree — run it in a dedicated tmux session rather than waiting on it directly:

```bash
cd /srv/home/aste0033/GitHub/mpnst-scOmics
tmux new -s leafrank_R1
# inside the tmux session:
micromamba run -n R Rscript scDNA/5_LeafRank/run_leafrank.R
# detach with Ctrl-b then d; reattach later with: tmux attach -t leafrank_R1
```

`config.R` defaults to `n_threads = 8` (conservative, since this cluster has no queuing system
and shares cores directly with other users — see global CLAUDE.md). Given the empirical timing
probe suggests the default may be quite slow, consider raising it for this one-off run if the
cluster isn't busy — e.g. run
`micromamba run -n R Rscript -e 'n_threads <- 20; source("scDNA/5_LeafRank/run_leafrank.R")'`
instead, inside the same tmux session (the paper's own benchmark used `n_threads=20`). Output
goes to `/srv/home/aste0033/projects/MPNST/scDNA/LeafRank/per_region/minseg_1.5/R1/R1_pFitness.rds`.

## 0. Open decisions — RESOLVED

- [x] R package installation into the `R` micromamba env: **approved and done** (step 1).
- [x] Starting group: **`group_mode = "per_region"` first** (per_cluster later), starting group
      `"R1"`. **Caveat found during path verification:** despite the name, `per_region` groups
      are NOT strict anatomical partitions — the `"R1"` per-region tree's
      `R1_final_cn_profiles.tsv` contains cells prefixed `P_`, `R1_`, `R2_`, and `R5_` (1571
      unique cells total). Treat `group_name` as an opaque tree identifier for now; see
      `config.R`'s comment.
- [x] Sampling probability `rho`: computed directly as `n_sampled_cells / total_population` =
      `1e3 / 1e10` = **1.0e-7** (per the user-supplied estimate "~1e3 cells from 10e9"). Note:
      this differs by ~100x from the "~10e-6" figure mentioned verbally alongside it — left as
      a derived quantity in `run_leafrank.R` (not hardcoded) so either input can be corrected
      directly if 1e10 isn't the intended total population.
- [x] MATLAB: tested directly (2026-08-03). `module load MATLAB/2022b-r9` succeeds and the
      `matlab` binary runs, but license checkout fails: `License checkout failed. License
      Manager Error -15 -- Unable to connect to the license server.` This cluster has no
      login/compute-node distinction (no queuing system — see global CLAUDE.md), so this is a
      genuine license-server connectivity issue to resolve, not something explained by "try a
      different node" — before step 3 can rely on `fmincon`, this needs troubleshooting (e.g.
      license file/path/network config) or the documented `nloptr` fallback.

## 1. Environment setup — DONE

- [x] Installed into the `R` micromamba env via
      `micromamba install -n R -c conda-forge -c bioconda r-brobdingnag r-expm bioconductor-ggtreeextra r-r.matlab -y`
      — all 4 were available on conda-forge/bioconda matching the existing R 4.5.3 build; no
      CRAN fallback needed. Verified all of `Brobdingnag`, `expm`, `ggtreeExtra`, `R.matlab`,
      plus the previously-present `ape`, `doParallel`, `foreach`, `tibble`, `ggtree`, `ggplot2`,
      `devtools` load via `requireNamespace()`.
- [x] Installed LeafRank itself via
      `micromamba run -n R Rscript -e 'devtools::install_local("/srv/home/aste0033/GitHub/LeafRank", upgrade = "never")'`
      — succeeded (slow: ~5 min, mostly spent tar-bundling the repo's large `Results_analysis/`
      example-results directory, not actual compilation). Verified `library(LeafRank)` loads and
      `LeafRank::LeafRank` / `get_full_pars` / `get_ultrametric_prepared` / `get_ultrametric` are
      all exported with the signatures found during research (confirms `n_threads`, not
      `num_threads`, as the real parameter name).
- [x] `module load MATLAB/2022b-r9` tested — see "MATLAB" bullet above; module/binary OK,
      license checkout currently fails. Note: this cluster has no login/compute-node
      distinction (no queuing system — recorded in the global CLAUDE.md), so this is a real
      connectivity issue to resolve, not a "wrong node type" problem.

## 2. Prepare input tree — DONE (implemented in `prepare_tree.R`)

- [x] Load `<group_name>_final_tree.new` via `ape::read.tree()`. **Confirmed current path**
      (verified 2026-08-03, supersedes the older Haixi-directory paths found during initial
      research — see `config.R`):
      `/srv/home/aste0033/projects/MPNST/scDNA/MEDICC2/default/<group_mode>/minseg_<minseg>/<group_name>/`
      e.g. `.../per_region/minseg_1.5/R1/R1_final_tree.new`, or for per_cluster later,
      `.../per_cluster/minseg_1.5/scDNA_22_kmeans_10X_DLP_1.5_R1_1/..._final_tree.new`.
- [x] Check `ape::is.binary(tree)`; if not strictly bifurcating, run `ape::multi2di(tree)` first —
      LeafRank's message-passing code (`up_m_des`, `up_m_sib`, `get_sibling`) hard-assumes
      exactly two children per internal node and will silently use only the first two on a
      polytomy (no error, just wrong results). **Verified on `R1`: already strictly bifurcating**
      (MEDICC2's own output is binary here), so `multi2di()` was a no-op — but the check runs
      every time in case a future group isn't.
- [x] Extract the tumor clade only (drop the diploid pseudo-tip / root trunk) — implemented as
      `extract_tumor_clade()` in `prepare_tree.R`, ported from the exact logic in
      `tree_to_plot_data()` (`CNA-PhyloAnalysis/code/2_phylogenies/base_trees.R:44-61`): find the
      root node (`setdiff(tree$edge[,1], tree$edge[,2])`), its non-diploid child, and
      `ape::extract.clade()` on that child (falling back to `ape::drop.tip(tree, "diploid")` if
      the non-diploid child is itself a single tip). **Verified on `R1`: 1572 raw tips (incl.
      diploid) → 1571 tumor tips.**
- [x] Derive the per-tip WGD status vector (`tips_WGD`): implemented in `prepare_tree.R` —
      `<group_name>_copynumber_events_df.tsv` rows with `type == "wgd"` name the internal node(s)
      where a WGD event occurred; all descendant tips of each such node (via
      `ape::extract.clade()`) are marked WGD+ (`1`), positionally aligned to
      `tree_tumor$tip.label` order (LeafRank's `get_ultrametric_prepared` indexes by position,
      not by name). Reuses the same logic pattern as `has_wgd` in
      `CNA-PhyloAnalysis/code/5_frequencies/near_clonal_cnas.R` (~line 327). **Verified on `R1`:
      found exactly 1 WGD event on node `internal_591`, matching `R1_summary.tsv`'s
      `wgd_status: WGD on branch internal_591` — but this node sits so close to the tumor-clade
      root that ALL 1571/1571 tips came back WGD+ (see the top-of-file note on what this implies
      for step 3).**
- [x] Note: `<group_name>_final_cn_profiles.tsv`'s `is_wgd` column is a **per-segment** flag
      (True only for the ~247Mb region directly affected by the WGD event itself, on the
      `internal_591`-type row), not a per-sample/per-tip WGD status column — `prepare_tree.R`
      correctly uses the event-based descendant propagation above instead, not this column.
- [x] Outputs saved to the project folder (not git), per the file-organization decision above:
      `/srv/home/aste0033/projects/MPNST/scDNA/LeafRank/per_region/minseg_1.5/R1/R1_tumor_tree.rds`
      and `..._tips_wgd.rds`.

## 3. Rescale to an ultrametric tree — DONE (implemented in `rescale_tree.R`, using the simplified path)

- [x] **Decision (2026-08-03, per user direction): treat this cancer as fully truncal-WGD
      throughout** — skip the two-rate WGD-aware model entirely (not just for `R1`, for every
      group). Implemented as plain constant-rate `ape::chronos(tree_tumor, model = "clock")`
      (ape's built-in shortcut for a strict molecular clock — no MATLAB needed). **Verified on
      `R1`**: 1571 tips, 263 zero-length branches handled fine, fit completed in ~261 sec (~4.4
      min), converged (with a benign "false convergence (8)" warning from the underlying
      optimizer — common on large `chronos` fits, checked not to matter: root-to-tip distance
      is exactly 1.0 for every tip as expected). Saved to
      `.../LeafRank/per_region/minseg_1.5/R1/R1_ultrametric_tree.rds`.
- [ ] **DEFERRED, not deleted** — the two-rate WGD-aware path (`library(R.matlab)` +
      `get_ultrametric_prepared(phy, tips_WGD, MAT_path=...)` + `module load MATLAB/2022b-r9` +
      `modified_chronos(...)` from `LeafRank/MATLAB/modified_chronos/modified_chronos.m` +
      `get_ultrametric(...)`) is kept here as a reference for if a future tumor/patient genuinely
      needs it — do not reuse `rescale_tree.R` as-is for that; write a separate script. Notes if
      that day comes: reconcile `modified_chronos.m`'s actual default output path (differs from
      the README's stated `MATLAB/ultra-tree.csv`); if the Optimization Toolbox license issue
      (step 0/1) is still unresolved, `nloptr::nloptr` (`NLOPT_LD_SLSQP`) is a documented
      fallback for the underlying `fmincon` call (smooth Poisson NLL, analytic gradient
      available in `poisson_WGD.m`/`poisson_WGD_grad.m`); sanity-check any such rescaling with
      the paper's own hanging-rootogram + runs-test diagnostic before trusting it.

## 4. Configure LeafRank parameters — DONE (implemented in `run_leafrank.R`)

- [x] Sampling probability `rho` = 1.0e-7, from `config.R` (see step 0).
- [x] Branching-process seed from `config.R`: **V=16** fitness types, `b_rates = 0.2*1.1^(0:15)`
      (geometric ×1.1, the paper's truncal-WGD HGSOC preset), `d_rates = rep(0.18, 16)`,
      `nu = 1e-4`, `d_t = 0.01`.
- [x] `get_full_pars(b_rates=b_rates, d_rates=d_rates, nu=nu, rho=rho, tree=phy, model="default")`
      run against the real `R1` ultrametric tree — fast (0.58 sec), returns an **unnamed** list
      `[b_rates, d_rates, nu, rho, tau]` (index by position, not name). Fitted **tau ≈ 4.48e-3**
      for `R1` (implying internal model time `T = 1/tau ≈ 223.2`).
- [x] **Important correction found while implementing this**: the README's own
      `timeTo <- ceiling(max(get_depth(phy))/10)*10` formula is only correct when
      `time_scale (tau) == 1` (true for the README's raw-simulation-time demo tree, false here).
      `node_time_to_present()` (`LeafRank/R/tree_utils.R:112-131`) confirms `time_scale` is a
      **divisor** applied to raw tree-branch depths to get internal model time — so for our
      tau-fitted, unit-normalized (root-to-tip=1) ultrametric tree, `T_vector` must use the
      generalized form `ceiling(max(get_depth(phy)) / time_scale / 10) * 10` instead, otherwise
      the model-time window would be ~223x too short (a "boundary saturation"-style failure,
      silently wrong rather than an error). Implemented this way in `run_leafrank.R`.
- [x] `nu` expanded to a length-V vector (`rep(nu, V)`) before the `LeafRank()` call itself, per
      the README's own demo pattern (distinct from `get_full_pars()`'s scalar `nu` input).
- [ ] Not used: the `model` argument's `'aggressive 2'`/`'agressive2'` preset strings (noted
      during earlier research as a real README bug) — we pass explicit `b_rates`/`d_rates`
      instead of using a named preset, so this doesn't apply to the current config, but keep the
      correct spelling (`'aggressive 2'`, with a space) in mind if a preset is used later.

## 5. Run LeafRank — implemented, NOT executed by Claude (run it yourself, see above)

- [x] `run_leafrank.R` calls `LeafRank(phy, outFile, rho, d_t, time_scale, b_rates, d_rates,
      nu_vec, T_vector, non_negativity_cutoff, n_threads)` — confirmed the parameter is
      **`n_threads`**, not `num_threads` as the README's own usage example incorrectly implies.
- [x] **Empirical timing check performed** (2026-08-03): subsampled `R1` to 100 and 250 tips
      (matching the paper's own benchmark scale) and ran `LeafRank()` with the same config
      (`n_threads=8`). Even the 100-tip case had not finished after >10 minutes (only printed
      LeafRank's internal `"start: E_list"` progress line) — stopped the probe rather than wait
      further, since this was already conclusive: the full 1571-tip `R1` tree is expected to
      take **hours**, well past the 60-minute run-it-yourself threshold. **Did not run the full
      job** — see "Running the real job" above for the exact tmux command.
- [x] `LeafRank()` itself only returns the string `"success!"` — the actual results
      (`phylo`, `meanFitness`, `upMessages`, `downMessages`, `marginalProb`, `input_params`) are
      `saveRDS()`-written to `outFile` = `<group_output_dir>/<group_name>_pFitness.rds` (the
      project folder — this is raw/intermediate output, not a small final summary, so it doesn't
      belong in the repo's `results/`).
- [x] **Real run completed** (2026-08-04, user-run in tmux per the handoff above): `n_threads=20`,
      wall-clock **432.1 min (~7.2 h)** for the full 1571-tip `R1` tree — confirms the earlier
      timing-probe conclusion that this is a many-hours job, not a >60-min-but-tractable one.
      Output verified structurally sound: `meanFitness` length 3141 (tips+internal, no
      NA/NaN/Inf), `marginalProb` is a proper 3141×16 probability matrix (every row sums to
      exactly 1). No warnings/errors in the captured tmux log. **But see step 6 — the result
      itself shows a real "boundary saturation" problem, so this output is not yet trustworthy
      for interpretation as-is.**

## 6. Diagnose & interpret

- [x] **Boundary saturation confirmed on the `R1` run** (2026-08-04) — this is exactly the
      failure mode the paper documents (Results, Fig. 3), found on the very first real run:
      - 1168/1571 tips (74.4%) have `meanFitness` within 1e-4 of the maximum possible value
        (0.8354496, the birth rate of fitness type 16 of 16).
      - 1371/1571 tips (87.3%) have fitness **type 16** (the top type) as their single
        most-probable type; the remaining tips are spread thinly across types 4-15 (types 1-3
        got zero tips).
      - Mean posterior mass on the top type across all tips: **86.9%** (vs. 2.5e-7 on the
        bottom type) — i.e. the model is highly confident nearly everyone is near-maximal
        fitness, which flattens the ranking exactly as the paper describes ("obscuring the
        identification of high-fitness clones when most lineages are pushed to the maximum").
      - **Likely root cause**: per the paper, boundary saturation "arises... when the
        time-scaling factor is either too small or too large, or when the birth/death rates are
        incongruent with the tree topology." Our `tau` (time_scale ≈ 4.48e-3) was auto-fit by
        `get_full_pars()` to match `rho = 1e-7`, which itself came from the rough
        `n_sampled_cells=1e3 / total_population=1e10` estimate flagged as uncertain back in step
        0 (recall the ~100x discrepancy vs. the "~10e-6" figure mentioned verbally). An
        extremely small `rho` forces `get_full_pars()` to solve for a `tau` implying a very
        long elapsed time window (`T = 1/tau ≈ 223`), which plausibly requires most lineages to
        have grown near the fastest possible rate to reach the sampled population size in that
        window — exactly the saturation pattern observed. **Revisit the `total_population`
        estimate (and/or `rho` directly) before the next run** rather than re-running as-is;
        this is a parameter-calibration issue, not a code bug (the paper's own Fig. 3
        prescribes exactly this diagnostic-then-recalibrate loop as normal LeafRank usage, not
        a special failure).
      - **Cost implication**: because a full `R1` run takes ~7.2 h, iterating on this
        calibration at full per-region scale is expensive — consider testing candidate
        `rho`/`total_population` values on a small subsample first (as was done for the timing
        probe) before committing another multi-hour run.
- [x] **Calibration probe run (2026-08-04)** — `scDNA/5_LeafRank/` scratch script (not committed
      to repo, throwaway): fit `tau` via `get_full_pars()` against the **real 1571-tip tree**
      (so the calibration is faithful) for 3 candidate `total_population` values, then ran the
      actual `LeafRank()` call on a fast 50-tip random subsample of `R1` for each (same
      `n_threads=8`, ~20-21 min per candidate). All three avoided saturation entirely:

      | `total_population` | `rho` | `tau` | meanFitness range | types used | % saturated (max/min) |
      |---|---|---|---|---|---|
      | 1e7 | 1.571e-4 | 4.717e-3 | 0.223-0.333 | 2-6   | 0% / 0% |
      | 1e8 | 1.571e-5 | 4.640e-3 | 0.242-0.389 | 3-8   | 0% / 0% |
      | 1e9 | 1.571e-6 | 4.566e-3 | 0.247-0.466 | 3-10  | 0% / 0% |

      All three: 35/50 distinct `meanFitness` values (good resolution), clean monotonic trend
      (higher assumed population → slightly lower `tau` → distribution shifts to higher/wider
      fitness types) — a sharp contrast with the original `total_population=1e10` run, which put
      74.4% of tips within 1e-4 of the ceiling. **The failure threshold sits somewhere between
      1e9 and 1e10** — not pinned down further yet. **Open decision for the user**: which
      `total_population` to commit to for the next full `R1` run (or a value in between/near
      1e10 to test how close to the original ~1e10 belief is still viable) — this is a
      biological judgment call (how large is this tumor/region's actual cell population), not
      something to guess further without input.
- [x] **User decision (2026-08-04): `total_population = 1e9`.** Updated in `config.R` as the new
      default. **Consistency bug found and fixed while updating this**: the calibration probe
      computed `rho` from the tree's REAL tip count (1571), but `config.R`'s `rho` was still
      being computed from the placeholder `n_sampled_cells = 1e3` — these would have given
      *different* `rho` values (1e-6 vs. the validated 1.571e-6) had this gone unnoticed.
      Fixed by having `run_leafrank.R` **override** `n_sampled_cells`/`rho` with the real tree's
      tip count immediately after loading it (`config.R`'s own `rho` is now only a
      before-any-tree-is-loaded placeholder). Verified directly: with this fix,
      `total_population=1e9` + the real 1571-tip `R1` tree reproduces `rho=1.5710e-6` and
      `tau=4.5655e-3`, exactly matching the validated calibration-probe numbers. **Ready to
      launch the real run** with the same tmux command as before (`config.R` already defaults
      to the corrected value, no override needed).
- [x] **Second full `R1` run completed (2026-08-14/15) — `total_population=1e9`, `V=16`.**
      426.1 min (~7.1 h), clean completion, no errors. **But the subsample-based calibration did
      NOT transfer to the full tree**: still 70.3% of tips within 1e-4 of max fitness (vs. 74.4%
      for the original `total_population=1e10` run — only a ~4-point improvement), 85.7% of tips
      with the top type as most-probable, 84.4% mean posterior mass on the top type. **Root
      cause identified**: the 50-tip random-subsample calibration probe is not a faithful proxy
      for the full tree's saturation behavior — dropping ~97% of tips collapses/merges
      intermediate branches, so the induced subtree has far fewer, longer branches spanning
      roughly the same root-to-tip time span. That artificially lowers the branching *density*
      per unit of model-time the fitted parameters need to explain, making lower growth rates
      look sufficient on the subsample when the real (much more densely branching) full tree
      still needs near-maximal growth to explain its topology. **Lesson: tip-subsampling changes
      branching density and is not a reliable proxy for boundary-saturation testing** — future
      calibration attempts should either use much larger (less distorting) subsamples or a
      different diagnostic entirely, not small random-tip subsamples.
- [x] **Decision (2026-08-15, per user direction): widen the birth-rate spectrum instead of
      further adjusting `rho`/`tau`.** Rationale: `rho`/`total_population` changes alone weren't
      fixing it (see above), and repeated full 7h+ runs to hunt for a working `rho` value is
      expensive; widening `V` gives the model headroom above the old ceiling directly. Checked
      first that this wouldn't blow up runtime: `get_full_pars()` (fast, <1 sec, run against the
      real 1571-tip tree) confirmed **V=40 breaks outright** (root-finder failure, "missing
      value where TRUE/FALSE needed") — there's a hard ceiling somewhere between 30 and 40.
      Empirically timed `V=16` vs. `V=24` on the 50-tip subsample: 1266 sec vs. 1775 sec (+40%
      for +50% more types — sub-linear-ish, not the ~2.25x a naive quadratic-in-V scaling would
      predict), both showing 0% saturation on the subsample (expected — the subsample-proxy
      issue above applies regardless of `V`, so this doesn't re-validate the fix, only confirms
      the config runs and estimates the cost). **Chose `V=24`** (max birth rate 1.791, net top
      growth rate 1.611 vs. the old 0.655) — a `V=30` timing run was started but aborted by the
      user once `V=24` was judged good enough (`V=30` had already been confirmed numerically
      viable via `get_full_pars()`: tau=9.7548e-3, T=102.5, so it remains a documented fallback
      if `V=24` still saturates). Updated `config.R` default to `V=24`; verified directly against
      the real 1571-tip tree that this reproduces the exact validated `tau=6.9448e-3`. **Cost
      estimate for the next full run: ~10h** (7.1h × ~1.4, extrapolating the subsample timing
      ratio) — **not yet empirically confirmed at full-tree scale, since (per the lesson above)
      the subsample can't validate whether saturation is actually resolved. This will only be
      known after the next full run completes.**

- [x] **Second lever added (2026-08-15, per user direction): lowered `total_population` one more
      order of magnitude, from `1e9` to `1e8`.** Rationale: the tested range (`1e7`-`1e9`) showed
      no sign of a lower saturation boundary (only the upper cliff between `1e9` and `1e10` is
      known), so there's room to move further from the original `1e10` guess in exchange for
      safety margin — stacked as a second, independent lever alongside `V=24` rather than relying
      on either alone. Verified against the real 1571-tip tree with `V=24`: `rho=1.571e-5`
      (10x the `1e9` value), `tau=7.0145e-3` (only ~1% different from `V=24`'s `tau` at `1e9`
      — consistent with the paper's own observation that `tau` is fairly insensitive to `rho`'s
      order of magnitude; most of the expected effect here comes from `V=24` itself, not this
      change). Updated `config.R` default to `total_population=1e8`.

**Status (2026-08-15): ready to launch the third full `R1` run, with `V=24` + `total_population=1e8`
(two independent, stacked levers) as the fix for boundary saturation. User has explicitly
accepted that this cannot be further validated on a subsample before committing the ~7-10h run
(see the subsample-proxy lesson above) — this is the first real test of both changes together.**
Same tmux pattern as before:
```bash
cd /srv/home/aste0033/GitHub/mpnst-scOmics
tmux new -s leafrank_R1_v4
micromamba run -n R Rscript -e 'n_threads <- 20; source("scDNA/5_LeafRank/run_leafrank.R")'
```
Expect ~10 h (uncertain — first empirical confirmation of `V=24` at full-tree scale). If this
run *still* saturates, next things to try, in order of how much has already been de-risked:
`V=30` (already confirmed numerically viable, not yet empirically timed/tested), or revisit
whether `d_rates` (currently flat 0.18 for all types) rather than `b_rates` is the more
appropriate lever, or accept the paper's own point that partial saturation still leaves the
non-saturated tail informative.

- [x] **Cheap V=24 directional check (2026-08-19), per user request**: before committing another
      ~11h run, re-ran the 50-tip subsample probe at `V=24` across `total_population` in
      `{1e8, 1e7, 1e6}` (same subsample/seed as all prior probes). Clean, monotonic result: lower
      population consistently shifted the subsample's fitness distribution to lower types (1e8:
      types 4-10, range 0.265-0.477; 1e7: types 3-8, range 0.244-0.398; 1e6: types 2-6, range
      0.223-0.321) with 0% saturation at either extreme for all three — confirms the *direction*
      matches what the two real full-tree V=16 runs showed (lower population → less top-end
      saturation), though as always the subsample can't confirm the *magnitude* on the full tree.
      Noted `1e6` was already nudging toward the *low* boundary (type 2, close to the floor) —
      recommended not going below `1e6`, and settled on keeping `total_population=1e8` (already
      configured) as the value to actually run at full scale.
- [x] **Third full `R1` run completed (2026-08-19/20) — `V=24`, `total_population=1e8`.**
      679.8 min (~11.3 h, a bit more than the ~10h estimate but the right order of magnitude), no
      warnings/errors, structurally clean (`marginalProb` rows all sum to 1, no NA/NaN/Inf).
      **Real, substantial improvement — the fix is working, just not completely:**

      | Run | Config | % tips at max fitness | Mean posterior mass on top type | Distinct meanFitness values (of 1571 tips) |
      |---|---|---|---|---|
      | 1 | V=16, pop=1e10 | 74.4% | 86.9% | — (538/3141 incl. internal nodes) |
      | 2 | V=16, pop=1e9  | 70.3% | 84.4% | 262 (16.7%) |
      | **3** | **V=24, pop=1e8** | **37.6%** | **58.7%** | **526 (33.5%)** |

      Saturation nearly halved vs. run 2, resolution roughly doubled. Top-type distribution now
      genuinely spans the full type range (1, 5-24 all represented — not clustered almost
      entirely on one type as in runs 1-2), though types 17-24 still hold the bulk of tips
      (~86%) and type 24 alone is the single largest bucket (929/1571 = 59.1% have it as their
      *most probable* type, though only 37.6% are within 1e-4 of the literal ceiling — i.e. many
      "leaning toward max" cells still carry real, differentiating posterior spread rather than
      being fully clipped). **Interpretation: both stacked levers (V=24, pop=1e8) meaningfully
      reduced boundary saturation but did not fully eliminate it.** Output at
      `.../LeafRank/per_region/minseg_1.5/R1/R1_pFitness.rds` (this run's output; overwrote run
      2's file, which was not preserved separately).
      **Open decision**: push further (e.g. `V=30`, or `total_population` below `1e8`,
      cross-referencing the directional check above which suggests not going below `~1e6`) for
      another ~11h+ run, or treat 37.6%/58.7% as an acceptable working resolution and move on to
      step 6's visualization/robustness-check work and step 7 documentation.
- [x] **New lever introduced (2026-08-20), per user direction: widen the birth-rate SPREAD at
      fixed V=24, instead of raising V further or changing `total_population` more.** Added a
      `growth_factor` parameter to `config.R` (`b_rates <- 0.2 * growth_factor^(0:(V-1))`,
      default was implicitly `1.1`). Rationale: unlike `V`, this doesn't add per-node
      computation (message-passing cost scales with the *number* of types, not their specific
      values), so it's a cheaper way to raise the ceiling than `V=30`. Checked numerical
      stability via `get_full_pars()` (fast, real 1571-tip tree) across `growth_factor` in
      `{1.10, ..., 1.20}` — all solved fine (unlike `V=40`'s outright failure); `tau` is highly
      sensitive to this parameter (much more than to `total_population`: factor=1.10 -> tau=
      7.09e-3, factor=1.20 -> tau=4.66e-2, a ~6.6x swing vs. `total_population`'s ~1% swing
      across three orders of magnitude). **Chose `growth_factor=1.14`** (max b_rate 4.07, a
      proportionate step similar in relative size to the `V=16`->`V=24` ceiling increase)
      rather than the more extreme values, which push the top growth rate to values (9-13) that
      seem disproportionate given the low end is still anchored at 0.2.
      **50-tip subsample check** (same seed as all prior probes, `V=24`, `total_population=1e7`):
      `factor=1.10` (baseline) -> types 3-8, range 0.244-0.398; `factor=1.14` -> types 7-12,
      range 0.437-0.866. Both 0% saturated at either extreme, similar resolution (36 vs. 35
      distinct values of 50) — clean, stable result, cells occupy the lower-middle of the wider
      range rather than clustering near the new ceiling. (One run of this probe was lost mid-way
      to a session interruption that killed the background process; relaunched successfully with
      `setsid`/`nohup`/`disown` so it survives future session restarts -- worth using that
      pattern for any future long-running probe.) As always, subsample results don't confirm
      full-tree saturation is resolved (see the tip-subsampling/branching-density lesson above)
      — this is a stability/sanity check, not a validation of the fix's magnitude.
      **`config.R` now defaults to `V=24`, `growth_factor=1.14`, `total_population=1e7`** —
      ready for the next full `R1` run (not yet launched).
- [x] **`growth_factor` widening REVERTED back to `1.10` (2026-08-20)**, per a sharp catch by
      the user: comparing the two subsample results *by relative position within the type range*
      (not just absolute saturation status) shows `factor=1.14`'s cells topped out at type 12/24
      (~50% of the way up) vs. `factor=1.10`'s type 8/24 (~33% of the way up) — i.e. despite the
      higher absolute ceiling, `factor=1.14` pushed cells proportionally *closer* to it, not
      further away. Root cause: widening `growth_factor` also raises `tau` substantially (far
      more sensitively than `total_population` does — `1.10`->`1.14` moved tau from `7.09e-3` to
      `1.56e-2`, a 2.2x jump vs. total_population's ~1% swing across 3 orders of magnitude),
      which shrinks the implied time budget `T=1/tau` and requires relatively higher growth
      rates to explain the same branch lengths — directly counteracting the raised ceiling.
      `factor=1.14`'s subsample position (33-50%) sits *higher* than the best validated
      full-tree config's implied subsample position (`V=24/pop=1e8/factor=1.10`, ~17-42%, which
      gave the real 37.6% improvement) — a red flag it could be a net step backward rather than
      forward, not confirmed either way without an actual (expensive) full run. **Decided not to
      risk an ~11h run to find out** — reverted `growth_factor` to `1.10` in `config.R`, keeping
      `total_population=1e7` (already independently supported by the earlier directional check)
      as the sole active change from the last full run. **Lesson for future spread-widening
      attempts** (via `V` or `growth_factor`): check relative type-position on the subsample,
      not just whether saturation-at-the-boundary is 0% — both levers raise `tau`/shrink `T` as
      a side effect, so "0% saturated on a 50-tip subsample" alone doesn't distinguish a real
      improvement from one that's merely not-yet-visible-at-this-scale.
      **`config.R` now defaults to `V=24`, `growth_factor=1.10` (default/unchanged),
      `total_population=1e7`** — verified against the real 1571-tip tree
      (`rho=1.571e-4, tau=7.0855e-3`, matching the earlier validated probe exactly). **Ready for
      the next full `R1` run** (not yet launched).
- [x] **Fourth full `R1` run completed (2026-08-20/21) — `V=24`, `growth_factor=1.10`,
      `total_population=1e7`.** 682.4 min (~11.4h), clean completion, no warnings/errors,
      structurally sound (`marginalProb` rows sum to 1, no NA/NaN/Inf, config confirmed exactly:
      `rho=1.571e-4, tau=7.0855e-3`). **Continued real improvement:**

      | Run | Config | % tips at max fitness | Mean posterior mass on top type | Distinct values (of 1571) |
      |---|---|---|---|---|
      | 1 | V=16, pop=1e10 | 74.4% | 86.9% | — |
      | 2 | V=16, pop=1e9  | 70.3% | 84.4% | 262 (16.7%) |
      | 3 | V=24, pop=1e8  | 37.6% | 58.7% | 526 (33.5%) |
      | **4** | **V=24, pop=1e7** | **32.9%** | **47.0%** | **590 (37.6%)** |

      Smaller step than the `V=16`->`V=24` jump (expected, `total_population` alone has a more
      modest effect per the earlier calibration data), but a real, consistent improvement:
      mean posterior mass on the top type dropped below 50% for the first time; top-type
      distribution spans nearly the full range (1, 5-24 all represented, type 24 down to
      752/1571=47.9% as most-probable choice from 929/1571=59.1% in run 3). Output at
      `.../LeafRank/per_region/minseg_1.5/R1/R1_pFitness.rds` (overwrote run 3's file).
      **Open decision, same as after run 3**: push further (e.g. `total_population=1e6`, the
      lowest value the directional check found still safely mid-range and not yet approaching
      the low boundary) for another ~11h+ run with diminishing-but-still-real expected returns,
      or treat the current 32.9%/47.0% resolution as sufficient and move to step 6
      (visualization/robustness checks) and step 7 (documentation).
- [x] **Fifth full `R1` run configured (2026-08-21), per user direction: `V=30` + `total_population=1e6`,
      both levers pushed simultaneously.** `total_population=1e6` is the value the earlier
      directional subsample check flagged as already nudging toward the LOW boundary (type 2 of
      24) -- likely close to the practical floor for this lever. **`V=30` is going into a full
      run for the first time without a saturation-specific subsample check**: it was only
      confirmed numerically stable via `get_full_pars()` (`tau=9.7548e-3` at pop=1e9 back when
      first explored) and its subsample *timing* run was aborted mid-way once `V=24` was judged
      sufficient at the time -- unlike `V=24` and `growth_factor=1.14`, its effect on
      saturation/relative-type-position has never actually been checked on the subsample. This
      run is therefore a real gamble on an unvalidated combination, accepted per explicit user
      direction. Verified config resolves against the real 1571-tip tree:
      `V=30, growth_factor=1.10, total_population=1e6 -> rho=1.571e-3, tau=9.9795e-3, T=100.2`.
      **Runtime estimate (rough, not empirical): ~13-15h** -- extrapolating the `V=16->V=24`
      subsample timing ratio (+40% for +8 types) proportionally to `V=24->V=30` (+6 types, a
      further ~+30%) on top of run 4's 682 min baseline. Not yet launched.
- [x] **Cluster courtesy check (2026-08-21)**: user asked to raise LeafRank's scheduling
      priority given "a few processes running on the cluster" -- checked and found this host is
      heavily oversubscribed (`uptime` load average 151.93-166.48 on only 56 cores, ~2.7-3x
      overcommitted), dominated by two other users' processes (~17 and ~14 cores each) plus the
      user's own concurrently-running `scMultiOmics/2_TreeAlign/2_run_treealign_clonelabel.py`
      (~12 cores). **Negative `nice` (true higher-than-default priority) is not permitted for
      this account** -- confirmed via `ulimit -e` (nice ceiling `0`) and a direct test
      (`nice -n -5` -> "Permission denied"); no admin/root access available. Per user direction,
      reversed the ask: launch LeafRank with **`nice -n 10`** instead (a real, positive
      niceness) so it yields CPU to other users under contention rather than trying to compete
      for scarce cycles, given how oversubscribed the host already is. This doesn't change
      `n_threads`/`config.R` -- niceness affects scheduling priority under contention, not how
      many cores the job can use when they're free, so it's applied at launch time only:
      `nice -n 10 micromamba run -n R Rscript -e 'n_threads <- 20; source("scDNA/5_LeafRank/run_leafrank.R")'`
      (inherited by LeafRank's internal parallel workers automatically, since child processes
      inherit niceness on Linux). Given the load situation, actual wall-clock time for this run
      may end up *longer* than the ~13-15h estimate above, not shorter.

- [ ] Plot the per-cell fitness-distribution heatmap (paper Fig. 3C pattern) to check for the
      two failure signatures the paper documents: **ranking instability** (sibling cells with
      wildly different fitness — fix by increasing `rho`, decreasing `nu`, or increasing `V`)
      and **boundary saturation** (most cells pushed to fitness extremes — fix by adjusting
      `tau` or the birth/death rates). *(Numeric confirmation of boundary saturation done above;
      the visual heatmap itself is still pending.)*
- [ ] Run the paper's own robustness check: leave-20%-out subsampling, 50 replicates, re-run
      LeafRank on each, check Spearman correlation with the full-tree ranking — this is the
      paper's recommended proxy for accuracy given no ground truth exists for a real tumor.
      **Cost warning**: given step 5's timing finding (a single full `R1` run likely takes
      hours), 50 replicates at 80% size (~1257 tips each) would be a very large undertaking —
      consider fewer replicates and/or running this on a smaller `per_cluster` tree once that
      stage is reached, rather than at full `per_region` scale.
- [ ] Visualize with `ggtree` + `ggtreeExtra` (circular tree + fitness ring), matching the
      README demo pattern.

## 7. Document

- [ ] Once a first successful run exists, write `scDNA/5_LeafRank/README.md` following this
      repo's per-stage documentation convention (purpose, inputs, outputs, parameter choices,
      any approved third-party patches) — fold this todo.md into it rather than leaving both
      as permanent separate files, per CLAUDE.md's stage-doc convention.
- [ ] Update `literature/README.md`'s LeafRank checklist row from "Not run" once this stage
      actually runs.
- [ ] Consider (separately, not part of this task) whether to update
      `CNA-PhyloAnalysis/next_steps.md` items 4/7/8 to point at this stage's output once it
      exists.
