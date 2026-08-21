#!/usr/bin/env Rscript
#
# Purpose: Finalize the multi-type branching-process parameters (step 4) and
#          run LeafRank itself (step 5) on the prepared/rescaled tree.
# Input:   <group_output_dir>/<group_name>_ultrametric_tree.rds (from rescale_tree.R)
# Output:  <group_output_dir>/<group_name>_pFitness.rds -- LeafRank's own
#          saveRDS() output: list(phylo, meanFitness, upMessages, downMessages,
#          marginalProb, input_params)
# Requires: invoke from the repository root: micromamba run -n R Rscript scDNA/5_LeafRank/run_leafrank.R
#           R packages: LeafRank (local install), ape, doParallel, foreach,
#           Brobdingnag, expm.
#
# TIMING WARNING: this is the expensive step. LeafRank's own paper cites ~30
# min for a 250-cell tree at n_threads=20; the R1 per-region tree has 1571
# tips (~6x larger). Do NOT just run this and wait -- see
# scDNA/5_LeafRank/todo.md for the empirical timing estimate on this hardware
# and the tmux command to use if it's expected to exceed ~60 min.
#
# NOTE ON time_scale / T_vector: the README's own demo computes
# `timeTo <- ceiling(max(get_depth(phy))/10)*10` directly, but that is only
# correct when time_scale (tau) == 1 (true for their raw-simulation-time demo
# tree). Our tree is normalized to root-to-tip = 1 by rescale_tree.R, and tau
# is fitted (not 1), so this script uses the generalized form
# `ceiling(max(get_depth(phy)) / time_scale / 10) * 10` instead --
# `node_time_to_present()` (LeafRank/R/tree_utils.R:112-131) confirms
# `time_scale` is a divisor applied to raw tree-branch depths to get internal
# model time, so T_vector must span that same divided-through range.

source("scDNA/5_LeafRank/config.R")

suppressPackageStartupMessages({
  library(ape)
  library(LeafRank)
})

log_msg <- function(fmt, ...) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), sprintf(fmt, ...)))
}

tree_file_ultra <- file.path(group_output_dir, paste0(group_name, "_ultrametric_tree.rds"))
if (!file.exists(tree_file_ultra)) {
  stop("Ultrametric tree not found: ", tree_file_ultra, " -- run rescale_tree.R first.")
}
phy <- readRDS(tree_file_ultra)
log_msg("Loaded ultrametric tree: %d tips", ape::Ntip(phy))

# Override config.R's placeholder n_sampled_cells/rho with the REAL tip count
# of the tree actually being run -- this must match what the calibration
# probe validated (it used the real 1571-tip count for R1, not the config
# placeholder), otherwise the fitted tau would differ from what was checked.
n_sampled_cells <- ape::Ntip(phy)
rho <- n_sampled_cells / total_population
log_msg("Recomputed rho from real tip count: rho = %.3e  (n_sampled_cells=%d, total_population=%.0e)",
        rho, n_sampled_cells, total_population)

### Step 4: finalize branching-process parameters ---------------------------
log_msg("Solving for tau via get_full_pars() (V=%d, rho=%.3e)", V, rho)
pars <- get_full_pars(b_rates = b_rates, d_rates = d_rates, nu = nu, rho = rho,
                       tree = phy, model = "default")
names(pars) <- c("b_rates", "d_rates", "nu", "rho", "tau")
time_scale <- pars$tau
log_msg("Fitted time_scale (tau) = %.6e  (implies internal model time T = 1/tau = %.2f)",
        time_scale, 1 / time_scale)

nu_vec <- rep(nu, V)

timeFrom <- 0
timeTo   <- ceiling(max(get_depth(phy)) / time_scale / 10) * 10
timeBy   <- timeTo / 200
T_vector <- seq(from = timeFrom, to = timeTo, by = timeBy)
log_msg("T_vector: %d points spanning [%g, %g]", length(T_vector), timeFrom, timeTo)

non_negativity_cutoff <- 0

out_file <- file.path(group_output_dir, paste0(group_name, "_pFitness.rds"))
dir.create(group_output_dir, recursive = TRUE, showWarnings = FALSE)

log_msg("Config: V=%d, d_t=%.3f, n_threads=%d", V, d_t, n_threads)
log_msg("Output will be written to: %s", out_file)

### Step 5: run LeafRank -----------------------------------------------------
log_msg("Starting LeafRank() -- this is the expensive step, see TIMING WARNING above")
t0 <- Sys.time()
result <- LeafRank(phy, out_file, rho, d_t, time_scale, b_rates, d_rates, nu_vec,
                    T_vector, non_negativity_cutoff, n_threads)
elapsed <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
log_msg("LeafRank() returned '%s' after %.1f min", result, elapsed)
log_msg("Results saved to: %s", out_file)
