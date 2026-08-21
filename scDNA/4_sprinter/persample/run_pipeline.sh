#!/usr/bin/env bash
# Run the full per-sample SPRINTER pipeline:
#   1. SPRINTER runs (12 jobs: 6 regions × 2 technologies)
#   2. Summarize results into sphase_summary_persample.tsv
#   3. Plot sphase_fractions_persample.png
#
# Run from the repo root:
#   bash scDNA/4_sprinter/persample/run_pipeline.sh

set -euo pipefail

REPO_ROOT="$(git -C "$(dirname "$0")" rev-parse --show-toplevel)"
cd "$REPO_ROOT"

mkdir -p scDNA/4_sprinter/results

bash scDNA/4_sprinter/persample/run_sprinter_persample.sh

micromamba run -n R Rscript scDNA/4_sprinter/persample/summarize_sphase_persample.R

micromamba run -n R Rscript scDNA/4_sprinter/persample/plot_sphase_persample.R
