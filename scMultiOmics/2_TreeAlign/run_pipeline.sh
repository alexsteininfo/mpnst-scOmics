#!/usr/bin/env bash
# Full TreeAlign pipeline: scDNA/scRNA input prep -> both TreeAlign modes ->
# comparison. Runs in TREEALIGN_MODE="test" by default (region R1 / chr21
# only); set TREEALIGN_MODE=full to run the full cohort once the test pass
# looks right (see scMultiOmics/2_TreeAlign/README.md for what each mode
# covers, and the exact `treealign` micromamba env creation commands).
#
# Run from the repo root:
#   bash scMultiOmics/2_TreeAlign/run_pipeline.sh                        # test
#   TREEALIGN_MODE=full bash scMultiOmics/2_TreeAlign/run_pipeline.sh    # full

set -euo pipefail

log() { echo "[$(date '+%H:%M:%S')] $*"; }
die() { echo "[$(date '+%H:%M:%S')] ERROR: $*" >&2; exit 1; }

REPO_ROOT="$(git -C "$(dirname "$0")" rev-parse --show-toplevel)"
cd "$REPO_ROOT"

export TREEALIGN_MODE="${TREEALIGN_MODE:-test}"
[[ "$TREEALIGN_MODE" == "test" || "$TREEALIGN_MODE" == "full" ]] \
  || die "TREEALIGN_MODE must be 'test' or 'full', got '$TREEALIGN_MODE'"

log "TREEALIGN_MODE=$TREEALIGN_MODE"
mkdir -p scMultiOmics/2_TreeAlign/results

log "Step 0: scDNA gene x cell CN matrix, clone table, tree"
micromamba run -n R Rscript scMultiOmics/2_TreeAlign/0_prepare_scdna_cn.R

log "Step 1: scRNA raw expression export"
singularity exec /mnt/iribhm/software/singularity/single-cell.sif \
  Rscript scMultiOmics/2_TreeAlign/1_prepare_scrna_expr.R

if ! micromamba env list | grep -qE '(^|[/ ])treealign( |$)'; then
  die "micromamba env 'treealign' not found. Create it first — see the" \
      "'Environment setup' commands in scMultiOmics/2_TreeAlign/README.md."
fi

log "Step 2: TreeAlign clone-label mode"
micromamba run -n treealign python scMultiOmics/2_TreeAlign/2_run_treealign_clonelabel.py

log "Step 3: TreeAlign tree mode"
micromamba run -n treealign python scMultiOmics/2_TreeAlign/3_run_treealign_tree.py

log "Step 4: compare clone-label mode, tree mode, and published scRNA cell_type"
micromamba run -n R Rscript scMultiOmics/2_TreeAlign/4_compare_results.R

log "Done ($TREEALIGN_MODE)."
