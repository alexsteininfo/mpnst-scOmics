#!/bin/bash
# run_medicc2_full_bootstrap.sh
#
# Run MEDICC2 with --events on every bootstrap TSV produced by
# create_bootstrap_datasets.R. Each replicate TSV gets one full MEDICC2
# run including ancestral reconstruction, yielding integer branch lengths
# comparable to the final tree (unlike MEDICC2's internal bootstrap).
#
# Input directory structure expected (output of create_bootstrap_datasets.R):
#   <BOOTSTRAP_DIR>/<method>/<sample_name>/boot_001.tsv ... metadata.tsv
#
# Output directory structure (mirrors input):
#   <OUTPUT_DIR>/<method>/<sample_name>/<boot_NNN>/
#     boot_NNN_final_tree.new          <- used in step 3 burden analysis
#     boot_NNN_pairwise_distances.tsv
#     boot_NNN_copynumber.tsv
#     ...

# ---------------------------------------------------------------------------
# Configuration — adjust paths before running
# ---------------------------------------------------------------------------

MINSEG="1.5"
FILTER_LEN=$(awk "BEGIN {printf \"%d\", $MINSEG * 1000000}")

BOOTSTRAP_DIR="/srv/home/aste0033/projects/MPNST/scDNA/MEDICC2/bootstrap_custom/datasets"
OUTPUT_DIR="/srv/home/aste0033/projects/MPNST/scDNA/MEDICC2/bootstrap_custom/output/minseg_${MINSEG}"
METHOD="cell-wise"   # "chr-wise" or "cell-wise"

MEDICC2="micromamba run -n python medicc2"

# ---------------------------------------------------------------------------
# Run
# ---------------------------------------------------------------------------

for method_dir in "$BOOTSTRAP_DIR/$METHOD/"; do
    method=$(basename "$method_dir")

    for sample_dir in "$method_dir"*/; do
        sample=$(basename "$sample_dir")

        for tsv_file in "$sample_dir"boot_*.tsv; do
            rep=$(basename "$tsv_file" .tsv)
            out="$OUTPUT_DIR/$method/$sample/$rep"
            mkdir -p "$out"

            echo "Running MEDICC2: $method / $sample / $rep"

            $MEDICC2 \
                "$tsv_file" \
                "$out" \
                --input-type t \
                --filter-segment-length "$FILTER_LEN" \
                --n-cores 28 \
                --events
        done
    done
done
