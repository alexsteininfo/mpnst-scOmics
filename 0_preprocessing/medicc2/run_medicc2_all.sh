#!/bin/bash

# Set minimum segment length (1.5 Mb or 2.5 Mb following original analysis)
MINSEG="1.5"
FILTER_LEN=$(awk "BEGIN {printf \"%d\", $MINSEG * 1000000}")

# Define input and output directories
INPUT_FILE="/mnt/iribhm/homes/aste0033/projects/CNA-PhyloAnalysis/Haixi/scDNA/MEDICC2_input/all_10X_DLP/minseg_${MINSEG}/MPNST_all_single_cells_${MINSEG}_ds_6039_seed100.tsv"
OUTPUT_DIR="/mnt/iribhm/homes/aste0033/projects/CNA-PhyloAnalysis/MEDICC2/10X_DLP/all/minseg_${MINSEG}"

mkdir -p "$OUTPUT_DIR"

medicc2 \
    "$INPUT_FILE" \
    "$OUTPUT_DIR" \
    --input-type t \
    --filter-segment-length "$FILTER_LEN" \
    --n-cores 24 \
    --events
