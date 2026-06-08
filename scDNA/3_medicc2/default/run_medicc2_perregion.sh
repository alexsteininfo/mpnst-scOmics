#!/bin/bash

# Set minimum segment length (1.5 Mb or 2.5 Mb following original analysis)
MINSEG="2.5"
FILTER_LEN=$(awk "BEGIN {printf \"%d\", $MINSEG * 1000000}")

# Define input and output directories
INPUT_DIR="/mnt/iribhm/homes/aste0033/projects/CNA-PhyloAnalysis/Haixi/scDNA/MEDICC2_input/per_region_10X_DLP/minseg_${MINSEG}"
OUTPUT_DIR="/mnt/iribhm/homes/aste0033/projects/CNA-PhyloAnalysis/MEDICC2/10X_DLP/per_region/minseg_${MINSEG}"

# Loop through each TSV file in the input directory and run MEDICC2
for tsv_file in "$INPUT_DIR"/*.tsv; do
    sample=$(basename "$tsv_file" .tsv)
    out="$OUTPUT_DIR/$sample"
    mkdir -p "$out"

    medicc2 \
        "$tsv_file" \
        "$out" \
        --input-type t \
        --filter-segment-length "$FILTER_LEN" \
        --n-cores 24 \
        --events
done
