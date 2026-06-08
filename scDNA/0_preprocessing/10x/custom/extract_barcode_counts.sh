#!/bin/bash
# Extract per-barcode read counts from BAMs for all FIT208 samples.
# Prefers possorted_bam.bam (post-BQSR); falls back to markdup.bam.
# BQSR only affects base quality scores, not barcode tags or alignments,
# so running on markdup.bam gives identical cell-calling results.
#
# Output: <OUTDIR>/<sample>/barcode_counts.txt
#         Two-column text file: read_count  barcode (sorted by count descending)
set -euo pipefail

module load SAMtools/1.18-GCC-12.3.0

OUTDIR="/srv/home/aste0033/projects/MPNST/scDNA/10X/bwa_output"

extract_counts() {
    local sample=$1
    local sampledir="${OUTDIR}/${sample}"
    local out="${sampledir}/barcode_counts.txt"

    # Prefer final BAM; fall back to markdup if BQSR not yet done
    local bam
    if [[ -f "${sampledir}/possorted_bam.bam" ]]; then
        bam="${sampledir}/possorted_bam.bam"
    elif [[ -f "${sampledir}/markdup.bam" ]]; then
        bam="${sampledir}/markdup.bam"
        echo "[${sample}] NOTE: using markdup.bam (BQSR pending — counts will be identical)"
    else
        echo "[${sample}] ERROR: no BAM found, skipping"
        return 1
    fi

    if [[ -f "$out" ]]; then
        echo "[${sample}] barcode_counts.txt already exists, skipping"
        return 0
    fi

    echo "[${sample}] Extracting barcode counts from $(basename "$bam")..."
    samtools view -F 4 "$bam" \
        | awk '{for(i=12;i<=NF;i++) if($i~/^CB:Z:/) {print substr($i,6); break}}' \
        | sort | uniq -c | sort -rn \
        > "$out"

    local n_barcodes
    n_barcodes=$(wc -l < "$out")
    echo "[${sample}] Done — ${n_barcodes} unique barcodes"
}

_pids=()
for sample in FIT208A3 FIT208A4 FIT208A5 FIT208A6 FIT208A7 FIT208A8; do
    extract_counts "$sample" &
    _pids+=($!)
done

for _pid in "${_pids[@]}"; do
    wait "$_pid" || echo "WARNING: job $_pid exited with error"
done

echo "Barcode count extraction complete."
