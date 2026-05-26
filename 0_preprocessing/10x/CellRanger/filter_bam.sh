#!/usr/bin/env bash
# Filter Cell Ranger DNA BAMs to reads whose CB:Z: tag is in the per-sample
# knee-point barcode whitelist (valid_barcodes_knee_wl.txt).
#
# Source:  $CELLRANGER_DIR/<sample>/outs/possorted_bam.bam
# Barcode: $BWA_DIR/<sample>/valid_barcodes_knee_wl.txt
# Output:  $OUTDIR/<sample>/possorted_bam.bam  (+.bai index)
#
# $OUTDIR is a new sibling folder next to cellranger_output/ and bwa_output/.
# The original CellRanger BAMs are not modified.
#
# Uses samtools -D CB:FILE to filter by tag whitelist (requires samtools >= 1.12).

set -euo pipefail

module load SAMtools/1.18-GCC-12.3.0

# ── Paths ──────────────────────────────────────────────────────────────────────

BASE="/srv/home/aste0033/projects/MPNST/scDNA/10X"

CELLRANGER_DIR="${BASE}/cellranger_output"   # source BAMs (Cell Ranger DNA output)
BWA_DIR="${BASE}/bwa_output"                 # barcode whitelists live here
OUTDIR="${BASE}/bwa_filtered"                # filtered BAMs written here

THREADS=12

SAMPLES=(FIT208A3 FIT208A4 FIT208A5 FIT208A6 FIT208A7 FIT208A8)

# ── Helpers ────────────────────────────────────────────────────────────────────

log() { echo "[$(date '+%H:%M:%S')] $*"; }
die() { echo "[$(date '+%H:%M:%S')] ERROR: $*" >&2; exit 1; }

# ── Per-sample filter ──────────────────────────────────────────────────────────

for sample in "${SAMPLES[@]}"; do

    src_bam="${CELLRANGER_DIR}/${sample}/outs/possorted_bam.bam"
    barcodes="${BWA_DIR}/${sample}/valid_barcodes_knee_wl.txt"
    sample_outdir="${OUTDIR}/${sample}"
    out_bam="${sample_outdir}/possorted_bam.bam"

    if [[ ! -f "$src_bam" ]]; then
        log "${sample}: source BAM not found — skipping"
        continue
    fi
    if [[ ! -f "$barcodes" ]]; then
        log "${sample}: barcode whitelist not found — skipping"
        continue
    fi
    if [[ -f "$out_bam" ]]; then
        log "${sample}: filtered BAM already exists — skipping"
        continue
    fi

    n_barcodes=$(wc -l < "$barcodes")
    log "${sample}: filtering to ${n_barcodes} whitelisted cells"

    mkdir -p "$sample_outdir"

    # Keep only reads whose CB:Z: tag matches a barcode in the whitelist.
    # -D CB:FILE filters by tag value; -bh preserves the header in BAM output.
    samtools view \
        -bh \
        -@ "$THREADS" \
        -D "CB:${barcodes}" \
        -o "$out_bam" \
        "$src_bam"

    samtools index -@ "$THREADS" "$out_bam"

    n_reads=$(samtools view -c -@ "$THREADS" "$out_bam")
    log "${sample}: done — ${n_reads} reads  →  ${out_bam}"

done

log "=== All samples complete ==="
log "Filtered BAMs: ${OUTDIR}/<sample>/possorted_bam.bam"
log "Update BAM_BASE in chisel_rdr_10x.sh to: ${OUTDIR}"
