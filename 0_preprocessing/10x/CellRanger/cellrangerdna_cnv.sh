#!/usr/bin/env bash
# Run Cell Ranger DNA CNV pipeline for all FIT208 10x scDNA samples.
#
# Cell Ranger DNA does in one command:
#   1. Barcode extraction and whitelist correction (737K-crdna-v1.txt)
#   2. BWA mem alignment to hg38
#   3. Duplicate marking
#   4. Cell calling (knee-point on corrected barcode counts)
#   5. CNV calling per cell
#
# Output per sample in OUTDIR/<sample>/outs/:
#   possorted_bam.bam           — position-sorted BAM with corrected CB tags
#   per_cell_summary_metrics.csv — per-cell metrics including cell/non-cell flag
#   cnv_data.h5                 — CNV profiles
#   web_summary.html            — run QC report

set -euo pipefail

# ── Environment ───────────────────────────────────────────────────────────────

source /srv/home/aste0033/CellRanger-DNA/setup_env.sh

# ── Parameters ────────────────────────────────────────────────────────────────

FASTQ_BASE="/mnt/iribhm/people/mtarabichi/MPNST/pvanloo_20250618/Haixi/10x_DNA"
REF="/srv/home/aste0033/projects/MPNST/scDNA/10X/reference/refdata-hg38"
OUTDIR="/srv/home/aste0033/projects/MPNST/scDNA/10X/cellranger_output"

CORES=32
MEM=256   # GB

# ── Helper ────────────────────────────────────────────────────────────────────

log() { echo "[$(date '+%H:%M:%S')] $*"; }
die() { echo "[$(date '+%H:%M:%S')] ERROR: $*" >&2; exit 1; }

run_sample() {
    local sample=$1
    local fastqs=$2   # comma-separated paths for multi-run samples

    local sample_outdir="${OUTDIR}/${sample}"
    local outs="${sample_outdir}/outs/possorted_bam.bam"

    if [[ -f "$outs" ]]; then
        log "${sample}: possorted_bam.bam already exists — skipping"
        return 0
    fi

    log "${sample}: starting Cell Ranger DNA CNV"
    log "  fastqs: ${fastqs}"

    # Cell Ranger creates <id>/ in the current directory
    cd "$OUTDIR"

    cellranger-dna cnv \
        --id="${sample}" \
        --fastqs="${fastqs}" \
        --sample="${sample}" \
        --reference="${REF}" \
        --localcores="${CORES}" \
        --localmem="${MEM}"

    [[ -f "$outs" ]] || die "${sample}: expected BAM not found after run: ${outs}"
    log "${sample}: done  →  ${outs}"
}

# ── Checks ────────────────────────────────────────────────────────────────────

[[ -d "$REF" ]] || die "Reference not found: ${REF}  (run build_cellrangerdna_ref.sh first)"
mkdir -p "$OUTDIR"

# ── Sample runs ───────────────────────────────────────────────────────────────
# FIT208A3 and FIT208A4 were sequenced once (SC18143_1 only).
# FIT208A5–A8 were re-sequenced across three runs; all three paths are
# passed as a comma-separated list so Cell Ranger merges them automatically.

RUN1="${FASTQ_BASE}/SC18143_1"
RUN2="${FASTQ_BASE}/SC18143_2"
RUN2B="${FASTQ_BASE}/sc18143_2b"

run_sample FIT208A3 "${RUN1}"
run_sample FIT208A4 "${RUN1}"

COMBO="${RUN1},${RUN2},${RUN2B}"
run_sample FIT208A5 "${COMBO}"
run_sample FIT208A6 "${COMBO}"
run_sample FIT208A7 "${COMBO}"
run_sample FIT208A8 "${COMBO}"

# ── Summary ───────────────────────────────────────────────────────────────────

log "=== All samples complete ==="
for sample in FIT208A3 FIT208A4 FIT208A5 FIT208A6 FIT208A7 FIT208A8; do
    metrics="${OUTDIR}/${sample}/outs/per_cell_summary_metrics.csv"
    if [[ -f "$metrics" ]]; then
        n_cells=$(awk -F',' 'NR>1 {n++} END {print n+0}' "$metrics")
        log "  ${sample}: ${n_cells} cells called"
    else
        log "  ${sample}: metrics file missing — check logs under ${OUTDIR}/${sample}/"
    fi
done
