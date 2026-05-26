#!/usr/bin/env bash
# Test run of Cell Ranger DNA CNV pipeline on a subsampled FIT208A3 dataset.
# Identical to cellrangerdna_cnv.sh except:
#   - Only FIT208A3 is processed (single run, simplest sample)
#   - FASTQs are subsampled to SUBSAMPLE_READS reads per file (~5% of real data)
#   - Output goes to a separate test directory
#
# Purpose: catch pipeline/API errors quickly without waiting 24+ hours.
# Runtime: ~30–60 min depending on cluster load.

set -euo pipefail

# ── Environment ───────────────────────────────────────────────────────────────

source /srv/home/aste0033/CellRanger-DNA/setup_env.sh

# ── Parameters ────────────────────────────────────────────────────────────────

FASTQ_BASE="/mnt/iribhm/people/mtarabichi/MPNST/pvanloo_20250618/Haixi/10x_DNA"
REF="/srv/home/aste0033/projects/MPNST/scDNA/10X/reference/refdata-hg38"
OUTDIR="/srv/home/aste0033/projects/MPNST/scDNA/10X/cellranger_test"
TEST_FASTQ_DIR="${OUTDIR}/test_fastqs"

CORES=32
MEM=256

# Number of reads to keep per FASTQ file (4 lines per read in FASTQ format).
# 500 000 reads × 3 lanes = 1.5 M total reads → ~675 reads/cell for FIT208A3
# (~2 220 cells). Below ~200 reads/cell the barcode knee-point detection fails
# (100% incorrect barcodes error from TRIM_READS).
SUBSAMPLE_READS=500000

# ── Helpers ───────────────────────────────────────────────────────────────────

log() { echo "[$(date '+%H:%M:%S')] $*"; }
die() { echo "[$(date '+%H:%M:%S')] ERROR: $*" >&2; exit 1; }

# ── Step 1: subsample FASTQs ──────────────────────────────────────────────────

[[ -d "$REF" ]] || die "Reference not found: ${REF}"
mkdir -p "$TEST_FASTQ_DIR"

SRC_DIR="${FASTQ_BASE}/SC18143_1"
MAX_LINES=$(( SUBSAMPLE_READS * 4 ))

# Always wipe and regenerate — avoids any stale/corrupt files from prior runs.
rm -rf "$TEST_FASTQ_DIR"
mkdir -p "$TEST_FASTQ_DIR"

for src in "${SRC_DIR}"/FIT208A3_*.fastq.gz; do
    fname=$(basename "$src")
    dst="${TEST_FASTQ_DIR}/${fname}"
    log "subsample: ${fname}  (${SUBSAMPLE_READS} reads)"
    # awk reads to EOF so zcat exits cleanly — no SIGPIPE, no broken gzip flush.
    zcat "$src" | awk -v n="$MAX_LINES" 'NR<=n' | gzip -c > "$dst"
    actual=$(zcat "$dst" | wc -l)
    log "  → ${actual} lines written ($(( actual / 4 )) reads)"
done

# ── Step 2: run Cell Ranger DNA ───────────────────────────────────────────────
# Always wipe the pipestance too — CellRanger resumes from cached chunk outputs
# even if the input FASTQs changed, so stale chunk results would persist otherwise.

mkdir -p "$OUTDIR"
rm -rf "${OUTDIR}/FIT208A3"

outs="${OUTDIR}/FIT208A3/outs/possorted_bam.bam"
if [[ -f "$outs" ]]; then
    log "FIT208A3: possorted_bam.bam already exists — skipping"
else
    log "FIT208A3: starting Cell Ranger DNA CNV (test, ${SUBSAMPLE_READS} reads/file)"
    cd "$OUTDIR"
    cellranger-dna cnv \
        --id="FIT208A3" \
        --fastqs="${TEST_FASTQ_DIR}" \
        --sample="FIT208A3" \
        --reference="${REF}" \
        --localcores="${CORES}" \
        --localmem="${MEM}"
    [[ -f "$outs" ]] || die "FIT208A3: expected BAM not found: ${outs}"
    log "FIT208A3: done  →  ${outs}"
fi

log "=== Test run complete ==="
log "Output: ${OUTDIR}/FIT208A3/outs/"
