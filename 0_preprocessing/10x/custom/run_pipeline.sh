#!/bin/bash
# End-to-end pipeline: test run → full alignment → cell calling.
# Steps run in strict order; any failure aborts the pipeline.
set -euo pipefail

SCRIPTS="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUTDIR="/srv/home/aste0033/projects/MPNST/scDNA/10X/bwa_output"
TEST_OUTDIR="/srv/home/aste0033/projects/MPNST/scDNA/10X/bwa_output_test"

log() { echo "[$(date '+%H:%M:%S')] $*"; }
die() { echo "[$(date '+%H:%M:%S')] ERROR: $*" >&2; exit 1; }

# ── Step 1: Test run (one lane of FIT208A3) ───────────────────────────────────
log "=== Step 1: Test run (FIT208A3, one lane) ==="
bash "${SCRIPTS}/fastq_to_bam_bwa_test.sh"

# Validate test output
TEST_BAM="${TEST_OUTDIR}/FIT208A3/markdup.bam"
[[ -f "$TEST_BAM" ]] || die "Test BAM not found: ${TEST_BAM}"
log "Test BAM exists: ${TEST_BAM}"

KEEP_RATE=$(grep -oP 'Kept\s+:\s+\K[0-9]+' "${TEST_OUTDIR}/FIT208A3/extract_barcodes.log" 2>/dev/null \
            | tail -1 || echo "")
log "Barcode keep count: ${KEEP_RATE:-unknown}"

# In raw mode every read is kept, so the 'Kept' line should equal total reads.
# In whitelist mode, expect ~25-30%. Either way, abort if log is empty.
[[ -s "${TEST_OUTDIR}/FIT208A3/extract_barcodes.log" ]] \
    || die "extract_barcodes.log is empty — extraction may have failed."

log "Test run passed. Proceeding to full run."

# ── Step 2: Full pipeline — all 6 samples (alignment + markdup only) ──────────
log "=== Step 2: Full pipeline (all 6 samples, alignment + markdup) ==="
bash "${SCRIPTS}/fastq_to_bam_bwa.sh"

# Validate all 6 markdup BAMs exist (BQSR is a separate step via run_bqsr_only.sh)
log "Validating output BAMs..."
for sample in FIT208A3 FIT208A4 FIT208A5 FIT208A6 FIT208A7 FIT208A8; do
    bam="${OUTDIR}/${sample}/markdup.bam"
    [[ -f "$bam" ]] || die "Missing markdup BAM: ${bam}"
    log "  ${sample}: OK"
done

# ── Step 3: Extract per-barcode read counts ────────────────────────────────────
log "=== Step 3: Extract barcode counts ==="
bash "${SCRIPTS}/extract_barcode_counts.sh"

for sample in FIT208A3 FIT208A4 FIT208A5 FIT208A6 FIT208A7 FIT208A8; do
    counts="${OUTDIR}/${sample}/barcode_counts.txt"
    [[ -s "$counts" ]] || die "Missing or empty barcode_counts.txt for ${sample}"
done
log "Barcode counts ready for all 6 samples."

# ── Step 4: Knee-point cell calling ───────────────────────────────────────────
log "=== Step 4: Knee-point cell calling ==="
eval "$(micromamba shell hook --shell bash)"
micromamba activate R
Rscript "${SCRIPTS}/knee_point_cell_calling.R"

# ── Step 5: Cell Ranger DNA barcodes + comparison ─────────────────────────────
log "=== Step 5: Cell Ranger DNA barcode extraction and comparison ==="
Rscript "${SCRIPTS}/extract_barcodes_from_rds.R"

log "=== Pipeline complete ==="
log "Check: ${OUTDIR}/cell_calling_comparison.tsv for Jaccard scores."
log "If Jaccard >= 0.8 for all samples, the pipeline is working correctly."
log "Run run_bqsr_only.sh when you are ready to apply BQSR (produces possorted_bam.bam)."
