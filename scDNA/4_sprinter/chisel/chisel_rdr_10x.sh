#!/usr/bin/env bash
# Run chisel_rdr on FIT208 10x scDNA BAMs to produce 50 kb bin read-count
# profiles, which are used as input for SPRINTER.
#
# Cell identity is read from the CB:Z: tag (corrected cell barcode).
# Valid cells are restricted to the whitelist-corrected knee-point barcodes
# (valid_barcodes_knee_wl.txt), bypassing the minimum-reads threshold filter
# so that low-coverage samples (e.g. FIT208A5) are not over-filtered.
#
# Output per sample: $OUTDIR/<sample>/rdr.tsv
# The three-column rename (CHROMOSOME→CHR, NORMAL→NORM_COUNT, RDR→RAW_RDR)
# needed for SPRINTER is done inline after each sample completes.
#
# Requires: micromamba environment "chisel" with chisel installed (Python 2)

set -euo pipefail

# ── Parameters ────────────────────────────────────────────────────────────────

BAM_BASE="/mnt/iribhm/homes/aste0033/projects/MPNST/scDNA/10X/bwa_output"
REF="/mnt/iribhm/genomes/hg38-ngs/hg38.fa"
OUTDIR="/mnt/iribhm/homes/aste0033/projects/MPNST/scDNA/CHISEL/10X"

BIN_SIZE=50000   # 50 kb — required to match SPRINTER rtscores.csv.gz
THREADS=12

SAMPLES=(FIT208A3 FIT208A4 FIT208A5 FIT208A6 FIT208A7 FIT208A8)

# ── Environment ───────────────────────────────────────────────────────────────

eval "$(micromamba shell hook --shell bash)"
micromamba activate chisel

log() { echo "[$(date '+%H:%M:%S')] $*"; }
die() { echo "[$(date '+%H:%M:%S')] ERROR: $*" >&2; exit 1; }

[[ -f "$REF" ]] || die "Reference FASTA not found: $REF"

# ── Per-sample run ────────────────────────────────────────────────────────────

for sample in "${SAMPLES[@]}"; do

    bam="${BAM_BASE}/${sample}/possorted_bam.bam"
    barcodes="${BAM_BASE}/${sample}/valid_barcodes_knee_wl.txt"
    sample_outdir="${OUTDIR}/${sample}"
    rdr_out="${sample_outdir}/rdr.tsv"

    [[ -f "$bam"      ]] || { log "${sample}: BAM not found — skipping"; continue; }
    [[ -f "$barcodes" ]] || { log "${sample}: barcodes file not found — skipping"; continue; }

    if [[ -f "$rdr_out" ]]; then
        log "${sample}: rdr.tsv already exists — skipping"
        continue
    fi

    mkdir -p "$sample_outdir"
    log "${sample}: starting chisel_rdr ($(wc -l < "$barcodes") barcodes)"

    # Pre-filter BAM to valid barcodes — chisel_rdr does not support --barcodes;
    # samtools -D CB:file keeps only reads whose CB tag value is in the whitelist.
    filtered_bam="${sample_outdir}/filtered.bam"
    if [[ ! -f "$filtered_bam" ]]; then
        log "${sample}: filtering BAM to valid barcodes"
        samtools view -b -D CB:"$barcodes" "$bam" > "$filtered_bam"
        samtools index "$filtered_bam"
    fi

    # chisel_rdr writes output to the current working directory
    cd "$sample_outdir"

    chisel_rdr \
        -t "$filtered_bam" \
        -r "$REF" \
        -b "$BIN_SIZE" \
        -j "$THREADS" \
        --cellprefix "CB:Z:"

    # Rename columns for SPRINTER compatibility
    # CHROMOSOME → CHR, NORMAL → NORM_COUNT, RDR → RAW_RDR
    if [[ -f "${sample_outdir}/rdr/rdr.tsv" ]]; then
        awk 'BEGIN{OFS="\t"} NR==1{
            for(i=1;i<=NF;i++){
                if($i=="CHROMOSOME") $i="CHR"
                if($i=="NORMAL")     $i="NORM_COUNT"
                if($i=="RDR")        $i="RAW_RDR"
            }
        } 1' "${sample_outdir}/rdr/rdr.tsv" > "$rdr_out"
        log "${sample}: done  →  ${rdr_out}"
    else
        log "${sample}: WARNING — expected rdr/rdr.tsv not found; check chisel output dir"
    fi

done

log "=== All 10x samples complete ==="
log "Output: ${OUTDIR}/<sample>/rdr.tsv"
log "These files are ready for SPRINTER (columns: CHR START END CELL NORM_COUNT COUNT RAW_RDR)"
