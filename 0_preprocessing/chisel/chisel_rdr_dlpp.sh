#!/usr/bin/env bash
# Run chisel_rdr on FIT208 DLP+ per-sample BAMs to produce 50 kb bin
# read-count profiles, which are used as input for SPRINTER.
#
# Cell identity is read from the RG:Z: tag (read group).
# Read group IDs follow the pattern: <chip>_<library>_<x>x_<y>y
# e.g. LES4677A1_124574_44x_26y
#
# Per-sample BAMs were produced by split_bam_by_sample.sh, which split each
# chip's merged.bam by the cell-to-region mapping in scDNA_DLP_fastq.csv:
#   LES4677A1 → R1.bam, R2.bam
#   LES4677A2 → R3.bam, R4.bam
#   LES4677A4 → R5.bam, P.bam
#
# Output per sample: $OUTDIR/<sample>/rdr.tsv
# Columns are renamed for SPRINTER: CHROMOSOME→CHR, NORMAL→NORM_COUNT, RDR→RAW_RDR
#
# Requires: micromamba environment "chisel" with chisel installed (Python 2)

set -euo pipefail

# ── Parameters ────────────────────────────────────────────────────────────────

# chip: directory under BAM_BASE where the per-sample BAM lives
declare -A CHIP=(
    ["R1"]="LES4677A1"
    ["R2"]="LES4677A1"
    ["R3"]="LES4677A2"
    ["R4"]="LES4677A2"
    ["R5"]="LES4677A4"
    ["P"]="LES4677A4"
)
SAMPLES=(R1 R2 R3 R4 R5 P)

BAM_BASE="/mnt/iribhm/homes/aste0033/projects/MPNST/scDNA/DLPp/bwa_output"
REF="/mnt/iribhm/genomes/hg38-ngs/hg38.fa"
OUTDIR="/mnt/iribhm/homes/aste0033/projects/MPNST/scDNA/CHISEL/DLPp"

BIN_SIZE=50000   # 50 kb — required to match SPRINTER rtscores.csv.gz
THREADS=12

# Minimum reads per read group to be considered a cell.
# DLP+ empty wells will have very few reads and are filtered out here.
# Adjust downward if genuine cells are being excluded.
MIN_READS=100000

# ── Environment ───────────────────────────────────────────────────────────────

eval "$(micromamba shell hook --shell bash)"
micromamba activate chisel

log() { echo "[$(date '+%H:%M:%S')] $*"; }
die() { echo "[$(date '+%H:%M:%S')] ERROR: $*" >&2; exit 1; }

[[ -f "$REF" ]] || die "Reference FASTA not found: $REF"

# ── Per-sample run ────────────────────────────────────────────────────────────

for sample in "${SAMPLES[@]}"; do

    chip="${CHIP[$sample]}"
    bam="${BAM_BASE}/${chip}/${sample}.bam"
    sample_outdir="${OUTDIR}/${sample}"
    rdr_out="${sample_outdir}/rdr.tsv"

    [[ -f "$bam" ]] || { log "${sample}: BAM not found (${bam}) — skipping"; continue; }

    if [[ -f "$rdr_out" ]]; then
        log "${sample}: rdr.tsv already exists — skipping"
        continue
    fi

    mkdir -p "$sample_outdir"

    n_rg=$(samtools view -H "$bam" | grep -c "^@RG" || true)
    log "${sample}: starting chisel_rdr (${n_rg} read groups, min ${MIN_READS} reads filter)"

    # chisel_rdr writes output to the current working directory
    cd "$sample_outdir"

    chisel_rdr \
        -t "$bam" \
        -r "$REF" \
        -b "$BIN_SIZE" \
        -m "$MIN_READS" \
        -j "$THREADS" \
        --cellprefix "RG:Z:"

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

log "=== All DLP+ samples complete ==="
log "Output: ${OUTDIR}/<sample>/rdr.tsv"
log "These files are ready for SPRINTER (columns: CHR START END CELL NORM_COUNT COUNT RAW_RDR)"
