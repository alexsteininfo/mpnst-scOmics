#!/usr/bin/env bash
# Run SPRINTER on FIT208 DLP+ treating each sample as a single clone.
#
# All cells in a sample are assigned to clone "1" via `dev --fixclones`,
# bypassing SPRINTER's copy-number clustering step entirely.
# This yields a sample-level (not clone-level) S-phase fraction.
#
# Requires the hg38 SPRINTER resources to be generated first by running:
#   bash liftover_dlpp.sh
#
# Input:  $CHISEL_DIR/<sample>/rdr.tsv   (from chisel_rdr_dlpp.sh)
# Output: $OUTDIR/<sample>/sprinter.output.tsv.gz

set -euo pipefail

# ── Paths ──────────────────────────────────────────────────────────────────────

CHISEL_DIR="/mnt/iribhm/homes/aste0033/projects/MPNST/scDNA/CHISEL/DLPp"
OUTDIR="/srv/home/aste0033/projects/MPNST/scDNA/SPRINTER/DLPp_oneclone"
RESOURCE_DIR="/srv/home/aste0033/projects/MPNST/scDNA/SPRINTER/resources_hg38"

REF="/mnt/iribhm/genomes/hg38-ngs/hg38.fa"

SAMPLES=(R1 R2 R3 R4 R5 P)
THREADS=12

# ── SPRINTER parameters ────────────────────────────────────────────────────────

MIN_READS=100000   # minimum mapped reads per cell
RT_READS=200       # target reads/bin for replication-timing analysis
CN_READS=1000      # target reads/bin for copy-number analysis
MAX_PLOIDY=4       # raise to 5-6 if whole-genome doubling is suspected

# ── Environment ────────────────────────────────────────────────────────────────

eval "$(micromamba shell hook --shell bash)"
micromamba activate sprinter

log() { echo "[$(date '+%H:%M:%S')] $*"; }
die() { echo "[$(date '+%H:%M:%S')] ERROR: $*" >&2; exit 1; }

[[ -f "$REF" ]] || die "Reference FASTA not found: $REF"

RTSCORES_HG38="${RESOURCE_DIR}/rtscores_hg38.csv.gz"
GAPS_HG38="${RESOURCE_DIR}/gaps_hg38.tsv"

[[ -f "$RTSCORES_HG38" ]] || die "hg38 RT scores not found — run liftover_dlpp.sh first"
[[ -f "$GAPS_HG38"     ]] || die "hg38 gaps file not found  — run liftover_dlpp.sh first"

# ── Per-sample SPRINTER run ────────────────────────────────────────────────────

for sample in "${SAMPLES[@]}"; do

    rdr="${CHISEL_DIR}/${sample}/rdr.tsv"
    sample_outdir="${OUTDIR}/${sample}"
    out="${sample_outdir}/sprinter.output.tsv.gz"

    if [[ ! -f "$rdr" ]]; then
        log "${sample}: rdr.tsv not found — skipping (run chisel_rdr_dlpp.sh first)"
        continue
    fi
    if [[ ! -s "$rdr" ]]; then
        log "${sample}: rdr.tsv is empty — skipping (chisel_rdr may have failed)"
        continue
    fi
    if [[ -f "$out" ]]; then
        log "${sample}: output already exists — skipping"
        continue
    fi

    n_cells=$(awk 'NR>1{print $4}' "$rdr" | sort -u | wc -l)
    log "${sample}: starting SPRINTER — ${n_cells} cells assigned to single clone"

    mkdir -p "$sample_outdir"

    # Assign every cell in this sample to clone "1"
    fixclones_tsv="${sample_outdir}/fixclones.tsv"
    awk 'NR>1{print $4}' "$rdr" | sort -u | awk 'BEGIN{print "CELL\tCLONE"} {print $0"\t1"}' \
        > "$fixclones_tsv"

    input_gz="${sample_outdir}/input.tsv.gz"
    gzip -c "$rdr" > "$input_gz"

    cd "$sample_outdir"

    sprinter "$input_gz" \
        --refgenome     "$REF" \
        --minreads      "$MIN_READS" \
        --rtreads       "$RT_READS" \
        --cnreads       "$CN_READS" \
        --maxploidy     "$MAX_PLOIDY" \
        -j              "$THREADS" \
        -o              "$out" \
        dev \
        --fixclones  "$fixclones_tsv" \
        --rtscores   "$RTSCORES_HG38" \
        --gapsfile   "$GAPS_HG38"

    rm -f "$input_gz"
    log "${sample}: done → ${out}"

done

log "=== All DLP+ SPRINTER (single-clone) runs complete ==="
log "Output: ${OUTDIR}/<sample>/sprinter.output.tsv.gz"
