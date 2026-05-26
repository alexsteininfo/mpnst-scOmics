#!/bin/bash
# Barcode extraction mode:
#   RAW_MODE=false (default) — 1-mismatch whitelist correction; ~27% keep rate.
#   RAW_MODE=true            — keep all reads with raw R1[0:16] as CB tag; ~4× larger BAMs.
#
# Use RAW_MODE=false with the scDNA-specific whitelist (737K-crdna-v1.txt).
# The GEX v2 whitelist (737K-august-2016.txt) only has 25% overlap with the
# scDNA whitelist — using it discards ~70% of true cells. Always use the
# scDNA whitelist for Chromium Single Cell CNV data.
RAW_MODE=false

set -euo pipefail

# ── Modules ───────────────────────────────────────────────────────────────────
# BWA-MEME is run via Apptainer (no module needed — binary is self-contained in the SIF).
# GATK/BQSR is a separate step; run run_bqsr_only.sh when needed.
module load SAMtools/1.18-GCC-12.3.0
module load Python/3.12.3-GCCcore-13.3.0

# ── Paths ─────────────────────────────────────────────────────────────────────
FASTQ_BASE="/mnt/iribhm/people/mtarabichi/MPNST/pvanloo_20250618/Haixi/10x_DNA"
GENOME="/mnt/iribhm/genomes/hg38-gatk/Homo_sapiens_assembly38.fa"
OUTDIR="/srv/home/aste0033/projects/MPNST/scDNA/10X/bwa_output"
# BWA-MEME Singularity image (cluster-wide install, same path Mathieu uses)
BWAMEME_SIF="/mnt/iribhm/software/singularity/bwameme.sif"
SCRIPTS="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# ── Parallelism ───────────────────────────────────────────────────────────────
# Number of samples to process simultaneously.
# Cores are divided evenly across parallel jobs.
# Recommended: PARALLEL=3 → 8 cores/sample, wall time ~8 h instead of ~2–3 days.
# Memory budget at PARALLEL=3: ~3 × (70 GB BWA + 64 GB GATK) peaks well below 500 GB
# because BWA and GATK never overlap within the same sample.
PARALLEL=4
CORES_PER_JOB=$((32 / PARALLEL))   # 24 total cores divided across jobs

# ── Barcode whitelist ─────────────────────────────────────────────────────────
# Use the scDNA-specific whitelist (737K-crdna-v1.txt), NOT the GEX v2 whitelist
# (737K-august-2016.txt). The two lists share only ~25% of barcodes. Using the
# wrong one discards ~75% of true scDNA cells at the barcode-correction step.
WHITELIST="/srv/home/aste0033/projects/MPNST/scDNA/10X/737K-crdna-v1.txt"
if [[ "$RAW_MODE" == "true" ]]; then
    echo "Barcode mode:   RAW (all reads kept; whitelist not used)"
else
    if [[ ! -f "$WHITELIST" ]]; then
        echo "ERROR: scDNA whitelist not found: $WHITELIST"
        echo "Download it from https://github.com/10XGenomics/cellranger-dna"
        echo "(lib/python/cellranger/barcodes/737K-crdna-v1.txt.gz)"
        exit 1
    fi
    echo "Barcode mode:   whitelist ($WHITELIST)"
fi
echo "Parallel jobs:  $PARALLEL  (${CORES_PER_JOB} cores each)"

# ── Checks ────────────────────────────────────────────────────────────────────
for tool in apptainer samtools python3; do
    command -v "$tool" &>/dev/null || { echo "ERROR: $tool not in PATH"; exit 1; }
done
[[ -f "$BWAMEME_SIF" ]] || { echo "ERROR: BWA-MEME SIF not found: $BWAMEME_SIF"; exit 1; }
echo "BWA-MEME: $(apptainer exec --cleanenv "$BWAMEME_SIF" bwa-meme 2>&1 | grep -i 'version\|Program' | head -1 || true)"
echo "samtools: $(samtools --version | head -1)"

mkdir -p "$OUTDIR"

# ── Helper ────────────────────────────────────────────────────────────────────
process_sample() {
    local sample=$1
    shift
    local fastq_dirs=("$@")

    local sampledir="${OUTDIR}/${sample}"
    mkdir -p "$sampledir"

    # Collect R1 and R2 files across all sequencing runs, sorted per directory
    # so that R1[i] and R2[i] are always paired (same lane, same run).
    local r1_files=() r2_files=()
    for dir in "${fastq_dirs[@]}"; do
        while IFS= read -r f; do r1_files+=("$f"); done \
            < <(find "$dir" -maxdepth 1 -name "${sample}_*_R1_*.fastq.gz" | sort)
        while IFS= read -r f; do r2_files+=("$f"); done \
            < <(find "$dir" -maxdepth 1 -name "${sample}_*_R2_*.fastq.gz" | sort)
    done

    echo "[${sample}] ${#r1_files[@]} R1 files, ${#r2_files[@]} R2 files"

    if [[ ${#r1_files[@]} -eq 0 ]]; then
        echo "[${sample}] WARNING: no R1 files found, skipping"
        return
    fi

    local markdup_bam="${sampledir}/markdup.bam"
    local rg="@RG\tID:${sample}\tSM:${sample}\tPL:ILLUMINA\tLB:${sample}"
    local st_threads=$(( CORES_PER_JOB / 2 > 1 ? CORES_PER_JOB / 2 : 1 ))

    # ── Steps 1–4: align + markdup ────────────────────────────────────────────
    echo "[${sample}] Aligning and marking duplicates..."
    local bc_args=()
    if [[ "$RAW_MODE" == "true" ]]; then
        bc_args=(--raw)
    else
        bc_args=(--whitelist "$WHITELIST")
    fi
    python3 "${SCRIPTS}/extract_barcodes.py" \
        --r1 "${r1_files[@]}" \
        --r2 "${r2_files[@]}" \
        "${bc_args[@]}" \
        2>"${sampledir}/extract_barcodes.log" \
    | apptainer exec --cleanenv \
        --bind /mnt/iribhm:/mnt/iribhm \
        "$BWAMEME_SIF" \
        bwa-meme mem -t "$CORES_PER_JOB" -M -C -R "$rg" "$GENOME" - \
        2>"${sampledir}/bwa.log" \
    | samtools fixmate -m -u - - \
    | samtools sort -u -@ "$st_threads" -T "${sampledir}/sort_tmp" \
    | samtools markdup --barcode-tag CB -@ "$st_threads" - "$markdup_bam"

    samtools index "$markdup_bam"

    echo "[${sample}] Barcode stats:"
    tail -5 "${sampledir}/extract_barcodes.log"
    echo "[${sample}] Done → ${markdup_bam}  (run run_bqsr_only.sh for BQSR)"
}

# ── Parallel dispatch ─────────────────────────────────────────────────────────
# Submits samples as background jobs in batches of $PARALLEL.
# Each batch waits for all jobs to finish before the next batch starts.
_pids=()
_submit() {
    process_sample "$@" &
    _pids+=($!)
    if [[ ${#_pids[@]} -ge $PARALLEL ]]; then
        for _pid in "${_pids[@]}"; do
            wait "$_pid" || echo "WARNING: job $_pid exited with error"
        done
        _pids=()
    fi
}

# FIT208A3 (Primary/R1) and FIT208A4 (R2): only sequenced in SC18143_1
_submit FIT208A3 "${FASTQ_BASE}/SC18143_1"
_submit FIT208A4 "${FASTQ_BASE}/SC18143_1"

# FIT208A5–A8 (R3, R4, R5, P): combine all three sequencing runs for max depth
for sample in FIT208A5 FIT208A6 FIT208A7 FIT208A8; do
    _submit "${sample}" \
        "${FASTQ_BASE}/SC18143_1" \
        "${FASTQ_BASE}/SC18143_2" \
        "${FASTQ_BASE}/sc18143_2b"
done

# Wait for the last batch (fewer than $PARALLEL samples)
for _pid in "${_pids[@]}"; do
    wait "$_pid" || echo "WARNING: job $_pid exited with error"
done

echo "All samples complete."
