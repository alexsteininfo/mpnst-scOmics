#!/bin/bash
# Test run: FIT208A3 (Primary/R1, ~2,220 cells) from SC18143_1 only.
# Use this to validate the full pipeline before running all 6 samples.
#
# To make it even faster, limit to one lane by setting TEST_LANE below
# (e.g. TEST_LANE="L004"). Leave empty to use all lanes.
TEST_LANE="L004"

# Mirror the RAW_MODE setting from fastq_to_bam_bwa.sh.
# true  = keep all reads with raw barcode (recovers CR DNA cells, ~4× larger BAM)
# false = 1-mismatch whitelist correction only
RAW_MODE=true

set -euo pipefail

# ── Modules ───────────────────────────────────────────────────────────────────
# BWA-MEME is run via Apptainer (no module needed — binary is self-contained in the SIF).
# GATK/BQSR is a separate step; run run_bqsr_only.sh when needed.
module load SAMtools/1.18-GCC-12.3.0
module load Python/3.12.3-GCCcore-13.3.0

# ── Paths ─────────────────────────────────────────────────────────────────────
FASTQ_BASE="/mnt/iribhm/people/mtarabichi/MPNST/pvanloo_20250618/Haixi/10x_DNA"
GENOME="/mnt/iribhm/genomes/hg38-gatk/Homo_sapiens_assembly38.fa"
OUTDIR="/srv/home/aste0033/projects/MPNST/scDNA/10X/bwa_output_test"
SCRIPTS="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CORES=24
BWAMEME_SIF="/mnt/iribhm/software/singularity/bwameme.sif"
# ── Barcode whitelist ─────────────────────────────────────────────────────────
# scDNA-specific whitelist; only ~25% overlap with the GEX v2 whitelist.
WHITELIST="/srv/home/aste0033/projects/MPNST/scDNA/10X/737K-crdna-v1.txt"
if [[ "$RAW_MODE" == "true" ]]; then
    echo "Barcode mode:   RAW (all reads kept; whitelist not used)"
else
    if [[ ! -f "$WHITELIST" ]]; then
        echo "ERROR: scDNA whitelist not found: $WHITELIST"
        exit 1
    fi
    echo "Barcode mode:   whitelist ($WHITELIST)"
fi

# ── Checks ────────────────────────────────────────────────────────────────────
for tool in apptainer samtools python3; do
    command -v "$tool" &>/dev/null || { echo "ERROR: $tool not in PATH"; exit 1; }
done
[[ -f "$BWAMEME_SIF" ]] || { echo "ERROR: BWA-MEME SIF not found: $BWAMEME_SIF"; exit 1; }
echo "BWA-MEME: $(apptainer exec --cleanenv "$BWAMEME_SIF" bwa-meme 2>&1 | grep -i 'version\|Program' | head -1 || true)"
echo "samtools: $(samtools --version | head -1)"

mkdir -p "$OUTDIR"

# ── Test sample: FIT208A3, SC18143_1 only ────────────────────────────────────
sample="FIT208A3"
sampledir="${OUTDIR}/${sample}"
mkdir -p "$sampledir"

# Find R1/R2 files, optionally filtered to a single lane
lane_pattern="${TEST_LANE:+*${TEST_LANE}*}"
lane_pattern="${lane_pattern:-*}"

r1_files=()
r2_files=()
while IFS= read -r f; do r1_files+=("$f"); done \
    < <(find "${FASTQ_BASE}/SC18143_1" -maxdepth 1 \
        -name "${sample}_${lane_pattern}_R1_*.fastq.gz" | sort)
while IFS= read -r f; do r2_files+=("$f"); done \
    < <(find "${FASTQ_BASE}/SC18143_1" -maxdepth 1 \
        -name "${sample}_${lane_pattern}_R2_*.fastq.gz" | sort)

if [[ ${#r1_files[@]} -eq 0 ]]; then
    echo "ERROR: no R1 files found for ${sample} in SC18143_1"
    exit 1
fi

echo "[${sample}] ${#r1_files[@]} R1 files, ${#r2_files[@]} R2 files"
printf '  %s\n' "${r1_files[@]}"

markdup_bam="${sampledir}/markdup.bam"
rg="@RG\tID:${sample}\tSM:${sample}\tPL:ILLUMINA\tLB:${sample}"

# ── Steps 1–4: align + markdup ────────────────────────────────────────────────
echo "[${sample}] Aligning and marking duplicates..."
bc_args=()
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
    bwa-meme mem -t "$CORES" -M -C -R "$rg" "$GENOME" - \
    2>"${sampledir}/bwa.log" \
| samtools fixmate -m -u - - \
| samtools sort -u -@ 8 -T "${sampledir}/sort_tmp" \
| samtools markdup --barcode-tag CB -@ 8 - "$markdup_bam"

samtools index "$markdup_bam"

echo "[${sample}] Done → ${markdup_bam}  (run run_bqsr_only.sh for BQSR)"
echo ""
echo "── Verification ─────────────────────────────────────────────────────────"
echo "Flagstat:"
samtools flagstat "$markdup_bam"
echo ""
echo "Barcode tag present:"
samtools view "$markdup_bam" | head -1 | grep -o "CB:Z:[A-Z]*" || echo "WARNING: no CB tag found"
echo ""
echo "Read group:"
samtools view -H "$markdup_bam" | grep "^@RG"
