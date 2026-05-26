#!/bin/bash
# DLP+ single-cell FASTQ → per-sample deduplicated BAMs
#
# Aligns only the cells listed in scDNA_DLP_fastq.csv (~55% of FASTQs per chip;
# the rest are chip-border or QC-failed wells), then merges cells directly into
# per-sample BAMs. No intermediate per-chip merged.bam is produced.
#
#   Chip       Samples     Run directory
#   LES4677A1  R1, R2      220401_A01366_0166_AH2CNFDRX2
#   LES4677A2  R3, R4      220228_A01366_0152_AHTTMYDRXY
#   LES4677A4  R5, P       220909_A01366_0277_AHKF5CDRX2
#
# Pipeline per chip:
#   1. BWA-MEME per cell → samtools fixmate → sort → markdup  [per-cell BAM]
#      Markdup runs per-cell so reads from different cells at the same genomic
#      position are never treated as PCR duplicates of each other.
#   2. samtools merge, grouped by sample from CSV             [R1.bam … P.bam]
#   3. Cleanup of per-cell BAMs
#
# Output: <OUTDIR>/<chip>/<sample>.bam  (R1.bam, R2.bam, R3.bam, R4.bam, R5.bam, P.bam)
# BQSR is a separate step — run run_bqsr_only.sh when needed.
#
# Tags in final BAMs:
#   RG:Z:<cell_id>  — per-read read-group tag identifying the source cell
#                     (cell_id = chip position, e.g. LES4677A1_124574_14x_22y).
#   SM:<chip>       — sample field set to chip ID; not updated during merge.
#   No per-read CB or UB tags are added here.  If alleleCounter is needed for
#   downstream genotyping, add UB post-hoc with the awk approach in modbam.R.
#
# Parallelism:
#   CHIP_PARALLEL chips processed concurrently (each chip has ~4 400 cells).
#   CELL_PARALLEL cells aligned simultaneously within one chip (CELL_THREADS each).
#   A rolling pool keeps exactly CELL_PARALLEL jobs running at all times —
#   as soon as one cell finishes, the next starts (no fixed-batch idle time).
#   Total cores ≈ CHIP_PARALLEL × CELL_PARALLEL × CELL_THREADS.
#   Example: 1 × 16 × 2 = 32 cores (recommended: one chip at a time).
#
# NFS performance note:
#   Each BWA-MEME process loads the hg38 index (~4 GB) from /mnt/iribhm/ on
#   startup.  With 16 parallel cells this is 64 GB of NFS reads per batch.
#   Set LOCAL_GENOME_DIR to a fast local path (e.g. /tmp or /dev/shm) to copy
#   the index once and avoid repeated NFS I/O.  Leave empty to use GENOME as-is.

set -euo pipefail

# ── Modules ───────────────────────────────────────────────────────────────────
# BWA-MEME is run via Apptainer (no module needed — binary is self-contained in the SIF).
# GATK/BQSR is a separate step; run run_bqsr_only.sh when needed.
module load SAMtools/1.12-GCC-10.2.0

# ── Paths ─────────────────────────────────────────────────────────────────────
FASTQ_BASE="/mnt/iribhm/people/mtarabichi/MPNST/pvanloo_20250618/Haixi/DLP_plus"
GENOME="/mnt/iribhm/genomes/hg38-gatk/Homo_sapiens_assembly38.fa"
OUTDIR="/srv/home/aste0033/projects/MPNST/scDNA/DLPp/bwa_output"
BWAMEME_SIF="/mnt/iribhm/software/singularity/bwameme.sif"
CSV="/mnt/iribhm/people/mtarabichi/MPNST/pvanloo_20250618/Haixi/metadata/scDNA_DLP_fastq.csv"

# ── Parallelism ───────────────────────────────────────────────────────────────
CHIP_PARALLEL=1    # chips processed concurrently (recommend 1: each chip is ~4 400 cells)
CELL_PARALLEL=16   # concurrent per-cell alignments within one chip (rolling pool)
CELL_THREADS=2     # BWA-MEME threads per cell  →  16 × 2 = 32 cores (alignment phase)
BATCH_SIZE=500     # max BAMs per samtools merge call (stay within fd limits)
MERGE_THREADS=16   # threads for samtools merge (alignment done; all 32 cores free)
LOCAL_GENOME_DIR="/tmp" # genome index cached here to avoid NFS I/O on every cell startup

# ── Chip → run directory mapping ──────────────────────────────────────────────
# Each chip was sequenced in exactly one run; derived from DLP_metadata.R.
declare -A CHIP_TO_RUN=(
    ["LES4677A1"]="${FASTQ_BASE}/220401_A01366_0166_AH2CNFDRX2"  # R1, R2
    ["LES4677A2"]="${FASTQ_BASE}/220228_A01366_0152_AHTTMYDRXY"  # R3, R4
    ["LES4677A4"]="${FASTQ_BASE}/220909_A01366_0277_AHKF5CDRX2"  # R5, P
)

# ── Parse CSV → per-chip, per-region cell list files ─────────────────────────
# Written to a temp dir so background process_chip subshells can read them.
# cell_id extracted from fastq basename by stripping _S<N>_... suffix.
# chip = first _-delimited field of cell_id (e.g. "LES4677A1").
RGLIST_DIR=$(mktemp -d /tmp/dlpp_rglists_XXXXXX)
trap 'rm -rf "$RGLIST_DIR"' EXIT

[[ -f "$CSV" ]] || { echo "ERROR: CSV not found: $CSV"; exit 1; }

awk -F',' '
NR==1       { next }
$2 ~ /_R2_/ { next }
{
    split($2, p, "/"); fname = p[length(p)]
    sub(/_S[0-9]+_.*/, "", fname)
    split(fname, c, "_"); chip = c[1]
    print chip "_" $3 " " fname
}' "$CSV" | while IFS=' ' read -r key cell_id; do
    echo "$cell_id" >> "${RGLIST_DIR}/${key}.txt"
done

echo "Cell-to-region lists built in ${RGLIST_DIR}:"
for f in "${RGLIST_DIR}/"*.txt; do
    [[ -f "$f" ]] || continue
    printf "  %-30s  %d cells\n" "$(basename "$f")" "$(wc -l < "$f")"
done

# ── Validation ────────────────────────────────────────────────────────────────
for tool in apptainer samtools; do
    command -v "$tool" &>/dev/null || { echo "ERROR: $tool not in PATH"; exit 1; }
done
[[ -f "$BWAMEME_SIF" ]] || { echo "ERROR: BWA-MEME SIF not found: $BWAMEME_SIF"; exit 1; }
echo "BWA-MEME: $(apptainer exec --cleanenv "$BWAMEME_SIF" bwa-meme 2>&1 | grep -i 'version\|Program' | head -1 || true)"
echo "SAMtools: $(samtools --version | head -1)"

mkdir -p "$OUTDIR"

# ── Optional: cache genome index on local disk to avoid repeated NFS reads ────
# Each BWA-MEME process loads the full index on startup.  With 16 parallel cells
# this is ~64 GB of NFS traffic per minute.  Copying to /tmp or /dev/shm once
# (~4 GB, takes ~30 s locally) eliminates this bottleneck for the whole run.
ALIGN_GENOME="$GENOME"
if [[ -n "$LOCAL_GENOME_DIR" ]]; then
    mkdir -p "$LOCAL_GENOME_DIR"
    local_fa="${LOCAL_GENOME_DIR}/$(basename "$GENOME")"
    if [[ ! -f "${local_fa}.bwt.2bit.64" ]]; then
        echo "Copying genome index to ${LOCAL_GENOME_DIR} (one-time, ~20 GB)..."
        # Copy only the files BWA-MEME mode3 needs; skip the huge .rev_comp/
        # .pos_packed/.suffixarray files which are for other MEME modes.
        for ext in "" .0123 .bwt.2bit.64 .pac .ann .amb .fai; do
            cp "${GENOME}${ext}" "${LOCAL_GENOME_DIR}/"
        done
    fi
    ALIGN_GENOME="$local_fa"
    echo "Using local genome index: ${ALIGN_GENOME}"
fi

# ── Per-cell alignment ────────────────────────────────────────────────────────
# R1 and R2 paths are resolved once by process_chip() to avoid 4 400 separate
# find calls on the NFS run directory.
align_cell() {
    local cell_id=$1 chip=$2 celldir=$3 r1=$4 r2=$5
    local out_bam="${celldir}/${cell_id}.bam"

    [[ -f "$out_bam" ]] && return 0   # skip if already aligned (supports re-runs)

    local rg="@RG\tID:${cell_id}\tSM:${chip}\tPL:ILLUMINA\tLB:${chip}"

    apptainer exec --cleanenv \
        --bind /tmp:/tmp --bind /mnt/iribhm:/mnt/iribhm \
        "$BWAMEME_SIF" \
        bwa-meme mem -t "$CELL_THREADS" -M -R "$rg" "$ALIGN_GENOME" "$r1" "$r2" \
        2>"${celldir}/${cell_id}.bwa.log" \
    | samtools fixmate -m -u - - \
    | samtools sort -u -T "${celldir}/${cell_id}.sort_tmp" \
    | samtools markdup - "$out_bam" \
        2>"${celldir}/${cell_id}.markdup.log"
    # No per-cell index: samtools merge reads BAMs sequentially and never needs .bai
}

# ── Batch merge helper ────────────────────────────────────────────────────────
# Merges BAMs in batches of BATCH_SIZE to stay within OS fd limits.
# With ~2 200 cells per sample a flat merge would exceed ulimit -n (1024).
merge_bams() {
    local out_bam=$1; shift
    local bams=("$@")
    local n=${#bams[@]}

    if (( n <= BATCH_SIZE )); then
        samtools merge -f -@ "$MERGE_THREADS" "$out_bam" "${bams[@]}"
    else
        local tmpdir="${out_bam%.bam}_merge_tmp"
        mkdir -p "$tmpdir"
        local batch_bams=() i=0
        while (( i < n )); do
            local batch=("${bams[@]:$i:$BATCH_SIZE}")
            local batch_out="${tmpdir}/batch_${i}.bam"
            samtools merge -f -@ 8 "$batch_out" "${batch[@]}"
            batch_bams+=("$batch_out")
            i=$(( i + BATCH_SIZE ))
        done
        samtools merge -f -@ "$MERGE_THREADS" "$out_bam" "${batch_bams[@]}"
        rm -rf "$tmpdir"
    fi
}

# ── Per-chip pipeline ─────────────────────────────────────────────────────────
process_chip() {
    local chip=$1
    local rundir="${CHIP_TO_RUN[$chip]}"
    local chipdir="${OUTDIR}/${chip}"
    local celldir="${chipdir}/per_cell_bams"

    # Collect regions for this chip from the pre-parsed CSV lists
    local regions=()
    for f in "${RGLIST_DIR}/${chip}_"*.txt; do
        [[ -f "$f" ]] || continue
        local r; r=$(basename "$f" .txt); r="${r#${chip}_}"
        regions+=("$r")
    done
    [[ ${#regions[@]} -gt 0 ]] || { echo "[${chip}] ERROR: no regions found in CSV"; return 1; }

    # Skip chip entirely if all per-sample BAMs already exist
    local all_done=1
    for region in "${regions[@]}"; do
        [[ -f "${chipdir}/${region}.bam" ]] || { all_done=0; break; }
    done
    if [[ $all_done -eq 1 ]]; then
        echo "[${chip}] all per-sample BAMs exist — skipping"
        return 0
    fi

    mkdir -p "$celldir"

    # Build a hash set of CSV-valid cell IDs for this chip (O(1) lookup below).
    # Cells absent from the CSV are border/QC-failed wells — skip alignment entirely.
    declare -A valid_cells=()
    for f in "${RGLIST_DIR}/${chip}_"*.txt; do
        [[ -f "$f" ]] || continue
        while IFS= read -r cid; do valid_cells[$cid]=1; done < "$f"
    done

    # Scan the run directory once to build cell_id → R1/R2 paths.
    # This avoids one find() call per cell (4 400 NFS round-trips → 1).
    local -a cell_ids=() r1_paths=() r2_paths=()
    local n_fastq=0 n_skipped=0
    while IFS= read -r r1; do
        n_fastq=$(( n_fastq + 1 ))
        local cid r2
        cid=$(basename "$r1" | sed 's/_S[0-9]*_.*//')
        if [[ ! -v valid_cells[$cid] ]]; then
            n_skipped=$(( n_skipped + 1 ))
            continue
        fi
        r2="${r1/_R1_/_R2_}"
        if [[ ! -f "$r2" ]]; then
            r2=$(find "$rundir" -maxdepth 1 -name "${cid}_*_R2_*.fastq.gz" | sort | head -1)
        fi
        [[ -f "$r2" ]] || { echo "[${chip}] WARNING: no R2 for ${cid}, skipping"; continue; }
        cell_ids+=("$cid")
        r1_paths+=("$r1")
        r2_paths+=("$r2")
    done < <(find "$rundir" -maxdepth 1 -name "*_R1_*.fastq.gz" | sort)

    local n_cells=${#cell_ids[@]}
    echo "[${chip}] ${n_cells} / ${n_fastq} cells selected for alignment (${n_skipped} non-CSV cells skipped)"
    echo "[${chip}] Aligning (rolling pool: ${CELL_PARALLEL} parallel, ${CELL_THREADS} threads each)..."

    # Rolling pool: as soon as one cell finishes, the next starts immediately.
    # This eliminates the idle time of fixed batches where all slots wait for
    # the slowest cell before the batch advances.
    local pids=() n_done=0
    for (( i=0; i<n_cells; i++ )); do
        # If at capacity, wait for any one job to finish then prune dead PIDs.
        if [[ ${#pids[@]} -ge $CELL_PARALLEL ]]; then
            wait -n 2>/dev/null || true
            local live=()
            for _p in "${pids[@]}"; do
                kill -0 "$_p" 2>/dev/null && live+=("$_p") || (( n_done++ )) || true
            done
            pids=("${live[@]}")
        fi
        align_cell "${cell_ids[$i]}" "$chip" "$celldir" "${r1_paths[$i]}" "${r2_paths[$i]}" &
        pids+=($!)
        # Progress report every 200 cells submitted.
        if (( (i+1) % 200 == 0 )); then
            echo "[${chip}] $(( i+1 )) / ${n_cells} cells submitted (${n_done} confirmed done)..."
        fi
    done
    # Drain remaining jobs.
    for _p in "${pids[@]}"; do wait "$_p" || true; done

    # ── Merge per sample ──────────────────────────────────────────────────────
    # For each region, collect the per-cell BAMs listed in the CSV rglist file,
    # then merge directly into the per-sample output BAM.
    local any_merged=0
    for region in "${regions[@]}"; do
        local out_bam="${chipdir}/${region}.bam"
        if [[ -f "$out_bam" ]]; then
            echo "[${chip}/${region}] already exists — skipping"
            continue
        fi

        local region_bams=()
        while IFS= read -r cell_id; do
            local bam="${celldir}/${cell_id}.bam"
            [[ -f "$bam" ]] && region_bams+=("$bam")
        done < "${RGLIST_DIR}/${chip}_${region}.txt"

        if [[ ${#region_bams[@]} -eq 0 ]]; then
            echo "[${chip}/${region}] WARNING: no aligned BAMs found for this region"
            continue
        fi

        echo "[${chip}/${region}] merging ${#region_bams[@]} cells → ${out_bam}"
        merge_bams "$out_bam" "${region_bams[@]}"
        samtools index "$out_bam"
        echo "[${chip}/${region}] done → ${out_bam}"
        any_merged=1
    done

    [[ $any_merged -eq 1 ]] && rm -rf "$celldir"
    echo "[${chip}] complete  (run run_bqsr_only.sh for BQSR)"
}

# ── Dispatch ──────────────────────────────────────────────────────────────────
echo "Chips to process: ${!CHIP_TO_RUN[*]}"
echo "Parallelism:      ${CHIP_PARALLEL} chip(s) × ${CELL_PARALLEL} cells × ${CELL_THREADS} threads"

_pids=()
for chip in "${!CHIP_TO_RUN[@]}"; do
    process_chip "$chip" &
    _pids+=($!)
    if [[ ${#_pids[@]} -ge $CHIP_PARALLEL ]]; then
        for _pid in "${_pids[@]}"; do
            wait "$_pid" || echo "WARNING: job $_pid exited with error"
        done
        _pids=()
    fi
done

for _pid in "${_pids[@]}"; do
    wait "$_pid" || echo "WARNING: job $_pid exited with error"
done

echo "All chips complete."
