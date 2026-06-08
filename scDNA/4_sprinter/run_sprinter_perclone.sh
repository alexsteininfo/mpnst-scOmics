#!/usr/bin/env bash
# Run SPRINTER S-phase estimation for each predefined clone and per-region
# healthy cell population.
#
# Clones and healthy populations are defined in the metadata file produced by
# scDNA/2_clustering/create_cell_metadata.R.
#
# For each group (18 tumour clones + 6 healthy-per-region = 24 groups), three
# SPRINTER runs are performed:
#   combined  — 10X + DLP+ cells together
#   10X       — 10X cells only
#   DLP       — DLP+ cells only
#
# Groups whose CHISEL RDR file is missing are silently skipped (with a log
# entry). This means 10X-only and combined runs are skipped until
# chisel_rdr_10x.sh has been run.
#
# Input:   scDNA/2_clustering/metadata_cells.tsv   (run create_cell_metadata.R first)
# Output:  $BASE_OUTDIR/<group>/<mode>/sprinter.output.tsv.gz
#
# Requires: micromamba environment "sprinter"

set -euo pipefail

# ── Paths ──────────────────────────────────────────────────────────────────────

METADATA="/srv/home/aste0033/GitHub/mpnst-scOmics/scDNA/2_clustering/metadata_cells.tsv"
BASE_OUTDIR="/srv/home/aste0033/projects/MPNST/scDNA/SPRINTER/perclone"
RESOURCE_DIR="/srv/home/aste0033/projects/MPNST/scDNA/SPRINTER/resources_hg38"
LOGFILE="${BASE_OUTDIR}/run_log.txt"

REF="/mnt/iribhm/genomes/hg38-ngs/hg38.fa"
THREADS=24

# ── SPRINTER parameters ────────────────────────────────────────────────────────

MIN_READS=100000   # minimum mapped reads per cell
RT_READS=200       # target reads/bin for replication-timing analysis
CN_READS=1000      # target reads/bin for copy-number analysis
MAX_PLOIDY=8       # high ploidy for MPNST (highly aneuploid)

# ── Environment ────────────────────────────────────────────────────────────────

eval "$(micromamba shell hook --shell bash)"
micromamba activate sprinter

log() { echo "[$(date '+%H:%M:%S')] $*"; }
die() { echo "[$(date '+%H:%M:%S')] ERROR: $*" >&2; exit 1; }

[[ -f "$METADATA"      ]] || die "Metadata not found: $METADATA — run create_cell_metadata.R first"
[[ -f "$REF"           ]] || die "Reference FASTA not found: $REF"

RTSCORES_HG38="${RESOURCE_DIR}/rtscores_hg38.csv.gz"
GAPS_HG38="${RESOURCE_DIR}/gaps_hg38.tsv"

[[ -f "$RTSCORES_HG38" ]] || die "hg38 RT scores not found — run liftover_dlpp.sh first"
[[ -f "$GAPS_HG38"     ]] || die "hg38 gaps file not found — run liftover_dlpp.sh first"

mkdir -p "$BASE_OUTDIR"

# ── Helper: run one SPRINTER job ───────────────────────────────────────────────
#
# Arguments:
#   $1  group_name   output directory name, e.g. "healthy_R1" or "clone_P_2"
#   $2  clone_label  value in the clone_label metadata column
#   $3  region       region filter for healthy groups; empty string for clones
#   $4  mode         "combined" | "10X" | "DLP"

run_sprinter() {
  local group_name="$1"
  local clone_label="$2"
  local region="$3"
  local mode="$4"

  local outdir="${BASE_OUTDIR}/${group_name}/${mode}"
  local out="${outdir}/sprinter.output.tsv.gz"

  if [[ -f "$out" ]]; then
    log "${group_name}/${mode}: output already exists — skipping"
    return 0
  fi

  local tmpdir
  tmpdir=$(mktemp -d)

  # Build metadata filter and write "barcode <TAB> chisel_rdr_path" for this group
  # Metadata columns (1-indexed, tab-separated, with header on row 1):
  #   1=barcode  2=region  3=technology  4=kmeans_cluster
  #   5=clone_label  6=cell_type  7=chisel_rdr_path
  local awk_prog
  if [[ "$clone_label" == "healthy" ]]; then
    # Healthy group: filter by clone_label=healthy AND region AND (optionally) technology
    if [[ "$mode" == "combined" ]]; then
      awk_prog='NR>1 && $5==cl && $2==rg'
    else
      awk_prog='NR>1 && $5==cl && $2==rg && $3==md'
    fi
    awk -F'\t' -v cl="$clone_label" -v rg="$region" -v md="$mode" \
      "$awk_prog { print \$1\"\t\"\$7 }" "$METADATA" > "$tmpdir/bc_path.txt"
  else
    # Tumour clone: filter by clone_label AND (optionally) technology
    if [[ "$mode" == "combined" ]]; then
      awk_prog='NR>1 && $5==cl'
    else
      awk_prog='NR>1 && $5==cl && $3==md'
    fi
    awk -F'\t' -v cl="$clone_label" -v md="$mode" \
      "$awk_prog { print \$1\"\t\"\$7 }" "$METADATA" > "$tmpdir/bc_path.txt"
  fi

  local n_cells
  n_cells=$(wc -l < "$tmpdir/bc_path.txt")

  if [[ "$n_cells" -eq 0 ]]; then
    log "${group_name}/${mode}: 0 cells — skipping"
    printf '%s | %s/%s | skipped (0 cells in metadata)\n' \
      "$(date '+%Y-%m-%d %H:%M:%S')" "$group_name" "$mode" >> "$LOGFILE"
    rm -rf "$tmpdir"
    return 0
  fi

  # Collect unique CHISEL RDR paths for this subset
  cut -f2 "$tmpdir/bc_path.txt" | sort -u > "$tmpdir/rdr_paths.txt"

  # Check all RDR files exist (missing = 10X CHISEL not yet run)
  local missing=0
  while IFS= read -r rdr_path; do
    if [[ ! -f "$rdr_path" ]]; then
      log "${group_name}/${mode}: CHISEL RDR not found: ${rdr_path} — skipping"
      printf '%s | %s/%s | skipped (missing RDR: %s)\n' \
        "$(date '+%Y-%m-%d %H:%M:%S')" "$group_name" "$mode" "$rdr_path" >> "$LOGFILE"
      missing=1
    fi
  done < "$tmpdir/rdr_paths.txt"

  if [[ "$missing" -eq 1 ]]; then
    rm -rf "$tmpdir"
    return 0
  fi

  # For each RDR file: filter rows to the target barcodes for that file
  > "$tmpdir/merged_rdr.tsv"
  while IFS= read -r rdr_path; do
    # Barcodes that map to this specific RDR file
    awk -F'\t' -v p="$rdr_path" '$2==p { print $1 }' "$tmpdir/bc_path.txt" \
      > "$tmpdir/bc_this_rdr.txt"
    # Keep only RDR rows whose CELL column (col 4) is in the target set
    awk 'NR==FNR { b[$1]=1; next } $4 in b' \
      "$tmpdir/bc_this_rdr.txt" "$rdr_path" >> "$tmpdir/merged_rdr.tsv"
  done < "$tmpdir/rdr_paths.txt"

  local n_rows
  n_rows=$(wc -l < "$tmpdir/merged_rdr.tsv")

  if [[ "$n_rows" -eq 0 ]]; then
    log "${group_name}/${mode}: no RDR rows after barcode filtering — skipping"
    printf '%s | %s/%s | skipped (empty RDR after filter)\n' \
      "$(date '+%Y-%m-%d %H:%M:%S')" "$group_name" "$mode" >> "$LOGFILE"
    rm -rf "$tmpdir"
    return 0
  fi

  # Build fixclones TSV: assigns every cell to clone "1"
  cut -f1 "$tmpdir/bc_path.txt" | sort -u | \
    awk 'BEGIN{print "CELL\tCLONE"} {print $0"\t1"}' > "$tmpdir/fixclones.tsv"

  # Gzip the merged RDR (SPRINTER expects compressed input)
  gzip -c "$tmpdir/merged_rdr.tsv" > "$tmpdir/input.tsv.gz"

  mkdir -p "$outdir"
  log "${group_name}/${mode}: starting SPRINTER (${n_cells} cells, ${n_rows} RDR rows)"

  local exit_code=0
  (
    cd "$outdir"
    sprinter "$tmpdir/input.tsv.gz" \
      --refgenome   "$REF" \
      --minreads    "$MIN_READS" \
      --rtreads     "$RT_READS" \
      --cnreads     "$CN_READS" \
      --maxploidy   "$MAX_PLOIDY" \
      -j            "$THREADS" \
      -o            "$out" \
      dev \
      --fixclones   "$tmpdir/fixclones.tsv" \
      --rtscores    "$RTSCORES_HG38" \
      --gapsfile    "$GAPS_HG38"
  ) || exit_code=$?

  if [[ "$exit_code" -eq 0 ]] && [[ -f "$out" ]]; then
    log "${group_name}/${mode}: done → ${out}"
    printf '%s | %s/%s | done\n' \
      "$(date '+%Y-%m-%d %H:%M:%S')" "$group_name" "$mode" >> "$LOGFILE"
  else
    log "${group_name}/${mode}: SPRINTER failed (exit ${exit_code})"
    printf '%s | %s/%s | failed (exit %s)\n' \
      "$(date '+%Y-%m-%d %H:%M:%S')" "$group_name" "$mode" "$exit_code" >> "$LOGFILE"
  fi

  rm -rf "$tmpdir"
}

# ── Collect groups from metadata ───────────────────────────────────────────────

# Healthy: unique regions where cell_type == "healthy"
mapfile -t HEALTHY_REGIONS < <(
  awk -F'\t' 'NR>1 && $6=="healthy" { print $2 }' "$METADATA" | sort -u
)

# Tumour clones: unique clone_labels where cell_type == "tumour"
mapfile -t CLONE_LABELS < <(
  awk -F'\t' 'NR>1 && $6=="tumour" { print $5 }' "$METADATA" | sort -u
)

log "Healthy regions: ${HEALTHY_REGIONS[*]}"
log "Tumour clones:   ${CLONE_LABELS[*]}"
log "Modes: combined 10X DLP"
log "Output: ${BASE_OUTDIR}"
log "Log: ${LOGFILE}"
echo ""

# ── Run all groups ─────────────────────────────────────────────────────────────

for region in "${HEALTHY_REGIONS[@]}"; do
  for mode in combined 10X DLP; do
    run_sprinter "healthy_${region}" "healthy" "$region" "$mode"
  done
done

for clone_label in "${CLONE_LABELS[@]}"; do
  group_name="clone_${clone_label}"
  for mode in combined 10X DLP; do
    run_sprinter "$group_name" "$clone_label" "" "$mode"
  done
done

log "=== All per-clone SPRINTER runs complete ==="
log "Outputs: ${BASE_OUTDIR}/<group>/<mode>/sprinter.output.tsv.gz"
log "Log:     ${LOGFILE}"
