#!/usr/bin/env bash
# Run SPRINTER S-phase estimation per region and technology.
#
# Each of the 12 runs (6 regions × 2 technologies) pools all cells from that
# region/technology — healthy, tumour clones, and unassigned — with their real
# clone assignments passed via --fixclones. Cluster-2 "removed" cells are
# excluded. This allows SPRINTER to perform within-run GC/RT bias correction
# across all clones simultaneously.
#
# Input:   scDNA/2_clustering/metadata_cells.tsv   (run create_cell_metadata.R first)
# Output:  $BASE_OUTDIR/<region>/<technology>/sprinter.output.tsv.gz
#
# Requires: micromamba environment "sprinter"

set -euo pipefail

# ── Paths ──────────────────────────────────────────────────────────────────────

METADATA="/srv/home/aste0033/GitHub/mpnst-scOmics/scDNA/2_clustering/metadata_cells.tsv"
BASE_OUTDIR="/srv/home/aste0033/projects/MPNST/scDNA/SPRINTER/persample"
RESOURCE_DIR="/srv/home/aste0033/projects/MPNST/scDNA/SPRINTER/resources_hg38"
LOGFILE="${BASE_OUTDIR}/run_log.txt"

REF="/mnt/iribhm/genomes/hg38-ngs/hg38.fa"
THREADS=24

# ── SPRINTER parameters ────────────────────────────────────────────────────────

MIN_READS=100000
RT_READS=200
CN_READS=1000
MAX_PLOIDY=8

# ── Environment ────────────────────────────────────────────────────────────────

eval "$(micromamba shell hook --shell bash)"
micromamba activate sprinter

log() { echo "[$(date '+%H:%M:%S')] $*"; }
die() { echo "[$(date '+%H:%M:%S')] ERROR: $*" >&2; exit 1; }

[[ -f "$METADATA"      ]] || die "Metadata not found: $METADATA — run create_cell_metadata.R first"
[[ -f "$REF"           ]] || die "Reference FASTA not found: $REF"

RTSCORES_HG38="${RESOURCE_DIR}/rtscores_hg38.csv.gz"
GAPS_HG38="${RESOURCE_DIR}/gaps_hg38.tsv"

[[ -f "$RTSCORES_HG38" ]] || die "hg38 RT scores not found — run resources/liftover_dlpp.sh first"
[[ -f "$GAPS_HG38"     ]] || die "hg38 gaps file not found — run resources/liftover_dlpp.sh first"

mkdir -p "$BASE_OUTDIR"

# ── Helper: run one SPRINTER job ───────────────────────────────────────────────
#
# Arguments:
#   $1  region       one of: P R1 R2 R3 R4 R5
#   $2  technology   one of: 10X DLP

run_sprinter_persample() {
  local region="$1"
  local tech="$2"

  local outdir="${BASE_OUTDIR}/${region}/${tech}"
  local out="${outdir}/sprinter.output.tsv.gz"

  if [[ -f "$out" ]]; then
    log "${region}/${tech}: output already exists — skipping"
    return 0
  fi

  local tmpdir
  tmpdir=$(mktemp -d)

  # Filter metadata: region == X AND technology == Y AND cell_type != "removed"
  # Metadata columns (1-indexed, tab-separated, header on row 1):
  #   1=barcode  2=region  3=technology  4=kmeans_cluster
  #   5=clone_label  6=cell_type  7=chisel_rdr_path
  awk -F'\t' -v rg="$region" -v tc="$tech" \
    'NR>1 && $2==rg && $3==tc && $6!="removed" { print $1"\t"$5"\t"$7 }' \
    "$METADATA" > "$tmpdir/bc_clone_path.txt"

  local n_cells
  n_cells=$(wc -l < "$tmpdir/bc_clone_path.txt")

  if [[ "$n_cells" -eq 0 ]]; then
    log "${region}/${tech}: 0 cells — skipping"
    printf '%s | %s/%s | skipped (0 cells in metadata)\n' \
      "$(date '+%Y-%m-%d %H:%M:%S')" "$region" "$tech" >> "$LOGFILE"
    rm -rf "$tmpdir"
    return 0
  fi

  # Collect unique CHISEL RDR paths for this subset
  cut -f3 "$tmpdir/bc_clone_path.txt" | sort -u > "$tmpdir/rdr_paths.txt"

  local n_rdr_files
  n_rdr_files=$(wc -l < "$tmpdir/rdr_paths.txt")
  if [[ "$n_rdr_files" -gt 1 ]]; then
    log "${region}/${tech}: WARNING — ${n_rdr_files} RDR paths found (expected 1); merging"
  fi

  # Check all RDR files exist
  local missing=0
  while IFS= read -r rdr_path; do
    if [[ ! -f "$rdr_path" ]]; then
      log "${region}/${tech}: CHISEL RDR not found: ${rdr_path} — skipping"
      printf '%s | %s/%s | skipped (missing RDR: %s)\n' \
        "$(date '+%Y-%m-%d %H:%M:%S')" "$region" "$tech" "$rdr_path" >> "$LOGFILE"
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
    awk -F'\t' -v p="$rdr_path" '$3==p { print $1 }' "$tmpdir/bc_clone_path.txt" \
      > "$tmpdir/bc_this_rdr.txt"
    awk 'NR==FNR { b[$1]=1; next } $4 in b' \
      "$tmpdir/bc_this_rdr.txt" "$rdr_path" >> "$tmpdir/merged_rdr.tsv"
  done < "$tmpdir/rdr_paths.txt"

  local n_rows
  n_rows=$(wc -l < "$tmpdir/merged_rdr.tsv")

  if [[ "$n_rows" -eq 0 ]]; then
    log "${region}/${tech}: no RDR rows after barcode filtering — skipping"
    printf '%s | %s/%s | skipped (empty RDR after filter)\n' \
      "$(date '+%Y-%m-%d %H:%M:%S')" "$region" "$tech" >> "$LOGFILE"
    rm -rf "$tmpdir"
    return 0
  fi

  # Build fixclones TSV with integer clone IDs.
  # SPRINTER assumes numeric CLONE labels internally (max()+1 arithmetic when no
  # normal clone is detected; np.isnan checks in plot code). String labels trigger
  # crashes in both code paths. Integer IDs are safe; downstream R scripts recover
  # original clone labels from metadata via barcode join, so the integers are
  # transparent to the analysis.
  awk '{print $2}' "$tmpdir/bc_clone_path.txt" | sort -u \
    | awk '{print $1"\t"NR}' > "$tmpdir/clone_map.txt"
  { printf 'CELL\tCLONE\n'
    awk 'NR==FNR { id[$1]=$2; next } { print $1"\t"id[$2] }' \
      "$tmpdir/clone_map.txt" "$tmpdir/bc_clone_path.txt"; } \
    > "$tmpdir/fixclones.tsv"

  # Gzip the merged RDR (SPRINTER expects compressed input)
  gzip -c "$tmpdir/merged_rdr.tsv" > "$tmpdir/input.tsv.gz"

  mkdir -p "$outdir"
  log "${region}/${tech}: starting SPRINTER (${n_cells} cells, ${n_rows} RDR rows)"

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
    log "${region}/${tech}: done → ${out}"
    printf '%s | %s/%s | done\n' \
      "$(date '+%Y-%m-%d %H:%M:%S')" "$region" "$tech" >> "$LOGFILE"
  else
    log "${region}/${tech}: SPRINTER failed (exit ${exit_code})"
    printf '%s | %s/%s | failed (exit %s)\n' \
      "$(date '+%Y-%m-%d %H:%M:%S')" "$region" "$tech" "$exit_code" >> "$LOGFILE"
  fi

  rm -rf "$tmpdir"
}

# ── Run all 12 combinations ────────────────────────────────────────────────────

REGIONS=(P R1 R2 R3 R4 R5)
TECHS=(10X DLP)

log "Regions:    ${REGIONS[*]}"
log "Technologies: ${TECHS[*]}"
log "Output:     ${BASE_OUTDIR}"
log "Log:        ${LOGFILE}"
echo ""

for region in "${REGIONS[@]}"; do
  for tech in "${TECHS[@]}"; do
    run_sprinter_persample "$region" "$tech"
  done
done

log "=== All per-sample SPRINTER runs complete ==="
log "Outputs: ${BASE_OUTDIR}/<region>/<technology>/sprinter.output.tsv.gz"
log "Log:     ${LOGFILE}"
