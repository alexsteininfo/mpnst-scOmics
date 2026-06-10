#!/usr/bin/env bash
# Lift SPRINTER's bundled hg19 resources (rtscores, gap track) to hg38.
#
# SPRINTER ships with hg19-coordinate replication timing scores and gap
# regions. Because the DLP+ libraries are aligned to hg38, these files
# must be lifted once before running any SPRINTER analysis scripts.
#
# Approach:
#   For each hg19 50 kb bin, the midpoint is lifted to hg38 via pyliftover
#   (pure-Python reimplementation of UCSC liftOver; no native binary needed).
#   The hg38 position is snapped to its containing 50 kb bin. Bins whose
#   midpoints do not lift are dropped. Where multiple hg19 bins fall into
#   the same hg38 bin, scores are averaged.
#
# Does NOT require micromamba / conda activation — uses the sprinter env
# Python binary directly via its full path.
#
# Output (written to $RESOURCE_DIR):
#   rtscores_hg38.csv.gz       — replication timing scores on hg38 50 kb bins
#   gaps_hg38.tsv              — mapping gap regions in hg38 coordinates
#   hg19ToHg38.over.chain.gz   — UCSC chain file (downloaded if absent)

set -euo pipefail

# ── Paths ──────────────────────────────────────────────────────────────────────

RESOURCE_DIR="/srv/home/aste0033/projects/MPNST/scDNA/SPRINTER/resources_hg38"
CHAIN_URL="https://hgdownload.cse.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz"
CHAIN_FILE="${RESOURCE_DIR}/hg19ToHg38.over.chain.gz"

RTSCORES_HG38="${RESOURCE_DIR}/rtscores_hg38.csv.gz"
GAPS_HG38="${RESOURCE_DIR}/gaps_hg38.tsv"

# Python binary from the sprinter micromamba environment (no activation needed)
PYTHON="/srv/home/aste0033/micromamba/envs/sprinter/bin/python3"

# ── Checks ─────────────────────────────────────────────────────────────────────

log() { echo "[$(date '+%H:%M:%S')] $*"; }
die() { echo "[$(date '+%H:%M:%S')] ERROR: $*" >&2; exit 1; }

[[ -x "$PYTHON" ]] || die "Python not found: $PYTHON"
"$PYTHON" -c "import pyliftover" 2>/dev/null || die "pyliftover not installed — run: $PYTHON -m pip install pyliftover"

mkdir -p "$RESOURCE_DIR"

# ── Download chain file ────────────────────────────────────────────────────────

if [[ ! -f "$CHAIN_FILE" ]]; then
    log "Downloading hg19ToHg38 chain file..."
    curl -fsSL "$CHAIN_URL" -o "$CHAIN_FILE"
    log "Chain file saved to ${CHAIN_FILE}"
else
    log "Chain file already present: ${CHAIN_FILE}"
fi

# ── Liftover ──────────────────────────────────────────────────────────────────

SPRINTER_RESOURCES="$("$PYTHON" -c \
    'import sprinter, os; print(os.path.join(os.path.dirname(sprinter.__file__), "resources"))')"

log "Lifting hg19 SPRINTER resources → hg38..."

"$PYTHON" - \
    "$SPRINTER_RESOURCES" \
    "$CHAIN_FILE" \
    "$RTSCORES_HG38" \
    "$GAPS_HG38" <<'PYEOF'
import sys, gzip, csv, os
from collections import defaultdict
from pyliftover import LiftOver

resources_dir, chain_file, out_rtscores, out_gaps = sys.argv[1:]

lo = LiftOver(chain_file)

# ── Liftover rtscores ─────────────────────────────────────────────────────────
print("  Lifting over rtscores.csv.gz hg19 → hg38...", flush=True)

rtscores_hg19 = os.path.join(resources_dir, 'rtscores.csv.gz')
with gzip.open(rtscores_hg19, 'rt') as f:
    reader = csv.reader(f)
    header = next(reader)   # ['chr', 'start', 'stop', 'A549', ...]
    rows   = list(reader)

# Lift midpoint of each hg19 50 kb bin to hg38, snap to 50 kb bin
hg38_bins = defaultdict(list)
for i, row in enumerate(rows):
    chrom = row[0] if row[0].startswith('chr') else 'chr' + row[0]
    start = int(float(row[1]))
    stop  = int(float(row[2]))
    mid   = (start + stop) // 2
    result = lo.convert_coordinate(chrom, mid)
    if not result:
        continue
    new_chrom, new_pos = result[0][0], result[0][1]
    hg38_bins[(new_chrom, (new_pos // 50000) * 50000)].append(i)

AUTOSOMES = {f'chr{i}' for i in range(1, 23)}

def chrom_sort_key(k):
    return int(k[0].replace('chr', ''))

n_score_cols = len(header) - 3
with gzip.open(out_rtscores, 'wt', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(header)
    for (chrom, bin_start) in sorted(
            ((c, s) for c, s in hg38_bins if c in AUTOSOMES), key=chrom_sort_key):
        row_idxs = hg38_bins[(chrom, bin_start)]
        scores   = []
        for col in range(3, 3 + n_score_cols):
            vals = []
            for i in row_idxs:
                try:
                    vals.append(float(rows[i][col]))
                except (ValueError, IndexError):
                    pass
            scores.append(str(sum(vals) / len(vals)) if vals else 'NA')
        writer.writerow([chrom, bin_start, bin_start + 50000] + scores)

n_lifted = sum(len(v) for v in hg38_bins.values())
n_unmap  = len(rows) - n_lifted
print(f"  rtscores: {len(rows)} hg19 bins → {len(hg38_bins)} hg38 bins "
      f"({n_unmap} unmapped/dropped)", flush=True)

# ── Build hg38 gaps file from UCSC centromere + gap tracks ───────────────────
# hg19 centromere/telomere coordinates do not exist in the chain file, so
# liftover produces almost nothing. The correct approach is to fetch the
# hg38 centromere and gap annotations directly from UCSC.
import urllib.request

print("  Building gaps_hg38.tsv from UCSC hg38 tracks...", flush=True)

AUTOSOMES = {str(i) for i in range(1, 23)}

rows_out = []

# 1. Centromeres — from the dedicated hg38 centromere track
cen_url = "https://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/centromeres.txt.gz"
with urllib.request.urlopen(cen_url) as resp:
    import io
    with gzip.open(io.BytesIO(resp.read()), 'rt') as f:
        for line in f:
            parts = line.rstrip('\n').split('\t')
            # columns: bin, chrom, chromStart, chromEnd, name
            chrom = parts[1].replace('chr', '')
            if chrom not in AUTOSOMES:
                continue
            rows_out.append((int(chrom), int(parts[2]), int(parts[3]), 'centromere'))

# 2. Telomeres and other structural gaps — from the gap track
gap_url = "https://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/gap.txt.gz"
keep_types = {'telomere', 'heterochromatin', 'short_arm'}
with urllib.request.urlopen(gap_url) as resp:
    with gzip.open(io.BytesIO(resp.read()), 'rt') as f:
        for line in f:
            parts = line.rstrip('\n').split('\t')
            # columns: bin, chrom, chromStart, chromEnd, ix, n, size, type, bridge
            chrom    = parts[1].replace('chr', '')
            gap_type = parts[7]
            if chrom not in AUTOSOMES or gap_type not in keep_types:
                continue
            rows_out.append((int(chrom), int(parts[2]), int(parts[3]), gap_type))

rows_out.sort()
with open(out_gaps, 'w') as f:
    f.write('CHR\tSTART_POS\tEND_POS\tTYPE\n')
    for chrom, start, end, gap_type in rows_out:
        f.write(f'{chrom}\t{start}\t{end}\t{gap_type}\n')

print(f"  gaps: {len(rows_out)} hg38 regions written (centromeres + telomeres)", flush=True)
print("  Done.", flush=True)
PYEOF

log "hg38 resources written to ${RESOURCE_DIR}:"
log "  ${RTSCORES_HG38}"
log "  ${GAPS_HG38}"
