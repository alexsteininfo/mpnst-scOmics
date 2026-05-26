#!/usr/bin/env bash
set -euo pipefail

source /srv/home/aste0033/CellRanger-DNA/setup_env.sh

FASTA="/mnt/iribhm/genomes/hg38-ngs/hg38.fa"
CONTIG_DEFS="/srv/home/aste0033/CellRanger-DNA/hg38_contig_defs.json"
REF_DIR="/srv/home/aste0033/projects/MPNST/scDNA/10X/reference"

mkdir -p "$REF_DIR"
cd "$REF_DIR"

cellranger-dna mkref "$FASTA" "$CONTIG_DEFS"

echo "Reference built at: ${REF_DIR}/refdata-hg38"
