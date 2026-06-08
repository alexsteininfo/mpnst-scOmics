# CHISEL — RDR computation for DLP+ data

CHISEL (`chisel_rdr`) counts reads per 50 kb genomic bin per cell from a multi-cell BAM file and computes the read-depth ratio (RDR) against a simulated normal sample. The resulting `rdr.tsv` is used as input for SPRINTER.

## Scripts

| Script | Purpose |
|--------|---------|
| `chisel_rdr_dlpp.sh` | Runs `chisel_rdr` on all six DLP+ per-sample BAMs (R1–R5, P) and renames output columns for SPRINTER compatibility |
| `chisel_rdr_10x.sh` | Equivalent script for 10x scDNA BAMs (pending 10x data) |

## Pipeline steps (`chisel_rdr_dlpp.sh`)

1. For each sample, locate the per-sample BAM at `bwa_output/<chip>/<sample>.bam` (produced by `fastq_to_bam_dlpp.sh`).
2. Run `chisel_rdr` with 50 kb bins, `--cellprefix RG:Z:` (DLP+ cells are identified by read-group tag), and a minimum of 100,000 reads per cell.
3. The tool simulates a diploid normal BAM via ART Illumina for mappability correction, then counts reads per bin per cell in the tumour BAM.
4. Rename the output columns for SPRINTER: `CHROMOSOME→CHR`, `NORMAL→NORM_COUNT`, `RDR→RAW_RDR`.

**Output:** `$CHISEL_DIR/<sample>/rdr.tsv` — tab-separated file with columns `CHR START END CELL NORM_COUNT COUNT RAW_RDR`.

## Patches applied to the `chisel` micromamba environment

Two bugs in the installed CHISEL package affect DLP+ data and were patched in-place. The patches are minimal and do not change algorithmic behaviour.

### 1. ART read-length cap — `chisel/bin/chisel_rdr.py` lines 251 and 290

**Problem:** `chisel_rdr` auto-detects the read length from the BAM (151 bp for these libraries) and passes it directly to ART Illumina. The `HS25` sequencer profile hard-caps at 150 bp, so ART exits with an error and no output is produced.

**Fix:** Cap the read length passed to ART at 150: `art_read = min(read, 150)`. A 1 bp difference is negligible for read-depth normalisation.

```python
# Before
cmd = '... -l {} ...'.format(..., read, ...)

# After
art_read = min(read, 150)  # HS25 profile supports max 150 bp
cmd = '... -l {} ...'.format(..., art_read, ...)
```

Applied to both the paired-end (`sequencing_paired`) and single-end (`sequencing_single`) code paths.

### 2. Cell barcode regex — `chisel/RDREstimator.py` line 213

**Problem:** CHISEL extracts per-read cell identifiers using the awk regex `[ACGT]+` after the tag prefix. DLP+ cell IDs take the form `LES4677A1_124574_44x_26y` (digits, underscores, mixed case) — none of which match `[ACGT]+`. As a result, zero cells were found despite ~688k reads per cell on average.

**Fix:** Broaden the character class to `[A-Za-z0-9_]+`.

```python
# Before
"awk '... /{}[ACGT]+{}/ ...'"

# After
"awk '... /{}[A-Za-z0-9_]+{}/ ...'"
```

This change is backward-compatible: standard 10x barcodes (ACGT-only) still match the wider class.
