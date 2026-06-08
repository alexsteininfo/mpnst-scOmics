# 10x scDNA-seq Processing: fastq → BAM

## Project

Patient **FIT208**, MPNST. Six spatial tumor regions sequenced with 10x Genomics
single-cell DNA (scDNA-seq). Goal: align per-lane Illumina fastqs to produce
per-cell BAM files.

**Write permissions:**
- Output: `/srv/home/aste0033/projects/MPNST/scDNA/10X/`
- Source fastqs: `/mnt/iribhm/people/mtarabichi/MPNST/pvanloo_20250618/Haixi/10x_DNA/` (read-only)
- Reference genomes: `/mnt/iribhm/genomes/` (read-only)

---

## Library structure and barcode encoding

### Read layout

10x Chromium Single Cell CNV (scDNA) uses a 3-read Illumina library:

```
R1  [16 bp barcode][8 bp UMI][genomic insert, if read length > 24 bp]
R2  [genomic insert — the sequence that gets aligned]
I1  [8 bp sample index — used by Illumina for demultiplexing only; discarded afterwards]
```

- **R1** is not aligned. It is parsed by `extract_barcodes.py` to extract the cell
  barcode (R1[0:16]) and UMI (R1[16:24]). The remaining bases (if present) are ignored.
- **R2** carries the genomic DNA sequence and is the only read passed to the aligner.
- Libraries are typically sequenced at 151 bp, giving R1 ≫ 24 bp. The extra bases
  on R1 beyond position 24 are adapter/primer sequence and are never used.

### Cell barcodes and the whitelist

Each 10x gel bead carries one of 737,280 known 16-mer DNA sequences printed on it —
the **barcode whitelist**. Every cell captured in a droplet inherits one whitelist
barcode from its bead. Any 16-mer in a sequenced R1 that matches a whitelist entry
(exactly or within 1 mismatch) is accepted as a real cell barcode.

**scDNA-specific whitelist:** `737K-crdna-v1.txt` (737,280 barcodes, stored at
`/srv/home/aste0033/projects/MPNST/scDNA/10X/737K-crdna-v1.txt`).

> **Do NOT use `737K-august-2016.txt`** (GEX v2 gel beads). The two lists share
> only ~25% of barcodes — using the GEX whitelist discards ~75% of true scDNA cells.

**Hamming-1 correction:** `extract_barcodes.py` builds a correction map at runtime
(48 neighbours per whitelist barcode = 16 positions × 3 alternative bases; `N` is
also included to rescue first-cycle quality drops). For each read, the raw barcode is
looked up: exact match → accept as-is; single mismatch → correct to whitelist barcode;
ambiguous (matches >1 whitelist barcode) or >1 mismatch → discard read.

### CB and UR BAM tags

After correction, the barcode and UMI are embedded in the R2 FASTQ header as a
comment (`CB:Z:<barcode>\tUR:Z:<umi>`). BWA-MEM/BWA-MEME's `-C` flag propagates
FASTQ comments into the SAM output, converting them into proper BAM optional fields:

```
CB:Z:ACGTACGTACGTACGT   — corrected cell barcode (16 bp, matches whitelist)
UR:Z:TTCCAAGG           — raw UMI (8 bp, uncorrected)
```

Every aligned read in the BAM therefore carries the identity of the cell it came
from. Downstream tools (`samtools markdup --barcode-tag CB`, cell calling scripts,
alleleCounter) use `CB` to group reads by cell.

**Important:** in **RAW_MODE** (see below), CB tags contain the uncorrected barcode
read directly from R1. These tags are not guaranteed to be in the whitelist.

### RAW_MODE vs whitelist mode

| Mode | Reads kept | CB tag | Use case |
|---|---|---|---|
| **Whitelist** (`RAW_MODE=false`) | ~8–12% of reads | Whitelist-corrected 16-mer | Production runs; clean cell calling |
| **RAW** (`RAW_MODE=true`) | 100% of reads | Raw R1[0:16] barcode | Exploratory; post-hoc whitelist correction |

The current BAMs were produced with `RAW_MODE=true`. Whitelist correction was applied
post-hoc with `apply_whitelist_to_counts.py` (see Cell calling section).

> The expected read keep rate in whitelist mode for this data has not been directly
> measured; ~8–12% of unique barcode types survive whitelist filtering.

### Duplicate marking

`samtools markdup --barcode-tag CB` treats reads from different cells as independent:
only true within-cell PCR duplicates (same CB + same genomic position) are flagged.
Without `--barcode-tag CB`, reads from different cells landing at the same position
would be incorrectly marked as duplicates of each other.

---

## Data structure

### Sequencing runs (all read-only)

| Folder | Date | Samples | Size | Notes |
|---|---|---|---|---|
| `SC18143_1` | Feb 2019 | FIT208A1–A8 | 411 GB | Initial run; A3/A4 deep, A5–A8 minimal |
| `SC18143_2` | Jul 2019 | FIT208A5–A8 | 83 GB | Moderate-depth re-sequencing |
| `sc18143_2b` | Sep 2019 | FIT208A5–A8 | 215 GB | High-depth re-sequencing; used in analysis |

**FIT208A1 and FIT208A2** appear only in `SC18143_1` and are absent from the
analysis metadata (`10X_DNA_metadata.R`). Likely excluded due to QC failure; fastqs
retained but not processed.

### Sample → region mapping

| Sample | Region | Expected cells | Fastq source(s) |
|---|---|---|---|
| FIT208A3 | R1 (Primary) | ~2,220 | SC18143_1 only |
| FIT208A4 | R2 | ~1,440 | SC18143_1 only |
| FIT208A5 | R3 | ~106 | SC18143_1 + SC18143_2 + sc18143_2b |
| FIT208A6 | R4 | ~363 | SC18143_1 + SC18143_2 + sc18143_2b |
| FIT208A7 | R5 | ~1,513 | SC18143_1 + SC18143_2 + sc18143_2b |
| FIT208A8 | P  | ~1,154 | SC18143_1 + SC18143_2 + sc18143_2b |

A5–A8 are called "combo" samples in the original analysis code.

---

## Processing pipeline

**Active pipeline:** BWA-MEME + custom barcode extraction.
**Aligner:** BWA-MEME (drop-in BWA-MEM replacement using learned indices, ~1.5× faster).
Run via Apptainer SIF: `/mnt/iribhm/software/singularity/bwameme.sif`.

**MEME index:** only at `hg38-gatk/Homo_sapiens_assembly38.fa`; the `hg38-ngs/hg38.fa`
index lacks `.bwt.2bit.64` / `.0123` files needed by BWA-MEME mode 3.

### Pipeline steps

```
R1 + R2  →  extract_barcodes.py              barcode extraction + Hamming-1 correction
         →  apptainer bwa-meme mem -C -R      alignment + CB/UR tags + read group
         →  samtools fixmate -m               mate score tags (required by markdup)
         →  samtools sort -T <tmpdir>         coordinate sort
         →  samtools markdup --barcode-tag CB per-cell duplicate removal
         →  markdup.bam                       ← intermediate output (alignment done)

         →  run_bqsr_only.sh (separate step)
         →  gatk BaseRecalibrator
         →  gatk ApplyBQSR
         →  possorted_bam.bam                ← final BAM (BQSR-corrected base qualities)
```

BQSR is intentionally separated: `markdup.bam` is suitable for CNV analysis, and
BQSR is only needed if SNV / mitochondrial variant calling is planned downstream.

### Scripts

| Script | Purpose |
|---|---|
| `fastq_to_bam_bwa.sh` | Main pipeline: fastq → markdup.bam (all 6 samples) |
| `fastq_to_bam_bwa_test.sh` | Test run on FIT208A3 lane 4 only |
| `extract_barcodes.py` | Barcode extraction helper called by the alignment scripts |
| `run_bqsr_only.sh` | BQSR on existing markdup.bam → possorted_bam.bam |
| `run_pipeline.sh` | Orchestration: test → full run → validation |

### Prerequisites

```bash
module load SAMtools/1.18-GCC-12.3.0
# BWA-MEME runs via Apptainer (no module needed):
apptainer exec --cleanenv \
    --bind /mnt/iribhm:/mnt/iribhm \
    /mnt/iribhm/software/singularity/bwameme.sif \
    bwa-meme --version
# BQSR only (run_bqsr_only.sh):
module load GATK/4.2.4.1-GCCcore-10.2.0-Java-1.8
```

Known variant sites for BQSR (from `/mnt/iribhm/genomes/hg38-gatk/`):
- `Homo_sapiens_assembly38.dbsnp138.vcf.gz`
- `Mills_and_1000G_gold_standard.indels.hg38.vcf.gz`
- `Homo_sapiens_assembly38.known_indels.vcf.gz`

> Do not use files from `hg38-ana/`: their `.tbi` indices use non-`chr` contig names
> that are incompatible with the `chr`-prefixed reference.

### NFS performance note

BWA-MEME loads the genome index (~20 GB) from NFS on every cell startup. The script
caches the index to `/tmp` on first run to eliminate repeated NFS reads.

### Run

```bash
# Full pipeline (fastq → markdup.bam), all 6 samples:
bash code/0_preprocessing/10x/custom/fastq_to_bam_bwa.sh

# BQSR on completed markdup.bam (separate step, when needed):
bash code/0_preprocessing/10x/custom/run_bqsr_only.sh
```

### Verification checklist

**Before running:**
- [ ] `apptainer exec --cleanenv /mnt/iribhm/software/singularity/bwameme.sif bwa-meme 2>&1 | head -1`
- [ ] `samtools --version`
- [ ] `ls /srv/home/aste0033/projects/MPNST/scDNA/10X/737K-crdna-v1.txt`
- [ ] `ls /mnt/iribhm/genomes/hg38-gatk/Homo_sapiens_assembly38.fa.bwt.2bit.64`

**After each sample (e.g. FIT208A3):**
- [ ] `ls .../FIT208A3/` contains `markdup.bam`, `markdup.bam.bai`, `bwa.log`, `markdup.log`
- [ ] `samtools quickcheck .../FIT208A3/markdup.bam` exits 0
- [ ] `samtools view .../FIT208A3/markdup.bam | head -1 | grep -o "CB:Z:[A-Z]*"` — CB tag present

---

## Cell calling

Cell identity is not encoded in the BAM filename — thousands of cells are pooled per
BAM. "Cell calling" decides which barcodes represent real cells vs empty droplets or
background. Two complementary methods are used and compared.

### Methods

| Method | Script | Output |
|---|---|---|
| **Knee-point** (de novo) | `knee_point_cell_calling.R` | `valid_barcodes_knee[_wl].txt` |
| **Cell Ranger DNA** (original) | `extract_barcodes_from_rds.R` | `valid_barcodes_cellranger.txt` |

The Kneedle algorithm plots read count vs barcode rank on a log-log scale and finds
the inflection point separating real cells from background. Manual thresholds
(derived from the last Cell Ranger DNA-called cell per sample) override the automatic
detection and are hardcoded in `knee_point_cell_calling.R`.

### Whitelist correction modes

Both `knee_point_cell_calling.R` and `extract_barcodes_from_rds.R` support a suffix
variable to switch between raw and whitelist-corrected count files:

```r
# In knee_point_cell_calling.R — set at top of script or pass before source():
COUNT_SUFFIX <- ""      # reads barcode_counts.txt,    writes valid_barcodes_knee.txt
COUNT_SUFFIX <- "_wl"   # reads barcode_counts_wl.txt, writes valid_barcodes_knee_wl.txt

# In extract_barcodes_from_rds.R:
KNEE_SUFFIX <- ""       # compares against valid_barcodes_knee.txt
KNEE_SUFFIX <- "_wl"    # compares against valid_barcodes_knee_wl.txt
```

Override from the command line:
```bash
Rscript - <<'EOF'
COUNT_SUFFIX <- "_wl"
source("code/0_preprocessing/10x/custom/knee_point_cell_calling.R")
EOF
```

### Post-hoc whitelist correction

If BAMs were produced with `RAW_MODE=true`, the `barcode_counts.txt` files contain
raw (uncorrected) barcodes. Apply whitelist correction post-hoc without re-aligning:

```bash
python3 code/0_preprocessing/10x/custom/apply_whitelist_to_counts.py \
    --whitelist /srv/home/aste0033/projects/MPNST/scDNA/10X/737K-crdna-v1.txt \
    --outdir    /srv/home/aste0033/projects/MPNST/scDNA/10X/bwa_output
```

Produces `barcode_counts_wl.txt` per sample (exact + Hamming-1 correction, ambiguous
barcodes discarded). Run `knee_point_cell_calling.R` with `COUNT_SUFFIX <- "_wl"` afterwards.

### Results: Jaccard comparison (raw BAMs, May 2026)

Current BAMs were aligned with `RAW_MODE=true` (no whitelist filtering).
Cell Ranger DNA ground truth from `MPNST_scDNA_nopcf_raw_mtx.rds`.

**Raw barcodes** (no whitelist):

| Sample | Cell Ranger | Knee-point | Overlap | Jaccard |
|---|---|---|---|---|
| FIT208A3 (R1) | 1,421 | 5,256 | 1,421 | 0.270 |
| FIT208A4 (R2) | 996 | 3,842 | 996 | 0.259 |
| FIT208A5 (R3) | 82 | 304 | 82 | 0.270 |
| FIT208A6 (R4) | 265 | 984 | 265 | 0.269 |
| FIT208A7 (R5) | 1,032 | 3,828 | 1,032 | 0.270 |
| FIT208A8 (P) | 612 | 2,446 | 612 | 0.250 |

**Whitelist-corrected post-hoc** (`barcode_counts_wl.txt`, `737K-crdna-v1.txt`):

| Sample | Cell Ranger | Knee-point | Overlap | Jaccard |
|---|---|---|---|---|
| FIT208A3 (R1) | 1,421 | 3,114 | 1,421 | **0.456** |
| FIT208A4 (R2) | 996 | 1,758 | 996 | **0.567** |
| FIT208A5 (R3) | 82 | 146 | 82 | **0.562** |
| FIT208A6 (R4) | 265 | 449 | 265 | **0.590** |
| FIT208A7 (R5) | 1,032 | 1,717 | 1,032 | **0.601** |
| FIT208A8 (P) | 612 | 1,480 | 612 | **0.414** |

`only_cr = 0` in both modes: all Cell Ranger DNA barcodes are recovered. The
remaining Jaccard gap (i.e. the `only_knee` cells) reflects that our knee threshold
accepts more cells than Cell Ranger DNA. A full re-alignment with `RAW_MODE=false`
and `737K-crdna-v1.txt` is expected to bring the Jaccard above 0.7.

### Run

```bash
# Step 1: extract per-barcode read counts from BAMs
bash code/0_preprocessing/10x/custom/extract_barcode_counts.sh

# Step 2: (optional) post-hoc whitelist correction
python3 code/0_preprocessing/10x/custom/apply_whitelist_to_counts.py \
    --whitelist /srv/home/aste0033/projects/MPNST/scDNA/10X/737K-crdna-v1.txt \
    --outdir    /srv/home/aste0033/projects/MPNST/scDNA/10X/bwa_output

# Step 3: knee-point cell calling (raw or whitelist-corrected)
eval "$(micromamba shell hook --shell bash)" && micromamba activate R
Rscript code/0_preprocessing/10x/custom/knee_point_cell_calling.R          # raw
# or:
Rscript - <<'EOF'
COUNT_SUFFIX <- "_wl"
source("code/0_preprocessing/10x/custom/knee_point_cell_calling.R")
EOF

# Step 4: Cell Ranger DNA extraction + Jaccard comparison
Rscript code/0_preprocessing/10x/custom/extract_barcodes_from_rds.R        # raw
# or:
Rscript - <<'EOF'
KNEE_SUFFIX <- "_wl"
source("code/0_preprocessing/10x/custom/extract_barcodes_from_rds.R")
EOF

# Step 5: (optional) visual validation overlay
Rscript code/0_preprocessing/10x/custom/validate_cell_calling.R
```

After step 3, inspect `knee_plots[_wl].pdf`. If the automatic knee misses the
inflection for any sample, adjust `MANUAL_THRESHOLDS` in `knee_point_cell_calling.R`
and re-run.

### Output files (all under `/srv/home/aste0033/projects/MPNST/scDNA/10X/bwa_output/`)

| File | Description |
|---|---|
| `<sample>/barcode_counts.txt` | Raw per-barcode read counts (from BAM CB tags) |
| `<sample>/barcode_counts_wl.txt` | Whitelist-corrected counts (post-hoc) |
| `<sample>/valid_barcodes_knee.txt` | Knee-point called cells (raw) |
| `<sample>/valid_barcodes_knee_wl.txt` | Knee-point called cells (whitelist-corrected) |
| `<sample>/valid_barcodes_cellranger.txt` | Original Cell Ranger DNA cells |
| `knee_plots.pdf` / `knee_plots_wl.pdf` | 6-panel knee plots |
| `cell_calling_summary.tsv` / `cell_calling_summary_wl.tsv` | Per-sample summary |
| `cell_calling_comparison.tsv` / `cell_calling_comparison_wl.tsv` | Jaccard table |

### Subset BAM to valid cells (optional)

```bash
samtools view -D CB:valid_barcodes_cellranger.txt -b \
    -o possorted_bam.cells.bam possorted_bam.bam
samtools index possorted_bam.cells.bam
```

---

## Full execution order

```
1. bash code/0_preprocessing/10x/custom/fastq_to_bam_bwa.sh
   → produces markdup.bam per sample

2. bash code/0_preprocessing/10x/custom/run_bqsr_only.sh          (when SNV calling needed)
   → produces possorted_bam.bam per sample

3. bash code/0_preprocessing/10x/custom/extract_barcode_counts.sh
   → can run on markdup.bam before BQSR

4. python3 code/0_preprocessing/10x/custom/apply_whitelist_to_counts.py ...
   → only needed if RAW_MODE=true was used in step 1

5. Rscript code/0_preprocessing/10x/custom/knee_point_cell_calling.R
   → inspect knee_plots.pdf; adjust MANUAL_THRESHOLDS if needed

6. Rscript code/0_preprocessing/10x/custom/extract_barcodes_from_rds.R
   → Jaccard comparison table
```
