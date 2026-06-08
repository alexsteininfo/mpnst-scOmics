# DLP+ Preprocessing Pipeline

## Overview

Converts per-cell DLP+ FASTQs directly into per-sample deduplicated BAMs. `fastq_to_bam_dlpp.sh` aligns each cell with BWA-MEME, then uses the cell-to-region mapping from `scDNA_DLP_fastq.csv` to merge cells straight into per-sample output BAMs.

| Chip | Samples | FASTQs (total) | Cells in CSV | Run directory |
|---|---|---|---|---|
| LES4677A1 | R1, R2 | 4 415 | 2 450 | `220401_A01366_0166_AH2CNFDRX2/` |
| LES4677A2 | R3, R4 | 4 347 | 2 401 | `220228_A01366_0152_AHTTMYDRXY/` |
| LES4677A4 | R5, P  | 4 388 | 2 450 | `220909_A01366_0277_AHKF5CDRX2/` |

~55% of sequenced cells are retained per chip. See [Cell filtering](#cell-filtering) for why.

---

## Input requirements

Per-cell FASTQ pairs in the three run directories under `DLP_plus/`. FASTQs are already demultiplexed by chip position — no barcode extraction is needed.

Filename convention: `{CHIP}_{FLOWCELL}_{X}x_{Y}y_{S#}_{LANE}_R{1|2}_001.fastq.gz`  
Example: `LES4677A1_124574_14x_22y_S958_L001_R1_001.fastq.gz`

The cell ID (`LES4677A1_124574_14x_22y`) is extracted automatically from the filename by trimming from `_S{N}_` onward.

### Per-sample split

Cell-to-region assignment comes from `scDNA_DLP_fastq.csv` (path:
`/mnt/iribhm/people/mtarabichi/MPNST/pvanloo_20250618/Haixi/metadata/scDNA_DLP_fastq.csv`).
The CSV has three columns (`data`, `fastq`, `region`); the cell ID is extracted from
the fastq basename by stripping the `_S<N>_...` suffix. Confirmed mapping:

| Sample | Chip | Cells in CSV |
|---|---|---|
| R1 | LES4677A1 | 1 225 |
| R2 | LES4677A1 | 1 225 |
| R3 | LES4677A2 | 1 176 |
| R4 | LES4677A2 | 1 225 |
| R5 | LES4677A4 | 1 225 |
| P  | LES4677A4 | 1 225 |

`fastq_to_bam_dlpp.sh` parses the CSV at startup to build per-region cell lists,
then merges cells directly into per-sample BAMs during the merge step — no split
pass over a combined chip BAM is needed.

Only CSV cells (~2 450 per chip) are aligned. The ~1 950 non-CSV cells per chip
(chip-border or QC-failed wells) are skipped entirely at the FASTQ scan stage.

---

## Cell filtering

~45% of sequenced cells are absent from `scDNA_DLP_fastq.csv`. The loss occurs in
two independent stages, verified by comparing FASTQ filenames (x/y chip coordinates)
against the CSV for LES4677A1 (representative of all three chips):

| Stage | Cells dropped (LES4677A1) | Cause |
|---|---|---|
| Chip border exclusion | ~1 010 (~23%) | x < 14 or x > 62 excluded from all analysis |
| Interior QC failures  | ~955 (~22%)   | Empty wells, low coverage, or other QC failures |
| **Total dropped**     | **~1 965**    | **~55% retention** |

**Stage 1 — chip border exclusion:**  
The physical chip runs from x=01–72, but the 13 outer columns on each side
(x=01–13 and x=63–72) are excluded from Haixi's analysis. Border wells have
unreliable cell capture due to edge effects in the microfluidic channel and are
standard practice to exclude in DLP+ pipelines.

**Stage 2 — interior QC filtering:**  
Of the ~3 400 interior cells (x=14–62), ~28% are still absent from the CSV.
These are most likely empty wells — physically present positions on the chip that
captured no cell and therefore produced a FASTQ file with very few reads. There is
no separate QC report from Haixi's pipeline available to verify the exact criterion,
but read-count thresholding is the standard approach for DLP+ empty-well removal.

The CSV is therefore Haixi's post-QC cell list from the original CAMP analysis,
not a raw demultiplexing output.

---

## Pipeline

```
FASTQ pairs (one pair per cell, ~4 400 cells per chip)
  └── cell IDs auto-discovered from R1 filenames (no metadata file needed)
  └── scDNA_DLP_fastq.csv → cell-to-region map (parsed at startup)
        └── BWA-MEM  -R "@RG\tID:<cell_id>\tSM:<chip>..."
              └── samtools fixmate → sort → markdup   [per-cell BAM]
                    └── samtools merge (grouped by region)  [R1.bam … P.bam]
                          └── GATK BaseRecalibrator
                                └── GATK ApplyBQSR   [possorted_bam.bam]
```

**Step 1 — per-cell alignment and deduplication**  
Each cell's R1/R2 FASTQ pair is aligned independently with BWA-MEM. The read group (`-R`) tags every read with `RG:Z:<cell_id>` (chip position) and `SM:<chip>` (placeholder until per-sample split). `samtools fixmate` adds mate-score annotations; `samtools markdup` removes coordinate-level PCR duplicates *within* that cell's reads only. Running markdup per-cell before the merge prevents reads from two different cells at the same genomic position from being falsely called duplicates of each other.

**Step 2 — merge**  
Per-cell BAMs are merged by sample (not by chip) using the cell-to-region mapping from `scDNA_DLP_fastq.csv`. `samtools merge` produces one coordinate-sorted BAM per sample; per-cell RG tags are preserved.

**Step 3 — BQSR**  
`GATK BaseRecalibrator` builds a recalibration table using dbSNP138, Mills indels, and known indels from the GATK hg38 bundle. `ApplyBQSR` applies the correction. Uses the same resource files as the 10x pipeline.

---

## Output

```
<OUTDIR>/
  LES4677A1/
    R1.bam / R1.bam.bai    ← per-sample BAM, merged directly from per-cell BAMs
    R2.bam / R2.bam.bai
  LES4677A2/
    R3.bam / R3.bam.bai
    R4.bam / R4.bam.bai
  LES4677A4/
    R5.bam / R5.bam.bai
    P.bam  / P.bam.bai
```

Intermediate per-cell BAMs are deleted on success. The BQSR step (`run_bqsr_only.sh`) is separate and optional.

---

## Parallelism

| Parameter | Value | Notes |
|---|---|---|
| `CHIP_PARALLEL` | 1 | One chip at a time; each chip is ~4 400 cells |
| `CELL_PARALLEL` | 16 | Concurrent per-cell alignments within one chip (rolling pool) |
| `CELL_THREADS` | 2 | BWA-MEME threads per cell |
| `MERGE_THREADS` | 16 | samtools merge threads (alignment phase is done; all cores free) |
| **Alignment phase** | **32 cores** | 16 cells × 2 BWA threads; samtools fixmate/markdup add negligible overhead |
| **Merge phase** | **≤16 cores** | I/O-bound on NFS; more threads offer diminishing returns |

Per-cell BAMs are not indexed (`.bai` files are not created). `samtools merge` reads
BAMs sequentially and requires no index; skipping indexing saves ~10–20 min per chip.

---

## Runtime estimate

| Chip | Cells aligned | Cells merged | Median R1 size | Estimated time |
|---|---|---|---|---|
| LES4677A1 | 2 450 | 2 450 | 21 MB | ~3–4 h |
| LES4677A2 | 2 401 | 2 401 | 14 MB | ~2–3 h |
| LES4677A4 | 2 450 | 2 450 | 22 MB | ~3–4 h |
| **Total (sequential)** | | | | **~8–11 h** |

---

## Tags in the final BAM

| Tag | Source | Meaning |
|---|---|---|
| `RG:Z:<cell_id>` | BWA `-R` flag | Identifies the source cell (chip position, e.g. `LES4677A1_124574_14x_22y`) |
| `SM:<chip>` | BWA `-R` flag | Chip ID; not updated during the per-sample split (headers are preserved as-is) |
| `MC:Z:` | samtools fixmate | Mate cigar string |
| `ms:i:` | samtools fixmate | Mate score for markdup |

No per-read `CB:Z:` or `UB:Z:` tags are added. The `RG` tag is sufficient for ASCAT.sc and most CNA tools. If `alleleCounter` is needed for downstream SNP genotyping, add an artificial `UB:Z:` field post-hoc using the awk approach in Haixi's `modbam.R`.

---

## Differences from the 10x pipeline (`fastq_to_bam_bwa.sh`)

| Aspect | 10x | DLP+ |
|---|---|---|
| Barcode structure | 16 bp CB + 10 bp UMI in R1 | None; cells pre-demultiplexed by chip position |
| Barcode extraction | `extract_barcodes.py` (whitelist correction) | Not needed |
| Alignment mode | Stream all reads per sample in one BWA call | BWA called once per cell |
| Duplicate marking | `samtools markdup --barcode-tag CB` on merged stream | `samtools markdup` per-cell (before merge) |
| Merge step | Not needed | `samtools merge` after per-cell alignment |
| Per-read CB / UMI tags | `CB:Z:` and `UR:Z:` on every read | None (cell identity carried by `RG:Z:`) |
| Output granularity | One BAM per sample | One BAM per chip → split into per-sample BAMs |
| Parallelism model | N samples × M cores/sample | 1 chip × M cells × K threads/cell |

---

## Differences from Haixi's original DLP+ pipeline (`align_bwamem.R` + `modbam.R`)

| Aspect | Original (CAMP) | This pipeline |
|---|---|---|
| Language | R `system()` calls | Pure bash |
| Reference genome | `/camp/lab/vanloop/working/mtarabichi/reffiles/genome.fa` | `/mnt/iribhm/genomes/hg38-ngs/hg38.fa` |
| SAMtools version | 1.3.1 | 1.12 |
| Deduplication tool | Sambamba `--remove-duplicates` | `samtools markdup` |
| Output granularity | Per-sample BAM (after filtering by barcode list) | Per-sample BAM, merged directly from per-cell BAMs using CSV mapping |
| Merge step | Implicit / manual | Explicit `samtools merge` |
| BQSR | Not applied | Applied (GATK BaseRecalibrator + ApplyBQSR) |
| UB field | Added via `modbam.R` awk hack (row number, not a real UMI) | Not added; add post-hoc if alleleCounter needed |
| Parallelism | R `mclapply` (20 cores) | Bash background jobs, configurable |

---

## Note on UMIs

DLP+ does **not** use UMIs. Physical cell isolation on the microfluidic chip (one cell per well) means library preparation starts from a single cell's genome with no need for molecular barcoding. PCR duplicates are identified by genomic coordinate only — the same approach as bulk WGS.

The `UB:Z:` tag in `modbam.R` is **not a true UMI**. It was created artificially (assigning each read its SAM row number as a unique identifier) because `alleleCounter` required a UB field. It played no role in deduplication and has no biological meaning.
