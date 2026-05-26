# CellRanger DNA 1.0.0 — Patch Log

**Installation:** `/srv/home/aste0033/CellRanger-DNA/cellranger-dna-1.0.0/`

## Background

CellRanger DNA 1.0.0 (2018, Python 2.7) was run using the Python 2.7 interpreter
bundled with CellRanger 3.1.0 (`/opt/local/easybuild/software/CellRanger/3.1.0/
miniconda-cr-cs/`). That interpreter ships with newer library versions than CellRanger DNA
1.0.0 was written against, causing a cascade of API incompatibilities that required
manual patches.

All edits are backward-compatible with the intended Python 2.7 semantics.

---

## Fix 1 — `has_key()` removed in pysam 0.14

**Error:** `AttributeError: 'AlignmentHeader' object has no attribute 'has_key'`

pysam 0.14 changed `AlignmentHeader` from a plain dict to a custom object. The dict
method `has_key(x)` was removed; the idiomatic replacement is `x in obj`.

**Files patched** (all `has_key(x)` → `x in obj` replacements):

| File | Notes |
|---|---|
| `tenkit/lib/python/tenkit/bam.py` | 5 occurrences |
| `tenkit/lib/python/tenkit/alarms.py` | |
| `tenkit/lib/python/tenkit/supernova.py` | |
| `tenkit/lib/python/tenkit/cache.py` | |
| `tenkit/lib/python/tenkit/dict_utils.py` | |
| `tenkit/lib/python/tenkit/hdf5.py` | |
| `lib/python/longranger/cnv/copy_number_tools.py` | |
| `lib/python/longranger/cnv/contig_manager.py` | |
| `lib/python/crdna/singlecell_dna_cnv/copy_number_tools.py` | |

---

## Fix 2 — Immutable `AlignmentHeader` (pysam 0.14)

**Error:** `TypeError: 'AlignmentHeader' object does not support item assignment`

pysam 0.14 made `AlignmentHeader` immutable. Code that mutated it directly
(e.g. `header['RG'] = ...`) must first call `.to_dict()` to get a mutable dict copy.

**Files patched:**

### `tenkit/lib/python/tenkit/bam.py`
```python
# Before (line ~70):
header = template.header

# After:
_hdr = template.header
header = _hdr.to_dict() if hasattr(_hdr, 'to_dict') else _hdr
```

### `lib/python/crdna/bam.py` — `BamTemplateShim.__init__`
```python
# Before:
self.header = template.header
if not keep_comments and 'CO' in self.header:
    del self.header['CO']

# After:
self.header = template.header.to_dict() if hasattr(template.header, 'to_dict') else template.header
self.references = template.references
if not keep_comments and 'CO' in self.header:
    self.header = {k: v for k, v in self.header.items() if k != 'CO'}
```

---

## Fix 3 — Pandas computation import

**Error:** `AttributeError: module 'pandas.core' has no attribute 'expressions'`

The internal pandas API for disabling numexpr moved between versions.

**File:** `tenkit/lib/python/tenkit/pandas/__init__.py`

```python
# Before:
from pandas.core.computation.expressions import set_use_numexpr
set_use_numexpr(False)

# After:
try:
    from pandas.core.computation.expressions import set_use_numexpr
    set_use_numexpr(False)
except (ImportError, AttributeError):
    pass
```

---

## Fix 4 — Missing matplotlib and statsmodels

**Error:** `ImportError: No module named matplotlib`

Neither `matplotlib` nor `statsmodels` is installed in CellRanger 3.1.0's Python 2.7
environment. The plotting utilities in CellRanger DNA 1.0.0 import them at module
level, causing every stage that imports `plot_utils` to crash at startup.

**File:** `lib/python/plot_utils.py`

All top-level `import matplotlib` and `from statsmodels...` calls were wrapped in
`try/except ImportError` blocks. A module-level flag `_HAS_MATPLOTLIB` is set so
that plotting functions can skip gracefully when the library is absent.

---

## Fix 5 — Boolean index mismatch in clustering

**Error:** `IndexError: boolean index did not match indexed array along dimension 0`

In `compute_heterogeneity()`, `windows` has shape `(ncells + ninodes,)` (it includes
internal tree nodes), but `cell_filter` is a boolean mask of length `ncells` only.

**File:** `lib/python/crdna/clustering.py`, line 122

```python
# Before:
blur_window = int(np.max(windows[cell_filter]) / 2)

# After:
blur_window = int(np.max(windows[:ncells][cell_filter]) / 2)
```

---

## Fix 6 — Missing `dlconverter` binary (DLOUPE_PREPROCESS)

**Error:** `OSError: [Errno 2] No such file or directory` when calling `dlconverter`

The `dlconverter` binary generates `.dloupe` Loupe Browser files. It is a
closed-source binary that is **not** included in the CellRanger DNA 1.0.0
pre-built tarball. The original `except` clause only caught `CalledProcessError`,
not `OSError`, so a missing binary crashed the pipeline rather than failing
gracefully. The `.dloupe` file is a visualization artifact; the analysis outputs
(`possorted_bam.bam`, `per_cell_summary_metrics.csv`, `cnv_data.h5`) are
produced by later stages.

**File:** `mro/stages/dloupe/dloupe_preprocess/__init__.py`

```python
# Before:
except subprocess.CalledProcessError, e:
    outs.output_for_dloupe = None
    martian.throw("Could not generate .dloupe file: \n%s" % e.output)

# After:
except (subprocess.CalledProcessError, OSError), e:
    outs.output_for_dloupe = None
    martian.log_info("dlconverter not available or failed, skipping .dloupe generation: %s" % str(e))
```

---

## Fix 7 — Wrong species key lookup in REPORT_SINGLECELL

**Symptom:** 0 cells in all downstream outputs despite correct cell detection.

**Root cause:** The hg38 reference built with `cellranger-dna mkref` has
`species_prefixes: ["chr"]` in `fasta/contig-defs.json`. This value is meant as a
chromosome-name prefix (chr1, chr2, …), **not** a species name. However,
`contig_manager.list_species()` returns it as-is, so `species_list = ["chr"]`.

`REPORT_SINGLECELL` then calls `args.cell_barcodes.get("chr", {})` to retrieve the
list of called cells. But `DETECT_CELL_BARCODES` stores cells under the key `""`
(empty string, returned by `species_from_contig` for standard `chrN` contigs that
have no underscore). The `"chr"` key is never present, so `cell_index` is empty and
every barcode gets `cell_id = "None"`. With 0 cells, `REPORT_BASIC_SUMMARY` produces
an empty `per_cell_summary_metrics.csv` and `MAKE_WEBSUMMARY` crashes.

**File:** `mro/stages/reporter/report_singlecell/__init__.py`

Three locations were patched to fall back to `""` when the species prefix is not a
key in `args.cell_barcodes`:

```python
# Cell index construction — before:
for sp in species_list:
    bc_list = args.cell_barcodes.get(sp, {}).keys()

# After:
for sp in species_list:
    cb_key = sp if sp in args.cell_barcodes else ""
    bc_list = args.cell_barcodes.get(cb_key, {}).keys()
```

```python
# barnyard_hits row — before:
bh_row.append( int(s in args.cell_barcodes and bc in args.cell_barcodes[s]) )

# After:
cb_key = s if s in args.cell_barcodes else ""
bh_row.append( int(cb_key in args.cell_barcodes and bc in args.cell_barcodes[cb_key]) )
```

```python
# barnyard row read counts and is_chr_cell_barcode — before:
[num_per_species[s] for s in species_list]
[len(contigs_per_species[s]) for s in species_list]
barnyard_row.append( int((speci in args.cell_barcodes) and (bc in args.cell_barcodes[speci])) )

# After:
[num_per_species[s if s in num_per_species else ""] for s in species_list]
[len(contigs_per_species[s if s in contigs_per_species else ""]) for s in species_list]
cb_key = speci if speci in args.cell_barcodes else ""
barnyard_row.append( int((cb_key in args.cell_barcodes) and (bc in args.cell_barcodes[cb_key])) )
```

**One-time data patch for FIT208A3:** Because `REPORT_SINGLECELL` had already
completed and rerunning it on 1.4 billion reads would take hours, the cached
`barnyard.csv` was patched directly with a Python script
(`/tmp/patch_barnyard.py`). The script loaded the 2220 cell barcodes from
`DETECT_CELL_BARCODES/_outs`, assigned `cell_id = "chr_cell_N"` (alphabetical
order, matching what the fixed code would produce), and wrote the file back.
The `_complete` markers for `REPORT_BASIC_SUMMARY`, `MAKE_ALARMS`, and
`MAKE_WEBSUMMARY` were then deleted so those stages rerun with the corrected data.

---

## Fix 8 — Empty `cell_count` crash in MAKE_WEBSUMMARY

**Error:** `IndexError: list index out of range` at `notcell_count.insert(0, cell_count[-1])`

When all barcodes have `cell_id == "None"` (triggered by Fix 7's root cause),
`cell_count` is an empty list and `cell_count[-1]` raises `IndexError`.

**File:** `mro/stages/reporter/make_websummary/__init__.py`

```python
# Before:
notcell_count.insert(0, cell_count[-1])
notcell_reads_depth.insert(0, cell_reads_depth[-1])

# After:
if cell_count:
    notcell_count.insert(0, cell_count[-1])
    notcell_reads_depth.insert(0, cell_reads_depth[-1])
```

---

---

## Fix 9 — VDR-deleted REPORT_BASIC summary.json (FIT208A3 one-time data patch)

**Error:** `IOError: [Errno 2] No such file or directory: '.../REPORT_BASIC/fork0/join-.../files/summary.json'`

**Root cause:** Martian's VDR (Volatile Data Removal) system automatically deleted
REPORT_BASIC's intermediate output files (575 files, 195 MB) after the pipeline had
apparently completed, to free disk space. This happened at `2026-05-15 14:58:05`.
When we subsequently cleared REPORT_BASIC_SUMMARY's `_complete` markers to force it
to rerun with the corrected `barnyard.csv`, REPORT_BASIC_SUMMARY could no longer read
its required input.

**What VDR deleted:** `REPORT_BASIC/fork0/join-.../files/summary.json` (and all chunk
outputs). This file contains sequencing quality statistics: total reads, Q30 base
fractions, correct barcode rate.

**Fix:** The REPORT_BASIC join step prints the complete statistics dict to its
`_stdout` log before writing the file. That log survived VDR. The file was
reconstructed directly from the log:

```bash
# Extract the JSON printed to stdout and write it back
python3 -c "
import sys, json
with open('REPORT_BASIC/fork0/join-.../\_stdout') as f:
    content = f.read()
start = content.find('{')
end = content.rfind('}') + 1
d = json.loads(content[start:end])
with open('REPORT_BASIC/fork0/join-.../files/summary.json', 'w') as f:
    json.dump(d, f, indent=4)
"
```

Reconstructed values (verified correct):
- `num_reads`: 1,414,897,512
- `correct_bc_rate`: 0.8620445308974436
- `r1_tot_bases`: 95,505,582,060 / `r1_q30_bases`: 83,093,831,413
- `r2_tot_bases`: 106,824,762,156 / `r2_q30_bases`: 84,411,191,597

**No code change** — this was a data-only repair for FIT208A3. For future samples,
VDR will not be a problem if the pipeline runs straight through without cache
manipulation.

---

## Summary table

| Fix | File | Error type |
|---|---|---|
| 1 | `tenkit/bam.py` + 8 others | `AttributeError: has_key` |
| 2 | `tenkit/bam.py`, `crdna/bam.py` | `TypeError: item assignment` |
| 3 | `tenkit/pandas/__init__.py` | `AttributeError: expressions` |
| 4 | `lib/python/plot_utils.py` | `ImportError: matplotlib` |
| 5 | `lib/python/crdna/clustering.py` | `IndexError: boolean index` |
| 6 | `mro/stages/dloupe/dloupe_preprocess/__init__.py` | `OSError: dlconverter` |
| 7 | `mro/stages/reporter/report_singlecell/__init__.py` | 0 cells in output |
| 8 | `mro/stages/reporter/make_websummary/__init__.py` | `IndexError: cell_count[-1]` |
| 9 | Data patch — REPORT_BASIC `summary.json` reconstructed from `_stdout` | `IOError: summary.json` (VDR) |
