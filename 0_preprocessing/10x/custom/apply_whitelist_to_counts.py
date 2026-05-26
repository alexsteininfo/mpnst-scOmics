#!/usr/bin/env python3
"""
Apply Hamming-1 whitelist correction to barcode_counts.txt files.

Reads the raw barcode count table produced by extract_barcode_counts.sh
(format: "  <count> <barcode>", space-padded, sorted by count descending),
corrects each barcode against the scDNA whitelist (exact match, then 1-mismatch),
sums counts for barcodes corrected to the same whitelist entry, and writes
barcode_counts_wl.txt in the same format.

This is equivalent to the correction done during alignment in whitelist mode,
applied post-hoc so the full pipeline does not need to be re-run.

Usage:
    python3 apply_whitelist_to_counts.py \
        --whitelist /path/to/737K-crdna-v1.txt \
        --outdir    /path/to/bwa_output \
        --samples   FIT208A3 FIT208A4 FIT208A5 FIT208A6 FIT208A7 FIT208A8
"""
import sys
import gzip
import argparse
import collections


def load_whitelist(path):
    opener = gzip.open if path.endswith('.gz') else open
    with opener(path, 'rt') as f:
        return set(line.strip().split()[0] for line in f if line.strip())


def build_correction_map(whitelist):
    """Map each 1-mismatch neighbour → canonical whitelist barcode.
    Ambiguous neighbours (matching >1 whitelist barcode) map to None."""
    corr = {}
    for bc in whitelist:
        for i in range(len(bc)):
            for base in 'NACGT':
                if base == bc[i]:
                    continue
                nb = bc[:i] + base + bc[i+1:]
                if nb in whitelist:
                    continue  # exact match exists; not a 1-mismatch neighbour
                if nb in corr:
                    corr[nb] = None  # ambiguous
                else:
                    corr[nb] = bc
    return corr


def correct_counts(in_path, out_path, whitelist, correction_map):
    corrected = collections.Counter()
    n_exact = n_hamming1 = n_discarded = 0

    with open(in_path) as fh:
        for line in fh:
            line = line.rstrip('\n')
            if not line.strip():
                continue
            parts = line.split()
            count, bc = int(parts[0]), parts[1]

            if bc in whitelist:
                corrected[bc] += count
                n_exact += 1
            elif bc in correction_map and correction_map[bc] is not None:
                corrected[correction_map[bc]] += count
                n_hamming1 += 1
            else:
                n_discarded += 1

    # Write in same format as uniq -c output (right-aligned count, space, barcode)
    sorted_bcs = sorted(corrected.items(), key=lambda x: -x[1])
    with open(out_path, 'w') as fh:
        for bc, cnt in sorted_bcs:
            fh.write(f"{cnt:>7} {bc}\n")

    total_in = n_exact + n_hamming1 + n_discarded
    kept_pct = 100 * (n_exact + n_hamming1) / total_in if total_in else 0
    print(f"  exact={n_exact:,}  hamming1={n_hamming1:,}  discarded={n_discarded:,}  "
          f"kept={kept_pct:.1f}%  unique_wl_barcodes={len(corrected):,}")
    return len(corrected)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--whitelist', required=True)
    ap.add_argument('--outdir',    required=True)
    ap.add_argument('--samples',   nargs='+',
                    default=['FIT208A3','FIT208A4','FIT208A5',
                             'FIT208A6','FIT208A7','FIT208A8'])
    args = ap.parse_args()

    import os
    print(f"Loading whitelist: {args.whitelist}")
    whitelist = load_whitelist(args.whitelist)
    print(f"  {len(whitelist):,} barcodes")
    print("Building 1-mismatch correction map...")
    correction_map = build_correction_map(whitelist)
    print(f"  {len(correction_map):,} neighbours mapped\n")

    for sample in args.samples:
        in_path  = os.path.join(args.outdir, sample, 'barcode_counts.txt')
        out_path = os.path.join(args.outdir, sample, 'barcode_counts_wl.txt')
        if not os.path.exists(in_path):
            print(f"[{sample}] SKIP: {in_path} not found")
            continue
        print(f"[{sample}]")
        correct_counts(in_path, out_path, whitelist, correction_map)
        print(f"  → {out_path}")


if __name__ == '__main__':
    main()
