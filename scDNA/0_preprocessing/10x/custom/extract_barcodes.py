#!/usr/bin/env python3
"""
Extract 10x cell barcodes from R1 reads and stream barcode-tagged R2 reads to stdout.

10x Chromium scDNA library structure:
  R1: [16 bp cell barcode][8 bp UMI][...]  -- identifies the cell
  R2: genomic DNA sequence                  -- what gets aligned

Default mode (--whitelist required):
  Barcode is corrected against the 10x whitelist (exact match first, then
  1-mismatch correction including N). Reads with uncorrectable barcodes are
  discarded. CB tag in output = corrected barcode.

Raw mode (--raw):
  Every read pair is kept. CB tag = raw R1[0:16] with no whitelist lookup.
  Use this to recover cells whose barcodes are >1 mismatch from the whitelist
  (e.g. Cell Ranger DNA rescued cells). Cell calling is done downstream by
  read-depth threshold on the raw barcode distribution.

Output FASTQ comment format (consumed by bwa mem -C):
  @READNAME CB:Z:<barcode>\tUR:Z:<raw_umi>

Usage (whitelist mode):
  python3 extract_barcodes.py \
      --r1 run1/S1_L001_R1.fastq.gz ... \
      --r2 run1/S1_L001_R2.fastq.gz ... \
      --whitelist 737K-august-2016.txt \
  | bwa mem -t 24 -C genome.fa -

Usage (raw mode):
  python3 extract_barcodes.py --raw \
      --r1 run1/S1_L001_R1.fastq.gz ... \
      --r2 run1/S1_L001_R2.fastq.gz ... \
  | bwa mem -t 24 -C genome.fa -
"""
import sys
import gzip
import argparse


def load_whitelist(path):
    opener = gzip.open if path.endswith('.gz') else open
    with opener(path, 'rt') as f:
        return set(line.strip().split()[0] for line in f if line.strip())


def build_correction_map(whitelist):
    """Return dict mapping each 1-mismatch neighbour to its canonical barcode.
    Ambiguous neighbours (matching >1 whitelist barcode) map to None.
    N is included as a possible corrupted base to rescue low-quality first-cycle
    base calls (a systematic artifact on some Illumina runs)."""
    corr = {}
    for bc in whitelist:
        for i in range(len(bc)):
            for base in 'NACGT':
                if base == bc[i]:
                    continue
                nb = bc[:i] + base + bc[i+1:]
                if nb in whitelist:
                    continue
                if nb in corr:
                    corr[nb] = None   # ambiguous: two valid corrections
                else:
                    corr[nb] = bc
    return corr


def open_fq(path):
    return gzip.open(path, 'rt') if path.endswith('.gz') else open(path, 'rt')


def iter_fastq(fh):
    while True:
        hdr = fh.readline()
        if not hdr:
            return
        seq  = fh.readline().rstrip()
        fh.readline()   # + line
        qual = fh.readline().rstrip()
        yield hdr.rstrip(), seq, qual


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--r1', nargs='+', required=True,
                    help='R1 fastq files (barcode+UMI read), one per lane/run')
    ap.add_argument('--r2', nargs='+', required=True,
                    help='R2 fastq files (genomic read), matching order with --r1')
    ap.add_argument('--whitelist', default=None,
                    help='10x barcode whitelist (.txt or .txt.gz); required unless --raw')
    ap.add_argument('--raw', action='store_true',
                    help='Raw mode: keep all reads, tag with unfiltered R1[0:bc_len] as CB. '
                         'No whitelist needed. Use for cell calling by read-depth threshold.')
    ap.add_argument('--bc-len', type=int, default=16,
                    help='Cell barcode length in R1 (default: 16)')
    ap.add_argument('--umi-len', type=int, default=10,
                    help='UMI length in R1 following barcode (default: 10 for 10x CNV)')
    args = ap.parse_args()

    if len(args.r1) != len(args.r2):
        sys.exit(f'ERROR: --r1 ({len(args.r1)} files) and --r2 ({len(args.r2)} files) '
                 f'must have the same number of files.')

    if not args.raw and args.whitelist is None:
        sys.exit('ERROR: --whitelist is required unless --raw is set.')

    bc_len  = args.bc_len
    umi_len = args.umi_len
    out     = sys.stdout

    if args.raw:
        sys.stderr.write('Raw mode: keeping all reads with unfiltered barcode as CB tag.\n')
        whitelist  = None
        correction = None
    else:
        sys.stderr.write(f'Loading whitelist from {args.whitelist} ...\n')
        whitelist  = load_whitelist(args.whitelist)
        correction = build_correction_map(whitelist)
        sys.stderr.write(f'  {len(whitelist):,} barcodes loaded\n')

    total = kept = n_exact = n_corrected = n_raw = 0

    for r1_path, r2_path in zip(args.r1, args.r2):
        sys.stderr.write(f'Processing {r1_path}\n')
        with open_fq(r1_path) as f1, open_fq(r2_path) as f2:
            for (h1, s1, _), (h2, s2, q2) in zip(iter_fastq(f1), iter_fastq(f2)):
                total += 1
                bc_raw  = s1[:bc_len]
                umi_raw = s1[bc_len:bc_len + umi_len]

                if args.raw:
                    bc_corr = bc_raw
                    n_raw += 1
                elif bc_raw in whitelist:
                    bc_corr = bc_raw
                    n_exact += 1
                elif bc_raw in correction and correction[bc_raw]:
                    bc_corr = correction[bc_raw]
                    n_corrected += 1
                else:
                    continue   # barcode not in whitelist; discard read

                kept += 1
                # Strip any existing comment from R2 header; add CB/UR tags as comment.
                # BWA -C appends FASTQ comments verbatim to SAM output as extra fields,
                # so tab-separated CB:Z:/UR:Z: values become proper SAM tags.
                read_name = h2[1:].split()[0]
                out.write(f'@{read_name} CB:Z:{bc_corr}\tUR:Z:{umi_raw}\n'
                          f'{s2}\n+\n{q2}\n')

    discarded = total - kept
    if args.raw:
        sys.stderr.write(
            f'Finished (raw mode): {total:,} read pairs processed\n'
            f'  All reads kept      : {kept:,} ({100*kept/max(total,1):.1f}%)\n'
        )
    else:
        sys.stderr.write(
            f'Finished: {total:,} read pairs processed\n'
            f'  Exact barcode match : {n_exact:,} ({100*n_exact/max(total,1):.1f}%)\n'
            f'  Corrected (1 mm)    : {n_corrected:,} ({100*n_corrected/max(total,1):.1f}%)\n'
            f'  Discarded           : {discarded:,} ({100*discarded/max(total,1):.1f}%)\n'
            f'  Kept                : {kept:,} ({100*kept/max(total,1):.1f}%)\n'
        )


if __name__ == '__main__':
    main()
