#!/usr/bin/env python3
"""
barcode_stats.py  —  Extract per-barcode region stats from a 10x BAM file.

Usage:
    python barcode_stats.py --bam outs/possorted_genome_bam.bam
    python barcode_stats.py --bam outs/possorted_genome_bam.bam --out my_stats.tsv
"""

import argparse
import csv
import sys
from collections import defaultdict

import pysam

def parse_args():
    parser = argparse.ArgumentParser(
        description="Compute per-barcode exonic/intronic/intergenic read fractions."
    )
    parser.add_argument(
        "--bam",
        required=True,
        help="Path to the Cell Ranger BAM file (must be sorted and indexed)."
    )
    parser.add_argument(
        "--out",
        default="barcode_stats.tsv",
        help="Output TSV path (default: barcode_stats.tsv)."
    )
    return parser.parse_args()


def process_bam(bam_path):
    EXCLUDE_FLAGS = 0x4 | 0x100 | 0x400 | 0x800
    counts = defaultdict(lambda: {"total": 0, "E": 0, "N": 0, "I": 0})
    reads_seen = 0
    reads_processed = 0

    print(f"Opening {bam_path} ...", file=sys.stderr)

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam:
            reads_seen += 1

            if reads_seen % 5_000_000 == 0:
                print(f"\r  {reads_seen:,} reads scanned, {reads_processed:,} counted ...",
                      end="", flush=True, file=sys.stderr)

            if read.flag & EXCLUDE_FLAGS:
                continue

            if not read.has_tag("CB"):
                continue
            cb = read.get_tag("CB")

            if not read.has_tag("RE"):
                continue
            re = read.get_tag("RE")

            counts[cb]["total"] += 1

            if re in ("E", "N", "I"):
                counts[cb][re] += 1

            reads_processed += 1

    print(f"\r  Done. {reads_seen:,} reads scanned, {reads_processed:,} counted.     ",
          file=sys.stderr)

    return counts


def write_tsv(counts, out_path):
    header = [
        "barcode", "total", "exonic", "intronic", "intergenic",
        "frac_exonic", "frac_intronic", "frac_intergenic",
    ]

    sorted_barcodes = sorted(counts.items(), key=lambda x: -x[1]["total"])

    with open(out_path, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(header)

        for cb, c in sorted_barcodes:
            tot = c["total"]
            e, n, ig = c["E"], c["N"], c["I"]
            writer.writerow([
                cb, tot, e, n, ig,
                f"{e  / tot:.4f}",
                f"{n  / tot:.4f}",
                f"{ig / tot:.4f}",
            ])

    print(f"Written {len(counts):,} barcodes → {out_path}", file=sys.stderr)


if __name__ == "__main__":
    args   = parse_args()
    counts = process_bam(args.bam)

    if not counts:
        print("ERROR: No valid reads found. Check your BAM path and tags.", file=sys.stderr)
        sys.exit(1)

    write_tsv(counts, args.out)
