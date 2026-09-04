#!/usr/bin/env python3
"""Find reads in a mouse interval that have a supplementary alignment on chrTau.

Split-read support across the insert-mouse junction appears as a primary
alignment on the mouse chromosome and a supplementary alignment (SA tag)
on chrTau, or vice versa. This script queries one side of that junction:
reads fetched from the specified mouse region that carry an SA tag referencing
chrTau.

Usage: python3 check_chrTau_crosslink.py <bam> <chrom> <start> <end>

  Coordinates are 0-based half-open (pysam fetch convention).
"""
import sys
import pysam

bam_path, chrom, start, end = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
bam = pysam.AlignmentFile(bam_path, "rb")

hits = []
for read in bam.fetch(chrom, start, end):
    if read.is_secondary or read.is_unmapped:
        continue
    if read.has_tag("SA"):
        sa = read.get_tag("SA")
        if "chrTau" in sa:
            hits.append((read.query_name, read.reference_start, read.reference_end, sa))

print(f"Reads in {chrom}:{start}-{end} with a supplementary alignment on chrTau: {len(hits)}")
for name, s, e, sa in hits:
    print(f"  {name}  mouse_pos={s}-{e}  SA={sa}")
