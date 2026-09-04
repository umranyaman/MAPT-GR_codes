#!/usr/bin/env python3
"""Check for large CIGAR deletions spanning the masked gap in flanking reads.

A deletion >=1000 bp whose coordinates overlap the gap region, seen in primary
reads, would indicate either a genuine structural variant or a mis-alignment
artefact at the boundary. Reads are fetched from a 2 kb window around the gap.

Usage: python3 check_flanking_deletion.py <bam> <chrom> <gap_start> <gap_end>

  gap_start / gap_end: 0-based coordinates of the masked mouse interval
                       (chr11:104,066,949-104,223,491 in the MAPT-GR reference).
"""
import sys
import pysam

bam_path, chrom, gap_start, gap_end = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
window = 2000
bam = pysam.AlignmentFile(bam_path, "rb")

found = False
for read in bam.fetch(chrom, max(0, gap_start - window), gap_end + window):
    if read.is_secondary or read.is_supplementary or read.is_unmapped:
        continue
    rpos = read.reference_start
    for op, length in read.cigartuples:
        if op == 2 and length >= 1000:  # D = deletion
            del_start = rpos
            del_end = rpos + length
            if del_start < gap_end and del_end > gap_start:
                print(f"  {read.query_name}: DELETION of {length}bp at {chrom}:{del_start}-{del_end} (MAPQ={read.mapping_quality})")
                found = True
        if op in (0, 7, 8, 2):
            rpos += length

if not found:
    print("  No large deletion (>=1000bp) found in any flanking read spanning this gap.")
