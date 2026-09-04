#!/usr/bin/env python3
"""Inspect reads near the MAPT-GR insert boundaries for alignment anomalies.

For each of the two boundary positions (start-of-gap and end-of-gap on the
mouse reference), reports reads with large soft/hard clips (>100 bp) whose
clip end is within 500 bp of the boundary, large insertions in the CIGAR
(>=50 bp), and split-read SA tags. These signatures indicate reads that span
the insert-mouse junction and were not cleanly resolved by the aligner.

Usage: python3 check_boundary_signature.py <bam> <chrom> <gap_start> <gap_end>

  gap_start / gap_end: 1-based boundary coordinates of the masked mouse interval
                       (chr11:104,066,949-104,223,491 in the MAPT-GR reference).
"""
import sys
import pysam

bam_path, chrom, b1, b2 = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
window = 500  # look this far either side of each boundary
bam = pysam.AlignmentFile(bam_path, "rb")

for label, bpos in [("start-of-gap", b1), ("end-of-gap", b2)]:
    print(f"=== reads near {label} (chr={chrom}, pos={bpos}, +/-{window}bp) ===")
    for read in bam.fetch(chrom, max(0, bpos - window), bpos + window):
        if read.is_secondary or read.is_unmapped:
            continue
        ref_end = read.reference_end
        clip_at_start = read.cigartuples[0][1] if read.cigartuples[0][0] in (4, 5) else 0
        clip_at_end = read.cigartuples[-1][1] if read.cigartuples[-1][0] in (4, 5) else 0
        large_insertions = [(op, length) for op, length in read.cigartuples if op == 1 and length >= 50]
        sa = read.get_tag("SA") if read.has_tag("SA") else None
        flags = []
        if abs(read.reference_start - bpos) < window and clip_at_start > 100:
            flags.append(f"SOFT/HARD-CLIP AT START ({clip_at_start}bp, ref_start={read.reference_start})")
        if ref_end and abs(ref_end - bpos) < window and clip_at_end > 100:
            flags.append(f"SOFT/HARD-CLIP AT END ({clip_at_end}bp, ref_end={ref_end})")
        if large_insertions:
            flags.append(f"LARGE INSERTION(S) in CIGAR: {large_insertions}")
        if sa:
            flags.append(f"has SA tag (split read): {sa}")
        if flags:
            print(f"  {read.query_name} (supplementary={read.is_supplementary}, MAPQ={read.mapping_quality}): " + "; ".join(flags))
    print()
