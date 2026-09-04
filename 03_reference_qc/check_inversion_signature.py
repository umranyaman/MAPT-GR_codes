#!/usr/bin/env python3
"""Check for inversion-breakpoint signatures in a BAM region.

An inversion breakpoint produces split reads where the primary alignment and
its supplementary alignment are on the same contig but opposite strands. This
script counts such events within a specified region. A cluster of these in a
narrow span is consistent with a real inversion; reads scattered across the
region are likely chimeric artefacts.

Also reports the forward/reverse strand balance of primary reads as a basic
sanity check.

Usage: python3 check_inversion_signature.py <bam> <chrom> <start> <end>

  Coordinates are 0-based half-open (pysam fetch convention).
"""
import sys
import pysam

bam_path, chrom, start, end = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
bam = pysam.AlignmentFile(bam_path, "rb")

n_fwd, n_rev = 0, 0
inversion_signature_hits = []

for read in bam.fetch(chrom, start, end):
    if read.is_secondary or read.is_supplementary or read.is_unmapped:
        continue
    if read.is_reverse:
        n_rev += 1
    else:
        n_fwd += 1

    if read.has_tag("SA"):
        sa = read.get_tag("SA")
        primary_strand = "-" if read.is_reverse else "+"
        for entry in sa.strip(";").split(";"):
            fields = entry.split(",")
            if len(fields) < 3:
                continue
            sa_rname, sa_pos, sa_strand = fields[0], fields[1], fields[2]
            if sa_rname == chrom and sa_strand != primary_strand:
                inversion_signature_hits.append((read.query_name, read.reference_start, primary_strand, sa_pos, sa_strand))

total = n_fwd + n_rev
print(f"=== {bam_path} :: {chrom}:{start}-{end} ===")
if total == 0:
    print("Primary reads: 0 (no coverage in region)")
else:
    print(f"Primary reads: {total} total, {n_fwd} forward ({100*n_fwd/total:.1f}%), {n_rev} reverse ({100*n_rev/total:.1f}%)")
print(f"\nSplit reads with opposite-strand SA on the same contig (inversion-breakpoint signature): {len(inversion_signature_hits)}")
if inversion_signature_hits:
    print("  read_name  primary_pos  primary_strand  ->  SA_pos  SA_strand")
    for name, pos, pstrand, sapos, sastrand in inversion_signature_hits:
        print(f"  {name}  {pos}  {pstrand}  ->  {sapos}  {sastrand}")
    positions = [pos for _, pos, _, _, _ in inversion_signature_hits]
    print(f"\n  position range of these hits: {min(positions)}-{max(positions)} (span={max(positions)-min(positions)}bp)")
    print("  -> tightly clustered in a small span = consistent with a real inversion breakpoint")
    print("  -> widely scattered across the region = not consistent with a real inversion (likely noise/chimeric artifacts)")
else:
    print("  -> none found: no evidence of an inversion-breakpoint signature in this region")
