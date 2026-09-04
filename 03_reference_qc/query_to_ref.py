#!/usr/bin/env python3
"""Map specific query (read) positions to their reference positions.

Useful for tracing where a known position within a read lands on the reference,
e.g. to check whether a read's internal coordinate corresponds to a boundary
or a known variant site. Hard-clip offsets are accounted for so that positions
are given in original query coordinates (as reported by Dorado/basecaller),
not the soft-clipped fragment seen by the aligner.

Usage: python3 query_to_ref.py <bam> <pos1> [pos2 ...]

  pos: 0-based query positions in original (hard-clip-corrected) coordinates.
  Iterates all records in the BAM; pass a region-filtered BAM to limit output.
"""
import sys
import pysam

bam_path = sys.argv[1]
targets = [int(x) for x in sys.argv[2:]]

bam = pysam.AlignmentFile(bam_path, "rb")

for read in bam.fetch(until_eof=True):
    if read.is_unmapped:
        continue
    pairs = read.get_aligned_pairs(matches_only=False)
    leading_hardclip = 0
    if read.cigartuples and read.cigartuples[0][0] == 5:
        leading_hardclip = read.cigartuples[0][1]

    qpos_to_rpos = {}
    for qpos, rpos in pairs:
        if qpos is not None:
            qpos_to_rpos[qpos + leading_hardclip] = rpos

    print(f"=== read {read.query_name} flag={read.flag} ref_name={read.reference_name} ref_start={read.reference_start} ===")
    for t in targets:
        rpos = qpos_to_rpos.get(t, "not in this record's query span")
        print(f"  query pos {t} -> reference pos {rpos}")
