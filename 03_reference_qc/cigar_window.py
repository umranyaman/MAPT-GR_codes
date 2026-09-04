#!/usr/bin/env python3
"""Report CIGAR indels within a query-coordinate window, with reference coordinates.

Operates in query space rather than reference space, which is necessary when
the region of interest is defined by its position within a read (e.g. an insert
whose query length is known) rather than by reference coordinates. Leading
soft/hard clips are accounted for so that query positions are relative to the
original read, not the clipped fragment.

Usage: python3 cigar_window.py <bam> <query_start> <query_end> [min_indel_len]

  query_start / query_end: 0-based positions in original query coordinates.
  min_indel_len: minimum indel length to report (default 1).
"""
import sys
import pysam

bam_path = sys.argv[1]
q_start = int(sys.argv[2])
q_end = int(sys.argv[3])
min_indel = int(sys.argv[4]) if len(sys.argv) > 4 else 1

bam = pysam.AlignmentFile(bam_path, "rb")

for read in bam.fetch(until_eof=True):
    if read.is_unmapped:
        continue
    rpos = read.reference_start
    qpos = 0
    for op, length in read.cigartuples:
        if op in (4, 5):
            qpos += length
            continue
        break

    print(f"=== read {read.query_name} flag={read.flag} ref_start={read.reference_start} "
          f"query_alignment_start(orig-coord)={qpos} ===")

    found_any = False
    for op, length in read.cigartuples:
        if op in (4, 5):
            continue
        if op in (0, 7, 8):
            qpos += length
            rpos += length
        elif op == 1:
            if q_start <= qpos <= q_end and length >= min_indel:
                print(f"  INSERTION: query {qpos}-{qpos+length} (len={length}) at reference pos {rpos}")
                found_any = True
            qpos += length
        elif op == 2:
            if q_start <= qpos <= q_end and length >= min_indel:
                print(f"  DELETION: at query pos {qpos}, reference {rpos}-{rpos+length} (len={length})")
                found_any = True
            rpos += length
        elif op == 3:
            rpos += length
    if not found_any:
        print(f"  (no indels >= {min_indel}bp found in query window {q_start}-{q_end} for this record)")
