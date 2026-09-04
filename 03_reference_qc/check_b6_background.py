#!/usr/bin/env python3
"""Check B6J vs B6N genotype at each Simon et al. 2013 Table 1 coding-variant
position (mm39-lifted), for one BAM, using pysam pileup base counting.

Reports mean MAPQ of the reads used at each position, and drops reads below
a minimum mapping quality (default MQ>=20) before counting bases, so a single
stray low-MQ misaligned read can't flip a majority-vote call at low-depth sites.

Usage: python3 check_b6_background.py <bam_path> <positions.tsv> [min_mapq]

positions.tsv columns (tab-separated, no header):
chrom  pos(1-based)  gene  B6J_base  B6N_base
"""
import sys
import pysam
from collections import Counter

bam_path = sys.argv[1]
pos_path = sys.argv[2]
min_mapq = int(sys.argv[3]) if len(sys.argv) > 3 else 20

bam = pysam.AlignmentFile(bam_path, "rb")

rows = []
with open(pos_path) as f:
    for line in f:
        chrom, pos, gene, b6j, b6n = line.rstrip("\n").split("\t")
        rows.append((chrom, int(pos), gene, b6j, b6n))

n_b6j = 0
n_b6n = 0
n_other = 0
n_nodata = 0
n_dropped_lowmq = 0

print(f"{'Gene':<12}{'Chrom:Pos':<20}{'Call_base':<11}{'DP_used':<9}{'DP_dropped_lowMQ':<18}{'MeanMQ':<8}{'B6J':<5}{'B6N':<5}{'Call'}")
for chrom, pos, gene, b6j, b6n in rows:
    counts = Counter()
    mapqs = []
    dropped = 0
    for pileupcolumn in bam.pileup(chrom, pos - 1, pos, truncate=True, min_base_quality=0):
        if pileupcolumn.pos != pos - 1:
            continue
        for pileupread in pileupcolumn.pileups:
            if pileupread.is_del or pileupread.is_refskip:
                continue
            mq = pileupread.alignment.mapping_quality
            if mq < min_mapq:
                dropped += 1
                continue
            base = pileupread.alignment.query_sequence[pileupread.query_position]
            counts[base.upper()] += 1
            mapqs.append(mq)

    total = sum(counts.values())
    n_dropped_lowmq += dropped
    mean_mq = f"{sum(mapqs)/len(mapqs):.1f}" if mapqs else "-"

    if total == 0:
        call = "NO_DATA" if dropped == 0 else "NO_DATA(all_reads_lowMQ)"
        top_base = "-"
        n_nodata += 1
    else:
        top_base, top_n = counts.most_common(1)[0]
        if top_base == b6j and top_base != b6n:
            call = "B6J"
            n_b6j += 1
        elif top_base == b6n and top_base != b6j:
            call = "B6N"
            n_b6n += 1
        elif top_base == b6j == b6n:
            call = "AMBIG(same base both strains)"
            n_other += 1
        else:
            call = f"OTHER({top_base} matches neither)"
            n_other += 1
        if top_n < total:
            call += f" [mixed: {dict(counts)}]"

    print(f"{gene:<12}{chrom+':'+str(pos):<20}{top_base:<11}{total:<9}{dropped:<18}{mean_mq:<8}{b6j:<5}{b6n:<5}{call}")

print(f"\n# Summary for {bam_path} (min_mapq={min_mapq})", file=sys.stderr)
print(f"# B6J matches: {n_b6j}, B6N matches: {n_b6n}, other/ambiguous: {n_other}, no usable coverage: {n_nodata}", file=sys.stderr)
print(f"# Total low-MQ reads dropped across all positions: {n_dropped_lowmq}", file=sys.stderr)
