#!/usr/bin/env python3
"""Find per-read H1->H2 switch breakpoints for reads spanning a candidate region.

Usage: python3 find_breakpoint.py <bam> <chrom> <region_start> <region_end> [span_margin]

Requires pysam. For each primary-aligned read whose alignment fully spans
[region_start - span_margin, region_end + span_margin] on the reference,
walk the aligned pairs (with_seq=True, so mismatches are lowercase in the
returned query sequence per pysam convention) and find the reference
position that maximizes (mismatch_rate_after - mismatch_rate_before) --
i.e. the best single change-point for that read.
"""
import sys
import pysam

bam_path = sys.argv[1]
ref_name = sys.argv[2]
region_start = int(sys.argv[3])
region_end = int(sys.argv[4])
margin = int(sys.argv[5]) if len(sys.argv) > 5 else 3000

bam = pysam.AlignmentFile(bam_path, "rb")

results = []

for read in bam.fetch(ref_name, region_start, region_end):
    if read.is_secondary or read.is_supplementary or read.is_unmapped:
        continue
    if read.reference_start > region_start - margin or read.reference_end is None or read.reference_end < region_end + margin:
        continue  # doesn't fully span with margin

    pairs = read.get_aligned_pairs(with_seq=True)
    # build a per-reference-position mismatch indicator, restricted to the window we care about
    positions = []
    mismatches = []
    for qpos, rpos, refbase in pairs:
        if rpos is None or qpos is None:
            continue  # skip indels for this scan
        if rpos < region_start - margin or rpos > region_end + margin:
            continue
        is_mismatch = refbase is not None and refbase.islower()  # pysam: lowercase ref base = mismatch
        positions.append(rpos)
        mismatches.append(1 if is_mismatch else 0)

    if len(positions) < 500:
        continue

    n = len(positions)
    prefix = [0] * (n + 1)
    for i in range(n):
        prefix[i + 1] = prefix[i] + mismatches[i]

    best_c = None
    best_score = -1
    # grid search candidate breakpoints (every 10th position) within the region of interest
    for i in range(1, n - 1):
        rpos = positions[i]
        if rpos < region_start or rpos > region_end:
            continue
        before_n = i
        after_n = n - i
        if before_n < 100 or after_n < 100:
            continue
        rate_before = prefix[i] / before_n
        rate_after = (prefix[n] - prefix[i]) / after_n
        score = rate_after - rate_before
        if score > best_score:
            best_score = score
            best_c = rpos

    if best_c is not None:
        results.append((read.query_name, best_c, best_score, read.reference_start, read.reference_end))

print(f"# {len(results)} reads spanning [{region_start - margin}, {region_end + margin}] with margin {margin}")
print("read_name\tbreakpoint_pos\tstep_score\tread_ref_start\tread_ref_end")
for r in sorted(results, key=lambda x: x[1]):
    print(f"{r[0]}\t{r[1]}\t{r[2]:.4f}\t{r[3]}\t{r[4]}")

if results:
    bps = [r[1] for r in results]
    bps.sort()
    median = bps[len(bps) // 2]
    print(f"\n# median breakpoint: {median}", file=sys.stderr)
    print(f"# range: {min(bps)} - {max(bps)}", file=sys.stderr)
