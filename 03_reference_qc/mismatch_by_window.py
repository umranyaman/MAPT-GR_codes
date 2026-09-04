#!/usr/bin/env python3
"""Compute per-window mismatch fraction from samtools mpileup output on stdin.

Parses the base string from each mpileup line: '.' and ',' are matches to the
reference (forward/reverse strand); any ACGTN letter is a mismatch. Indels are
skipped (they shift the count but do not contribute a base to the mismatch
tally). Output is a TSV with one row per non-overlapping window.

Usage:
  samtools mpileup -f ref.fa -r CONTIG bam.bam | \\
      python3 mismatch_by_window.py [window_size] > out.tsv

  window_size: bin size in reference bases (default 500).

Output columns: window_start, total_bases, mismatch_bases, mismatch_fraction
"""
import sys
import collections

WINDOW = int(sys.argv[1]) if len(sys.argv) > 1 else 500


def parse_pileup_line(line):
    parts = line.rstrip("\n").split("\t")
    pos = int(parts[1])
    depth = int(parts[3])
    bases = parts[4] if len(parts) > 4 else ""
    i = 0
    n = len(bases)
    matches = 0
    mismatches = 0
    while i < n:
        c = bases[i]
        if c == "^":
            i += 2
            continue
        if c == "$":
            i += 1
            continue
        if c in ".,":
            matches += 1
            i += 1
            continue
        if c in "ACGTNacgtn":
            mismatches += 1
            i += 1
            continue
        if c in "+-":
            j = i + 1
            num = ""
            while j < n and bases[j].isdigit():
                num += bases[j]
                j += 1
            length = int(num) if num else 0
            i = j + length
            continue
        if c == "*":
            i += 1
            continue
        i += 1
    return pos, matches, mismatches


win_mismatch = collections.defaultdict(int)
win_total = collections.defaultdict(int)

for line in sys.stdin:
    if not line.strip():
        continue
    pos, matches, mismatches = parse_pileup_line(line)
    w = (pos // WINDOW) * WINDOW
    win_mismatch[w] += mismatches
    win_total[w] += matches + mismatches

print("window_start\ttotal_bases\tmismatch_bases\tmismatch_fraction")
for w in sorted(win_total.keys()):
    t = win_total[w]
    m = win_mismatch[w]
    frac = m / t if t > 0 else 0.0
    print(f"{w}\t{t}\t{m}\t{frac:.4f}")
