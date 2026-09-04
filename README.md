# MAPT-GR_codes

Analysis scripts for long-read whole-genome sequencing of MAPT Genomic Replacement (MAPT-GR) mouse lines.

## Overview

Four lines sequenced on Oxford Nanopore PromethION 24 (R10.4.1 flow cells):
- H1.0 (JAX #035398) — MAPT H1 haplotype
- H2.1 (JAX #033668) — MAPT H2 haplotype
- IVS10+16 (JAX #036664) — H1 + intronic splice-site mutation (rs63751011)
- N279K (JAX #035794) — H1 + N279K pathogenic mutation

Tissue: olfactory bulb, 7-month-old adult males.

## Reference

Custom hybrid: GRCm39 + 190,497 bp MAPT-GR insert from GRCh38 (NC_000017.11) as `chrTau`, with the replaced mouse interval (chr11:104,066,949–104,223,491) masked.

## Scripts

### 01_merge/
`rebuild_merged_by_sample.slurm` — merges raw bam_pass chunks from run1 and topup sequencing into a single unaligned modBAM per sample.

### 02_alignment_methylation/
`alignment_by_sample.slurm` — aligns merged modBAMs to the hybrid reference (minimap2 v2.30, `-ax map-ont`), preserving 5mCG base modification tags, and runs modkit pileup for CpG methylation analysis.

### 03_reference_qc/
Python scripts for reference validation and genotype QC:
- `check_b6_background.py` — B6J/B6N genotype at coding-variant positions
- `find_breakpoint.py` — H1/H2 haplotype switch breakpoint detection per read
- `mismatch_by_window.py` — windowed mismatch fraction from mpileup
- `cigar_window.py` — CIGAR-level indel inspection
- `check_boundary_signature.py`, `check_chrTau_crosslink.py`, `check_flanking_deletion.py`, `check_inversion_signature.py` — GR insert structural validation
- `query_to_ref.py` — coordinate mapping utility

## Environment

micromamba environment (`modbamtools`): samtools ≥1.22, minimap2 ≥2.30, modkit, htslib, pysam

## Data

Raw sequencing data deposited at NCBI SRA: SUB16460615
