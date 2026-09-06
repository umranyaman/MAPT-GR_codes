# MAPT-GR_codes

Analysis scripts for long-read whole-genome sequencing of MAPT Genomic Replacement (MAPT-GR) mouse lines, as reported in:

> Birkle TJY, Foiani M, Yaman Ü, Damoc LS, Nirujogi RS, Tsefou E, Geary B, Avdic-Belltheus A, Watamura N, Faber K, Santambrogio A, Patel S, Davies H, Vendruscolo M, Duff KE. Divergent tau phosphorylation across MAPT mutations despite similar 4R tau shifts in MAPT-GR mice. (2026)

These scripts produced the sequencing depth statistics and GR locus coverage data reported in **Table 1**, and the B6J/B6N background and H1/H2 subhaplotype genotype calls in **Table S1**.

## Lines sequenced

Four lines on Oxford Nanopore PromethION 24 (R10.4.1 flow cells), olfactory bulb, 7-month-old adult males:

- H1.0 (JAX #035398) — MAPT H1 haplotype
- H2.1 (JAX #033668) — MAPT H2 haplotype
- IVS10+16 (JAX #036664) — H1 + intronic splice-site mutation (rs63751011)
- N279K (JAX #035794) — H1 + N279K pathogenic mutation

## Reference

Custom hybrid: GRCm39 + 190,497 bp MAPT-GR insert from GRCh38 (NC_000017.11) as `chrTau`, with the replaced mouse interval (chr11:104,066,949–104,223,491) masked.

## Scripts

### 01_merge/
`rebuild_merged_by_sample.slurm` — merges raw bam_pass chunks from run1 and topup sequencing into a single unaligned modBAM per sample. H1.0 required three barcodes due to dropout in run1.

### 02_alignment_methylation/
`alignment_by_sample.slurm` — aligns merged modBAMs to the hybrid reference (minimap2 v2.30, `-ax map-ont`), preserving 5mCG base modification tags, and runs modkit pileup for CpG methylation. The aligned BAMs and flagstat outputs underlie the Table 1 sequencing depth and GR locus coverage statistics.

### 03_reference_qc/
Python scripts for reference validation and genotype QC. Together these scripts validated the GR insert integrity and produced the genotype data in Table S1.

- `check_b6_background.py` — pysam pileup-based B6J/B6N genotype calls at 34 coding-variant positions (Simon et al. 2013); feeds Table S1
- `find_breakpoint.py` — per-read H1/H2 haplotype switch breakpoint detection by mismatch-rate change-point scan; feeds Table S1 subhaplotype assignment
- `mismatch_by_window.py` — windowed mismatch fraction from samtools mpileup; used for reference QC
- `cigar_window.py` — CIGAR-level indel inspection in query coordinates; used for reference QC
- `check_boundary_signature.py` — clip and split-read signatures near insert boundaries; GR locus validation
- `check_chrTau_crosslink.py` — split reads spanning the chrTau–mouse junction; GR locus validation
- `check_flanking_deletion.py` — large CIGAR deletions at the masked gap; GR locus validation
- `check_inversion_signature.py` — opposite-strand SA tag signatures indicating inversion breakpoints; GR locus validation
- `query_to_ref.py` — maps query (read) coordinates to reference positions; coordinate utility

## Environment

micromamba environment (`modbamtools`): samtools ≥1.22, minimap2 ≥2.30, modkit, htslib, pysam

## Data

Raw sequencing data are deposited at NCBI SRA under BioProject PRJNA1524514 (BioSamples SAMN62928088–SAMN62928091). Files are aligned modBAMs mapped to the custom hybrid reference described above.
