#!/usr/bin/env bash
#SBATCH --job-name=ont_dna_align_modkit
#SBATCH --output=logs/ont_dna_align_modkit_%j.out
#SBATCH --error=logs/ont_dna_align_modkit_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=50
#SBATCH --mem=100G
#SBATCH --time=48:00:00
#SBATCH --partition=debug

# Aligns per-sample unaligned modBAMs to the MAPT-GR hybrid reference and
# produces CpG methylation pileups via modkit.
#
# Input:  merged_by_sample/<sample>/<sample>_merged.bam  (from 01_merge)
# Output per sample (under dna_align_modkit_v3/<sample>/):
#   aligned_sorted_<sample>.bam / .bai  -- coordinate-sorted aligned modBAM
#   aln_mod_<sample>.bed.gz / .tbi     -- modkit pileup (traditional preset)
#   aligned_sorted_<sample>.flagstat.txt
#   aligned_sorted_<sample>.idxstats.txt
#
# Reference: GRCm39 + chrTau (NC_000017.11, 190,497 bp MAPT-GR insert) with
#            the replaced mouse interval (chr11:104,066,949-104,223,491) masked.
# Aligner:   minimap2 v2.30  -ax map-ont  -y (pass MM/ML tags through)
#
# Requires: samtools >=1.22, minimap2 >=2.30, modkit, bgzip, tabix
#           (micromamba modbamtools env)

set -euo pipefail

cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")}"
WORKDIR=$(pwd)
mkdir -p logs

echo "Started: $(date)"

source ~/.bashrc || true
eval "$(micromamba shell hook --shell=bash)"
micromamba activate modbamtools

THREADS=50

REF_FA="/home/DSUmranYaman/Tim/D_genome/umran/Updated_hybrid_genome/build/new_hybrid_genome_extended_chrTau.maskfix.fa"
MERGED_ROOT="/home/DSUmranYaman/MAPTGRLRS25/TimBirkle-ONT-WGS-run/merged_by_sample"
OUTDIR="$WORKDIR/dna_align_modkit_v3"

mkdir -p "$OUTDIR"

echo "REF_FA      : $REF_FA"
echo "MERGED_ROOT : $MERGED_ROOT"
echo "OUTDIR      : $OUTDIR"

[[ -f "$REF_FA" ]] || { echo "Reference not found: $REF_FA"; exit 1; }
[[ -d "$MERGED_ROOT" ]] || { echo "Input directory not found: $MERGED_ROOT"; exit 1; }

samtools faidx "$REF_FA"

for sdir in "$MERGED_ROOT"/*/; do
    [[ -d "$sdir" ]] || continue

    sample=$(basename "$sdir")
    echo "=============================="
    echo "Processing $sample"

    input_bam="$sdir/${sample}_merged.bam"

    if [[ ! -f "$input_bam" ]]; then
        echo "Skipping $sample: no ${sample}_merged.bam found"
        continue
    fi

    echo "Input BAM: $input_bam"

    echo "Verifying source BAM integrity before touching it..."
    samtools quickcheck -v "$input_bam" 2>/dev/null || true   # "no targets" false positive on uBAM, ignore
    if ! samtools flagstat "$input_bam" > /dev/null; then
        echo "FATAL: $input_bam failed a full read (truncated/corrupt). Skipping $sample."
        continue
    fi
    echo "Source BAM OK."

    out_sd="$OUTDIR/$sample"
    mkdir -p "$out_sd"
    cd "$out_sd"

    merged_sorted="merged_sorted_${sample}.bam"
    aligned_bam="aligned_${sample}.mod.bam"
    aligned_sorted="aligned_sorted_${sample}.bam"
    pileup_bed="aln_mod_${sample}.bed"

    echo "Sorting merged BAM..."
    samtools sort -@ "$THREADS" "$input_bam" > "$merged_sorted"
    samtools index -@ "$THREADS" "$merged_sorted"

    echo "Converting to FASTQ and aligning..."
    samtools fastq -@ "$THREADS" -T MM,ML "$merged_sorted" \
        | minimap2 -ax map-ont -t "$THREADS" -y "$REF_FA" - \
        > "$aligned_bam"

    echo "Sorting aligned BAM..."
    samtools sort -@ "$THREADS" "$aligned_bam" > "$aligned_sorted"
    samtools index -@ "$THREADS" "$aligned_sorted"

    echo "Running modkit pileup..."
    modkit pileup "$aligned_sorted" "$pileup_bed" \
        --ref "$REF_FA" \
        --threads "$THREADS" \
        --preset traditional

    echo "Compressing and indexing BED..."
    bgzip -f "$pileup_bed"
    tabix -f -p bed "${pileup_bed}.gz"

    samtools flagstat -@ "$THREADS" "$aligned_sorted" > "aligned_sorted_${sample}.flagstat.txt"
    samtools idxstats "$aligned_sorted" > "aligned_sorted_${sample}.idxstats.txt"

    echo "Done: $out_sd/${pileup_bed}.gz"
    cd "$WORKDIR"
done

echo "Finished: $(date)"
