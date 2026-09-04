#!/usr/bin/env bash
#SBATCH --job-name=rebuild_merged_by_sample
#SBATCH --output=logs/rebuild_merged_%j.out
#SBATCH --error=logs/rebuild_merged_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --partition=debug

# Merges raw bam_pass chunks from two PromethION runs (run1 + topup) into one
# unaligned modBAM per sample, preserving 5mCG base modification tags.
#
# Input:  per-barcode bam_pass directories from both runs
# Output: merged_by_sample/<sample>/<sample>_merged.bam  (unaligned modBAM)
#         merged_by_sample/rebuild_manifest.tsv           (provenance record)
#
# H1.0 used three barcodes (bc01 in run1 + bc25 + bc26 in topup) because it
# underperformed in run1 due to barcode dropout; all other samples used one
# barcode per run.
#
# Requires: samtools >=1.22 (micromamba modbamtools env)

#set -euo pipefail

cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")}"
mkdir -p logs

echo "Started: $(date)"

source ~/.bashrc || true
eval "$(micromamba shell hook --shell=bash)"
micromamba activate modbamtools

RUN1_BAMPASS="/home/DSUmranYaman/MAPTGRLRS25/TimBirkle-ONT-WGS-run/251125_TB_NBD_run1/TB_Pool_1-4/20251125_1533_3F_PBI60991_5103ca6f/bam_pass"
TOPUP_BAMPASS="/home/DSUmranYaman/MAPTGRLRS25/TimBirkle-ONT-WGS-run/260317_TB_NBD_topupseq_run1/TB_pool1_top-upseq/20260317_1301_3H_PBK63983_3405bcb9/bam_pass"
OUTROOT="/home/DSUmranYaman/MAPTGRLRS25/TimBirkle-ONT-WGS-run/merged_by_sample"

mkdir -p "$OUTROOT"

# Sample -> contributing (run_dir barcode) pairs, per the LRS facility's
# barcode sheets (run1 email 2025-12-18; topup email 2026-xx from Hermione).
# H1.0 was deliberately re-sequenced twice in the topup (bc25 + bc26) because
# it underperformed in run1 (barcode01 dropout).
declare -A SAMPLE_SOURCES
SAMPLE_SOURCES["H1.0"]="$RUN1_BAMPASS/barcode01 $TOPUP_BAMPASS/barcode25 $TOPUP_BAMPASS/barcode26"
SAMPLE_SOURCES["H2.1"]="$RUN1_BAMPASS/barcode02 $TOPUP_BAMPASS/barcode27"
SAMPLE_SOURCES["IVS10+16"]="$RUN1_BAMPASS/barcode03 $TOPUP_BAMPASS/barcode28"
SAMPLE_SOURCES["N279K"]="$RUN1_BAMPASS/barcode04 $TOPUP_BAMPASS/barcode29"

MANIFEST="$OUTROOT/rebuild_manifest.tsv"
echo -e "sample\tsource_dir\tn_chunk_files\ttotal_bytes" > "$MANIFEST"

for sample in "${!SAMPLE_SOURCES[@]}"; do
    echo "=============================="
    echo "Sample: $sample"

    out_dir="$OUTROOT/$sample"
    mkdir -p "$out_dir"
    out_bam="$out_dir/${sample}_merged.bam"

    # Gather every raw chunk BAM from every contributing barcode dir,
    # verifying each source dir exists and each chunk file is intact
    # before we touch anything.
    chunk_files=()
    for src_dir in ${SAMPLE_SOURCES[$sample]}; do
        [[ -d "$src_dir" ]] || { echo "FATAL: missing source dir $src_dir"; exit 1; }

        n=0
        bytes=0
        for f in "$src_dir"/*.bam; do
            [[ -f "$f" ]] || continue
            echo "  checking $f"
            samtools quickcheck -v "$f" 2>/dev/null || true   # expected "no targets" noise on uBAM, ignore
            if ! samtools flagstat "$f" > /dev/null 2>"/tmp/rebuild_err_$$"; then
                echo "FATAL: $f failed a full read (truncated/corrupt raw chunk)."
                cat "/tmp/rebuild_err_$$"
                rm -f "/tmp/rebuild_err_$$"
                exit 1
            fi
            rm -f "/tmp/rebuild_err_$$"
            chunk_files+=("$f")
            n=$((n+1))
            bytes=$((bytes + $(stat -c%s "$f")))
        done
        echo -e "$sample\t$src_dir\t$n\t$bytes" >> "$MANIFEST"
    done

    echo "  Concatenating ${#chunk_files[@]} chunk files -> $out_bam"
    samtools cat -o "$out_bam" "${chunk_files[@]}"

    echo "  Verifying merged output..."
    samtools quickcheck -v "$out_bam" 2>/dev/null || true   # same uBAM false-positive as above
    if ! samtools flagstat "$out_bam" > "$out_dir/${sample}_merged.flagstat.txt"; then
        echo "FATAL: merged output for $sample failed a full read — rebuild did not produce a clean file."
        exit 1
    fi

    echo "  OK: $out_bam"
    cat "$out_dir/${sample}_merged.flagstat.txt"

    chunk_files=()
done

echo ""
echo "All samples rebuilt successfully. Manifest: $MANIFEST"
echo "Finished: $(date)"
