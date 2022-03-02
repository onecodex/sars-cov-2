#!/bin/bash

# Example command: `bash covid19_call_variants.sh
# reference/nCoV-2019.reference.fasta
# data/twist-target-capture/RNA_control_spike_in_10_6_100k_reads.fastq.gz
# reference/artic-v1/ARTIC-V1.bed`

# shellcheck disable=SC1091
source activate jobscript-env # noqa

set -eou pipefail

# argv
: "${reference:=${1}}"
: "${input_fastq:=${2}}"
: "${primer_bed_file:=${3}}"

# defaults

: "${threads:=4}"
: "${prefix:=results}"
: "${min_quality:=20}"

echo "reference=${reference}"
echo "input_fastq=${input_fastq}"
echo "primer_bed_file=${primer_bed_file}"
echo "threads=${threads}"
echo "min_quality=${min_quality}"

# minimap2

# Trim polyA tail for alignment (33 bases)
seqtk trimfq -e 33 "${reference}" > trimmed-reference.fasta

# shellcheck disable=SC2086
echo "[1] mapping reads with minimap2"
minimap2 \
  -K 20M \
  -a \
  -x sr \
  -t "${threads}" \
  trimmed-reference.fasta \
  "${input_fastq}" \
  | samtools \
    view \
    -u \
    -h \
    -q $min_quality \
    -F \
    4 - \
  | samtools \
    sort \
    --threads "${threads}" \
    - \
  > "${prefix}.sorted.bam"

# Trim with ivar
echo "[2] Trimming with ivar"
samtools index "${prefix}.sorted.bam"

# samtools flagstat "${prefix}.sorted.bam"
ivar \
  trim \
  -e \
  -q 0 \
  -i "${prefix}.sorted.bam" \
  -b "${primer_bed_file}" \
  -p "${prefix}.ivar"

echo "[3] Sorting and indexing trimmed BAM"
samtools \
  sort \
  -@ "${threads}" \
  "${prefix}.ivar.bam" \
  > "${prefix}.sorted.bam"

samtools \
  index \
  "${prefix}.sorted.bam"

echo "[4] Generating pileup"
bcftools \
  mpileup \
  --annotate FORMAT/AD,INFO/AD \
  --fasta-ref "${reference}" \
  --max-depth 200 \
  --no-BAQ \
  "${prefix}.sorted.bam" \
  | bcftools call \
    --variants-only \
    --multiallelic-caller \
    --output-type z \
    --output "${prefix}.raw.vcf.gz"

# save raw vcf
bcftools view \
  < "${prefix}.raw.vcf.gz" \
  > "${prefix}.raw.vcf"

# filter out low-quality variants
bcftools view \
  --exclude "QUAL<150" \
  --output-type z \
  < "${prefix}.raw.vcf.gz" \
  > "${prefix}.vcf.gz"

# bcftools index requires a .vcf.gz file
# in the special indexed gzip format (can't use regular gzip)
bcftools index "${prefix}.vcf.gz"

# We want to generate a consensus sequence with:
# 1. well-supported deletions marked with '-'
# 2. low-coverage regions masked with N

# vcf is 1-based while bedgraph is 0-based; thus the somewhat convoluted path to mask.bed
bcftools query -f'%CHROM\t%POS0\t%END\n' ${prefix}.vcf.gz > variants.bed
# get low coverage sites (<10) in bedgraph format
bedtools genomecov -bga -ibam ${prefix}.sorted.bam | awk '$4 < 10' > low_coverage_sites.bed
bedtools subtract -a low_coverage_sites.bed -b variants.bed > mask.bed
# generate consensus and rename header
bcftools consensus \
  --fasta-ref "${reference}" \
  --mark-del '-' \
  -m mask.bed \
  "${prefix}.vcf.gz" \
  | sed \
    '/>/ s/$/ | One Codex consensus sequence/' \
    > "${prefix}.consensus.fa"

zcat "${prefix}.vcf.gz" > "${prefix}.vcf"

# Move some files around and clean up
mv "${prefix}.consensus.fa" "consensus.fa"
mv "${prefix}.vcf" "variants.vcf"
mv "${prefix}.sorted.bam" "covid19.bam"
mv "${prefix}.sorted.bam.bai" "covid19.bam.bai"
mv "${prefix}.raw.vcf" "variants.raw.vcf"
rm "${prefix}"*

echo "[ ] Finished!"
