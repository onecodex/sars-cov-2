#!/bin/bash

# Example command: `bash covid19_call_variants.sh
# reference/nCoV-2019.reference.fasta
# data/twist-target-capture/RNA_control_spike_in_10_6_100k_reads.fastq.gz
# reference/artic-v1/ARTIC-V1.bed`

# shellcheck disable=SC1091
source activate report # noqa

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
  --max-depth 0 \
  --count-orphans \
  --no-BAQ \
  --min-BQ 0 \
  "${prefix}.sorted.bam" \
  | bcftools call \
    --variants-only \
    --multiallelic-caller \
    --output-type z \
    --output "${prefix}.raw.vcf.gz"


# filter out low-quality variants
bcftools view \
  --exclude "QUAL<150" \
  --output-type z \
  < "${prefix}.raw.vcf.gz" \
  > "${prefix}.vcf.gz"

# bcftools index requires a .vcf.gz file
# in the special indexed gzip format (can't use regular gzip)
bcftools index "${prefix}.vcf.gz"

bcftools consensus \
  --fasta-ref "${reference}" \
  "${prefix}.vcf.gz" \
  | sed \
    '/>/ s/$/ | One Codex consensus sequence/' \
    > "${prefix}.consensus.fa"

zcat "${prefix}.vcf.gz" > "${prefix}.vcf"

# prepare for vcf annotation with snpEff
ls -lash /usr/local/bin/snpEff
mkdir /usr/local/bin/snpEff/data
mkdir /usr/local/bin/snpEff/data/NC_045512.2
mv /nCoV-2019.reference.fasta /usr/local/bin/snpEff/data/NC_045512.2/sequences.fa
mv /nCoV-2019.reference.gbk /usr/local/bin/snpEff/data/NC_045512.2/genes.gbk

# match chromosome name between the genbank file, fasta files, and vcf
# MN908947.3 and NC_045512.2 are identical genomes (Refseq vs assembly numbers)
sed -i "s|MN908947.3|NC_045512.2|" /usr/local/bin/snpEff/data/NC_045512.2/sequences.fa
sed -i "s|MN908947.3|NC_045512.2|" ${prefix}.vcf

# build custom snpeff database
echo "NC_045512.2.genome : nCoV-2019 ARTIC V3" >> /usr/local/bin/snpEff/snpEffect.config
java -Xmx4g -jar /usr/local/bin/snpEff/snpEff.jar build -c /usr/local/bin/snpEff/snpEffect.config -noGenome -genbank -v NC_045512.2
# run snpeff annotation on vcf
java -Xmx4g -jar /usr/local/bin/snpEff/snpEff.jar ann NC_045512.2 -verbose -config /usr/local/bin/snpEff/snpEffect.config -fastaProt ${prefix}.snpeff.vcf.faa -csvStats ${prefix}.snpeff.vcf.stats ${prefix}.vcf > ${prefix}.snpeff.vcf
# extract fields of interest from annotated vcf
java -Xmx4g -jar /usr/local/bin/snpEff/SnpSift.jar extractFields ${prefix}.snpeff.vcf POS REF ALT ANN[0].GENE ANN[0].GENEID ANN[0].EFFECT ANN[0].HGVS_P > ${prefix}.snpeff.vcf.extracted.tsv

# Move some files around and clean up
mv "${prefix}.consensus.fa" "consensus.fa"
mv "${prefix}.vcf" "variants.vcf"
mv "${prefix}.sorted.bam" "covid19.bam"
mv "${prefix}.sorted.bam.bai" "covid19.bam.bai"
mv "${prefix}.snpeff.vcf.extracted.tsv" "variants.snpeff.tsv"
rm "${prefix}"*

#echo "[ ] Finished!"
