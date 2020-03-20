#!/bin/bash

# Example command: `bash covid19_call_variants.sh reference/nCoV-2019.reference.fasta data/twist-target-capture/RNA_control_spike_in_10_6_100k_reads.fastq.gz reference/artic-v1/ARTIC-V1.bed`
set -eou pipefail

reference=$1
input_fastq=$2
bed_file=$3
threshold=600
threads=$(grep -c processor /proc/cpuinfo)
prefix="results"

mapping_mode=$(head -n 400 "$input_fastq" | awk -v "thresh=$threshold" 'NR % 4 == 2 {s+= length}END {if (s/(NR/4) > thresh) {print "map-ont"} else {print "sr"}}')

# Minimap2
MINIMAP_OPTS="-K 20M -a -x $mapping_mode -t $threads"
echo "[1] Running minimap with options: $MINIMAP_OPTS"
minimap2 $MINIMAP_OPTS "$reference" "$input_fastq"  | samtools view -u -h -F 4 - | samtools sort -@ "$threads" - > "${prefix}.sorted.bam"

# Get only mapped reads
echo "[2] Removing unmapped reads"

# Trim with ivar
echo "[2] Trimming with ivar"
samtools index "${prefix}.sorted.bam"
# samtools flagstat "${prefix}.sorted.bam"
ivar trim -e -i "${prefix}.sorted.bam" -b "$bed_file" -p "${prefix}.ivar"

# Trim with fastp
echo "[3] Sorting and indexing trimmed BAM"
samtools sort -@ "$threads" "${prefix}.ivar.bam" > "${prefix}.sorted.bam"
samtools index "${prefix}.sorted.bam"
# samtools flagstat "${prefix}.sorted.bam"

# Call variants with ivar
# TODO: Investigate extra variants for Twist capture data:
# MN908947.3	29750	C	T	0	0	0	4	3	33	0.666667	6	4.93577e-08	TRUE
# MN908947.3	29754	C	G	0	0	0	2	0	32	1	2	9.57854e-05	TRUE
# MN908947.3	29759	G	A	0	0	0	57	51	35	1	57	1.86298e-58	TRUE
echo "[4] Generating variants TSV"
samtools mpileup -f "$reference" -d 1000000 -A -B -Q 0 "${prefix}.sorted.bam" | ivar variants -p "${prefix}.ivar" -t 0.6

# Generate consensus sequence with ivar
echo "[5] Generating consensus sequence"
samtools mpileup -f "$reference" -d 1000000 -A -B -Q 0 "${prefix}.sorted.bam" | ivar consensus -p "${prefix}.ivar" -m 1 -t 0.6 -n N

# Remove excess files
# TODO: This was for Snippy output. Rename to include reference accession and ivar params
sed '/>/ s/$/ | One Codex consensus sequence/' < "${prefix}.ivar.fa" > consensus.fa
mv "${prefix}.ivar.tsv" "variants.tsv"
rm "${prefix}.ivar.bam" "${prefix}.ivar.qual.txt" "${prefix}.sorted.bam" "${prefix}.sorted.bam.bai" "${prefix}.ivar.fa"
