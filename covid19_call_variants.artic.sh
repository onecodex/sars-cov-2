#!/bin/bash

# To run this script:
# ./covid19_call_variants.artic.sh <samplename>

# This pipeline needs the fastq file of "pass" reads

set -e

#### Default parameters

: "${medaka_model:=r941_min_high_g360}" # Default Medaka model is for MinION (or GridION) R9.4.1 flowcells using the fast Guppy basecaller version 3.6.0, high accuracy base calling
: "${min_read_length:=400}"
: "${max_read_length:=600}" # Filters out all reads above 600 bp
: "${min_read_quality:=7}" # Filters out all reads below a quality of 7
: "${normalized_coverage:=200}" # Normalize to coverage of 200x to reduce runtime
: "${threads:=6}" # Number of threads for running the minimap2/medaka pipeline
# shellcheck disable=SC2154
: "${prefix=report}"

#### 1. Length and quality filtering

echo "[1] trimming/filtering reads with seqkit"
# shellcheck disable=SC1091
source activate report

# shellcheck disable=SC2002
cat "${1}" \
  | seqkit \
    seq \
      --min-len ${min_read_length} \
      --max-len ${max_read_length} \
      --min-qual ${min_read_quality} \
      --out-file "${prefix}.filtered.fastq"

# shellcheck disable=SC2086
echo "total reads after filtering: $(wc -l ${prefix}.filtered.fastq)"

# shellcheck disable=SC1091
source deactivate

#### 2. Alignment, variant calling, and consensus creation
# minimap2 alignment, Medaka for consensus creation and variant calling (experimental)

echo "[2] running ARTIC pipeline"
conda run -n artic \
  artic minion \
  --medaka \
  --medaka-model "${medaka_model}" \
  --normalise "${normalized_coverage}" \
  --threads "${threads}" \
  --scheme-directory ./artic-ncov2019/primer_schemes \
  --read-file "${prefix}.filtered.fastq" \
  --no-longshot \
  --strict nCoV-2019/V3 \
  "${prefix}"

echo "[3] generating variants.tsv"

#### 3. Variant annotation
gunzip "${prefix}.pass.vcf.gz"

# Replace the vcf chromosome name to the verison that snpEff will recognize
# (GenBank sequence NC_045512.2, which is identical to RefSeq MN908947.3)
echo "[3] cleaning up"
### Move files around and clean up
mv "${prefix}.pass.vcf" variants.vcf
mv "${prefix}.consensus.fasta" consensus.fa
mv "${prefix}.sorted.bam" covid19.bam
mv "${prefix}.sorted.bam.bai" covid19.bam.bai
rm -rf "${prefix}.*"

##### 6. Run generate_tsv.py

echo "[ ] finished!"
