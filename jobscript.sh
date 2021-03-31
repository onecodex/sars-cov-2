#!/bin/bash

set -euo pipefail

: "${THREADS=4}"
: "${ONE_CODEX_REPORT_FILENAME='report.pdf'}"

sample_filename="${1}"
instrument_vendor="Oxford Nanopore"

# generates the following files:
# variants.vcf
# covid19.bam (sorted+bai)
# consensus.fa
if [ "${instrument_vendor}" == "Oxford Nanopore" ]; then
  covid19_call_variants.artic.sh "${sample_filename}"
else
  covid19_call_variants.sh "${sample_filename}"
fi

# needed by report
echo "Getting depth using samtools"
conda run -n report \
  samtools depth covid19.bam > snps.depth

# Count total mapped reads (can we get this from summing snps.depth)
conda run -n report samtools view -F 2308 covid19.bam | wc -l > total_mapped_reads.txt

# call strains, generates:
# nextclade.tsv
# nextclade.json
# pangolin.csv
post_process_variants.sh consensus.fa

# render notebook

echo "Generating notebook!"

RESULTS_DIR="$(pwd)" \
  conda run -n report jupyter \
      nbconvert \
      --execute \
      --to onecodex_pdf \
      --ExecutePreprocessor.timeout=-1 \
      --output="${ONE_CODEX_REPORT_FILENAME}" \
      --output-dir="." \
      /repo/report.ipynb

echo "Finished!"
