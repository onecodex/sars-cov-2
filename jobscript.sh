#!/bin/bash

set -euo pipefail

sample_filename="${1}"

: "${INSTRUMENT_VENDOR:=Illumina}"
: "${ONE_CODEX_REPORT_FILENAME:=report.pdf}"

echo "--- sample_filename=${sample_filename}"
echo "--- INSTRUMENT_VENDOR=${INSTRUMENT_VENDOR}"

# generates the following files:
# variants.vcf
# covid19.bam (sorted+bai)
# consensus.fa
if [ "${INSTRUMENT_VENDOR}" == "Oxford Nanopore" ]; then
  covid19_call_variants.artic.sh "${sample_filename}"
else
  covid19_call_variants.sh \
    /share/nCoV-2019.reference.fasta \
    "${sample_filename}" \
    /share/ARTIC-V3.bed
fi

# needed by report
echo "Getting depth using samtools"
conda run -n report \
  samtools depth covid19.bam > snps.depth

# Count total mapped reads (can we get this from summing snps.depth)
conda run -n report samtools view -F 2308 covid19.bam | wc -l > total_mapped_reads.txt

# call strains

# Assign NextClade clade
echo "Assigning NextClade Clade"
nextclade --input-fasta consensus.fa --output-tsv nextclade.tsv --output-json nextclade.json

# Assign Pango lineage
echo "Updating Pango Database"
#conda run -n pangolin pangolin --update # being intensively updated

# TODO: copy pangolin database data somewhere.
echo "Assigning Pango Lineage"
conda run -n pangolin pangolin consensus.fa --outfile pangolin.csv

ls -lash /

# render notebook
cp /report.ipynb .
cp /annot_table.orfs.txt .
cp /share/low_complexity_regions.txt .
cp /share/aa_codes.txt .

echo "Generating notebook!"

RESULTS_DIR="$(pwd)" \
SAMPLE_PATH="${sample_filename}" \
PYTHONWARNINGS="ignore" \
conda run -n report jupyter \
      nbconvert \
      --execute \
      --to onecodex_pdf \
      --ExecutePreprocessor.timeout=-1 \
      --output="${ONE_CODEX_REPORT_FILENAME}" \
      --output-dir="." \
      report.ipynb

echo "Finished!"
