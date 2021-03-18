#1/bin/bash

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

# call strains

# Assign NextClade clade
echo "Assigning NextClade Clade"
nextclade --input-fasta consensus.fa --output-tsv nextclade.tsv --output-json nextclade.json

# Assign Pango lineage
echo "Updating Pango Database"
conda run -n pangolin pangolin --update # being intensively updated

# TODO: copy pangolin database data somewhere.
echo "Assigning Pango Lineage"
conda run -n pangolin pangolin consensus.fa --outfile pangolin.csv -t "${THREADS}"

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
