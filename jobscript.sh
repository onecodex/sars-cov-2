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

echo "Annotating VCF file using snpEff"
# prepare for vcf annotation with snpEff
if [ ! -d "/usr/local/bin/snpEff/data" ]; then
  mkdir /usr/local/bin/snpEff/data
  mkdir /usr/local/bin/snpEff/data/NC_045512.2
  mv /nCoV-2019.reference.fasta /usr/local/bin/snpEff/data/NC_045512.2/sequences.fa
  mv /nCoV-2019.reference.gbk /usr/local/bin/snpEff/data/NC_045512.2/genes.gbk
fi

# match chromosome name between the genbank file, fasta files, and vcf
# MN908947.3 and NC_045512.2 are identical genomes (Refseq vs assembly numbers)
sed -i "s|MN908947.3|NC_045512.2|" /usr/local/bin/snpEff/data/NC_045512.2/sequences.fa
sed -i "s|MN908947.3|NC_045512.2|" variants.vcf

# build custom snpeff database
echo "NC_045512.2.genome : nCoV-2019 ARTIC V3" >> /usr/local/bin/snpEff/snpEffect.config
java -Xmx4g -jar /usr/local/bin/snpEff/snpEff.jar build -c /usr/local/bin/snpEff/snpEffect.config -noGenome -genbank -v NC_045512.2
# run snpeff annotation on vcf
java -Xmx4g -jar /usr/local/bin/snpEff/snpEff.jar ann NC_045512.2 -verbose -config /usr/local/bin/snpEff/snpEffect.config -fastaProt variants.snpeff.vcf.faa -csvStats variants.snpeff.vcf.stats variants.vcf > variants.snpeff.vcf
# if Illumina, extract fields of interest from annotated vcf into a tsv
if [ "${INSTRUMENT_VENDOR}" == "Illumina" ]; then
  java -Xmx4g -jar /usr/local/bin/snpEff/SnpSift.jar extractFields variants.snpeff.vcf POS REF ALT DP4 ANN[0].EFFECT ANN[0].HGVS_P > variants.snpeff.vcf.extracted.tsv
fi
# edit the effects field
sed -i 's/_/ /g' variants.snpeff.vcf.extracted.tsv
sed -i 's/&/; /g' variants.snpeff.vcf.extracted.tsv

mv variants.snpeff.vcf.extracted.tsv variants.snpeff.tsv

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
