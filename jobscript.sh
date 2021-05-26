#!/bin/bash

set -euo pipefail

sample_filename="${1}"

: "${INSTRUMENT_VENDOR:=Illumina}"
: "${ONE_CODEX_REPORT_FILENAME:=report.pdf}"

echo "--- sample_filename=${sample_filename}"
echo "--- INSTRUMENT_VENDOR=${INSTRUMENT_VENDOR}"

########### First, attempt to update PangoLEARN. If update fails, the job fails

numtries=10
for i in $(seq 1 $numtries); do
	sleeptime=$((10+i*5))
	[ "$i" -gt 1 ] && sleep $sleeptime
	conda run -n pangolin pangolin --update && s=0 && break || s=$?
done
{ echo "Failed to update pangoLEARN after $numtries tries" ; exit $s; }


# generates the following files:
# variants.vcf
# covid19.bam (sorted+bai)
# consensus.fa
if [ "${INSTRUMENT_VENDOR}" == "Oxford Nanopore" ]; then
  covid19_call_variants.ont.sh "${sample_filename}"
else
  covid19_call_variants.sh \
    /reference/nCoV-2019.reference.fasta \
    "${sample_filename}" \
    /reference/ARTIC-V3.bed
fi

echo "Annotating VCF file using snpEff"
# match chromosome name between the genbank file, fasta files, and vcf
# MN908947.3 and NC_045512.2 are identical genomes (Refseq vs assembly numbers)
sed -i "s|MN908947.3|NC_045512.2|" variants.vcf

# build custom snpeff database
java -Xmx4g -jar /usr/local/bin/snpEff/snpEff.jar build -c /reference/snpEffect.config -noGenome -gff3 -v NC_045512.2

# run snpeff annotation on vcf
# snpeff expects your sequence.fa and genes.gbk file to be in ../data/NC_045512.2 relative to snpEffect.config
java -Xmx4g -jar /usr/local/bin/snpEff/snpEff.jar ann NC_045512.2 -verbose -config /reference/snpEffect.config -fastaProt variants.snpeff.vcf.faa -csvStats variants.snpeff.vcf.stats variants.vcf > variants.snpeff.vcf

# run bcftools csq to link consecutive SNPs on the same codon (BCSQ field)
conda run -n report bcftools csq --force --phase a -f /reference/NC_045512.2.reference.fasta -g /reference/data/NC_045512.2/genes.gff -Ov variants.snpeff.vcf -o variants.snpeff.csq.vcf

# Extract fields of interest from annotated vcf into a tsv, treating SR and DP4 as essentially identical information
# In the Medaka-generated vcf for ONT data:
# SR="Depth of spanning reads by strand which best align to each allele (ref fwd, ref rev, alt1 fwd, alt1 rev, etc.)"
# In the bcftools-generated vcf for Illumina data:
# DP4="Number of high-quality ref-forward , ref-reverse, alt-forward and alt-reverse bases"
if [ "${INSTRUMENT_VENDOR}" == "Oxford Nanopore" ]; then
  java -Xmx4g -jar /usr/local/bin/snpEff/SnpSift.jar extractFields variants.snpeff.csq.vcf POS REF ALT SR ANN[0].EFFECT ANN[0].HGVS_P BCSQ > variants.snpeff.tsv
  sed -i 's/SR/allele reads by strand/' variants.snpeff.tsv
else
  java -Xmx4g -jar /usr/local/bin/snpEff/SnpSift.jar extractFields variants.snpeff.csq.vcf POS REF ALT DP4 ANN[0].EFFECT ANN[0].HGVS_P BCSQ > variants.snpeff.tsv
  sed -i 's/DP4/allele reads by strand/' variants.snpeff.tsv
fi

# edit the effects field to make it more readable in the table
sed -i 's/_/ /g' variants.snpeff.tsv
sed -i 's/&/; /g' variants.snpeff.tsv

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
cp /reference/annot_table.orfs.txt .
cp /reference/low_complexity_regions.txt .
cp /reference/aa_codes.txt .

echo "Generating notebook!"

RESULTS_DIR="$(pwd)" \
SAMPLE_PATH="${sample_filename}" \
PYTHONWARNINGS="ignore" \
GIT_DIR="/.git" \
GIT_WORK_TREE="/" \
conda run -n report jupyter \
      nbconvert \
      --execute \
      --to onecodex_pdf \
      --ExecutePreprocessor.timeout=-1 \
      --output="${ONE_CODEX_REPORT_FILENAME}" \
      --output-dir="." \
      report.ipynb

echo "Finished!"
