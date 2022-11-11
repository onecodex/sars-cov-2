#!/bin/bash

set -euo pipefail

sample_filename="${1}"

: "${INSTRUMENT_VENDOR:=Illumina}"
: "${ONE_CODEX_REPORT_FILENAME:=report.pdf}"
: "${ARTIC_PRIMER_VERSION:=4.1}"
: "${CONSENSUS_MASK_DEPTH:=10}"

PRIMER_BEDFILE="/primer_schemes/nCoV-2019/V${ARTIC_PRIMER_VERSION}/nCoV-2019.scheme.bed"
INSERT_BEDFILE="/primer_schemes/nCoV-2019/V${ARTIC_PRIMER_VERSION}/nCoV-2019.insert.bed"

echo "--- sample_filename=${sample_filename}"
echo "--- INSTRUMENT_VENDOR=${INSTRUMENT_VENDOR}"

# generates the following files:
# variants.vcf
# covid19.bam (sorted+bai)
# consensus.fa
if [ "${INSTRUMENT_VENDOR}" == "Oxford Nanopore" ]; then
  covid19_call_variants.ont.sh \
	  "${sample_filename}" \
	  "${ARTIC_PRIMER_VERSION}"
else
  covid19_call_variants.sh \
    /reference/nCoV-2019.reference.fasta \
    "${sample_filename}" \
    "${PRIMER_BEDFILE}" \
    "${CONSENSUS_MASK_DEPTH}"
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
conda run -n jobscript-env bcftools csq --force --phase a -f /reference/NC_045512.2.reference.fasta -g /reference/data/NC_045512.2/genes.gff -Ov variants.snpeff.vcf -o variants.snpeff.csq.vcf

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
conda run -n jobscript-env \
  samtools depth -a covid19.bam > snps.depth

# Count total mapped reads
conda run -n jobscript-env samtools view -F 2308 covid19.bam | wc -l > total_mapped_reads.txt

# Generate per-insert depth stats and boxplot
conda run -n jobscript-env insert_coverage_stats.py ${INSERT_BEDFILE} covid19.bam

# call strains

# Assign NextClade clade
echo "Assigning NextClade Clade"

#adding new arg for dataset 'sars-cov-2', 'run' and removing --input-fasta since no longer supported
/usr/local/bin/nextclade run    --input-dataset /usr/local/bin/data/sars-cov-2    --output-tsv nextclade.tsv --output-json nextclade.json   consensus.fa

# Assign Pango lineage

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

#shellcheck disable=SC1000-SC9999
RESULTS_DIR="$(pwd)" SAMPLE_PATH="${sample_filename}" PYTHONWARNINGS="ignore" GIT_DIR="/.git" GIT_WORK_TREE="/" INSTRUMENT_VENDOR="${INSTRUMENT_VENDOR}" ONE_CODEX_REPORT_FILENAME="${ONE_CODEX_REPORT_FILENAME}" ARTIC_PRIMER_VERSION="${ARTIC_PRIMER_VERSION}" conda run -n jobscript-env jupyter nbconvert --execute --to onecodex_pdf --ExecutePreprocessor.timeout=-1 --output="${ONE_CODEX_REPORT_FILENAME}" --output-dir="." report.ipynb

echo "Removing unnecessary files"
rm -f aa_codes.txt \
	annot_table.orfs.txt \
	low_complexity_regions.txt \
	low_coverage_sites.bed \
	mask.bed \
	report.ipynb \
	total_mapped_reads.txt \
	trimmed-reference.fasta \
	variants.bed \
	variants.raw.vcf \
	variants.snpeff.csq.vcf.gz.tbi \
	variants.snpeff.tsv \
	variants.snpeff.vcf \
	variants.snpeff.vcf.faa \
	variants.snpeff.vcf.stats \
	variants.snpeff.vcf.stats.genes.txt \
	variants.vcf

mv variants.snpeff.csq.vcf variants.vcf

echo "Finished!"
