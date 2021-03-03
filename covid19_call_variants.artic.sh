#!/bin/bash

# To run this script:
# ./covid19_call_variants.artic.sh <samplename>

# This pipeline needs:
# the fastq file of "pass" reads;
# sequencing_summary.txt


samplename=$1
mkdir ${samplename}
mv *.fastq ${samplename}/


#### Default parameters

: "${medaka_model:=r941_min_high_g360}" # Default Medaka model is for MinION (or GridION) R9.4.1 flowcells using the fast Guppy basecaller version 3.6.0, high accuracy base calling
: "${max_read_length:=600}" # Filters out all reads above 600 bp
: "${min_read_quality:=7}" # Filters out all reads below a quality of 7
: "${normalized_coverage:=200}" Normalize to coverage of 200x to reduce runtime
: "${threads:=4}" # Number of threads for running the minimap2/medaka pipeline


#### 1. Length and quality filtering

source activate artic-ncov2019

# Guppyplex aggregates pre-demultiplexed reads (fastq). We use it here to filter by size (removing obviously chimeric reads) and quality (default cutoff score is 7).
# For all *fastq file in <directory>, this creates a single, aggregated fastq file called "guppyplexed_<directory>.fastq" in the current directory.
# We are assuming that only "pass" reads are in the base-called fastq files. If not, then remove --skip-quality-check (takes much longer).

artic guppyplex --skip-quality-check --directory ${samplename} --prefix guppyplexed --min-length 400 --max-length ${max_read_length} --quality ${min_read_quality}
mv guppyplexed_${samplename}.fastq ${samplename}/

cd ${samplename}/


#### 2. Alignment, variant calling, and consensus creation
# minimap2 alignment, Medaka for consensus creation and variant calling (experimental)

artic minion --medaka --medaka-model ${medaka_model} --normalise ${normalized_coverage} --threads ${threads} --scheme-directory ../artic-ncov2019/primer_schemes --read-file guppyplexed_${samplename}.fastq --no-longshot --strict nCoV-2019/V3 ${samplename}

source deactivate

source activate pangolin


# Move the reference fasta and gtf
cd /
mkdir /snpEff/data
mkdir /snpEff/data/MN908947.3
mv nCoV-2019.reference.fasta /snpEff/data/MN908947.3/sequences.fa
mv nCoV-2019.reference.gtf /snpEff/data/MN908947.3/genes.gtf

# Add the custom genome/gtf to a config file
cd snpEff
echo "MN908947.3.genome : nCoV-2019 ARTIC V3" >> snpEffect.config

# Build the custom database
java -Xmx4g -jar snpEff.jar build -c snpEffect.config -noGenome -gtf22 -v MN908947.3

# Run snpEff annotation
java -Xmx4g -jar snpEff.jar ann MN908947.3 -v -c snpEffect.config -fastaProt ../${samplename}/${samplename}.pass.snpeff.vcf.faa -csvStats ../${samplename}/${samplename}.pass.snpeff.vcf.stats ../${samplename}/${samplename}.pass.vcf > ../${samplename}/${samplename}.pass.snpeff.vcf

# Generate a bcf file with allele frequency information
cd ../${samplename}
bcftools mpileup --annotate FORMAT/AD,INFO/AD --fasta-ref ../artic-ncov2019/primer_schemes/nCoV-2019/V3/nCoV-2019.reference.fasta --max-depth 0 --count-orphans --no-BAQ --min-BQ 0 ${samplename}.sorted.bam | bcftools call -m -o ${samplename}.sorted.bcf
# In the bcf file, AD = comma-separated list of reads mapping to each allele, DP = total reads mapping to the variant position


#### 4. Pull out AA mutations, allele frequencies

# Pull out AA mutations (for generate_tsv.py)
grep "Variant" ${samplename}.pass.snpeff.vcf.faa > ${samplename}.pass.snpeff.vcf.faa.headers

# Remove header lines from the vcf and bcf files (for generate_tsv.py)
grep -v "^#" ${samplename}.pass.snpeff.vcf > ${samplename}.pass.snpeff.vcf.noheaders
grep -v "^#" ${samplename}.sorted.bcf > ${samplename}.sorted.bcf.noheaders

# Run generate_tsv.py
mv /generate_tsv.py /${samplename}
python generate_tsv.py ${samplename}


#### 5. Classify lineages

# Assign NextClade clade
nextclade --input-fasta ${samplename}.consensus.fasta --output-tsv ${samplename}.consensus.fasta.nextclade.tsv

# Assign Pango lineage
pangolin --update # being intensively updated
pangolin ${samplename}.consensus.fasta --outfile ${samplename}.consensus.fasta.pangolin.csv -t 10

source deactivate
