#!/bin/bash

# To run this script:
# ./covid19_call_variants.artic.sh <samplename>

# This pipeline needs the fastq file of "pass" reads

set -euo pipefail

#samplename=$(echo $1 | cut -d. -f1) # remove ".fastq.gz" from input file to get the sample
work_dir=${PWD}
echo ${work_dir}

mv $1 sample.fastq

#### Default parameters

: "${medaka_model:=r941_min_high_g360}" # Default Medaka model is for MinION (or GridION) R9.4.1 flowcells using the fast Guppy basecaller version 3.6.0, high accuracy base calling
: "${max_read_length:=600}" # Filters out all reads above 600 bp
: "${min_read_quality:=7}" # Filters out all reads below a quality of 7
: "${normalized_coverage:=200}" # Normalize to coverage of 200x to reduce runtime
: "${threads:=6}" # Number of threads for running the minimap2/medaka pipeline


#### 1. Length and quality filtering

source activate artic-ncov2019

# Guppyplex aggregates pre-demultiplexed reads (fastq). We use it here to filter by size (removing obviously chimeric reads) and quality (default cutoff score is 7).
# For all *fastq file in <directory>, this creates a single, aggregated fastq file called "guppyplexed_<directory>.fastq" in the current directory.
# We are assuming that only "pass" reads are in the base-called fastq files. If not, then remove --skip-quality-check (takes much longer).

artic guppyplex --directory ./ --prefix guppyplexed --min-length 400 --max-length ${max_read_length} --quality ${min_read_quality}

mv guppyplexed_.fastq sample_guppyplexed.fastq


#### 2. Alignment, variant calling, and consensus creation
# minimap2 alignment, Medaka for consensus creation and variant calling (experimental)

artic minion --medaka --medaka-model ${medaka_model} --normalise ${normalized_coverage} --threads ${threads} --scheme-directory /artic-ncov2019/primer_schemes --read-file sample_guppyplexed.fastq --no-longshot --strict nCoV-2019/V3 sample

source deactivate


#### 3. Variant annotation

source activate report

# Replace the vcf chromosome name to the verison that snpEff will recognize (GenBank sequence NC_045512.2, which is identical to RefSeq MN908947.3)
gunzip sample.pass.vcf.gz
sed -i 's/MN908947.3/NC_045512.4/g' sample.pass.vcf

# Run the annotation
cd /snpEff
java -Xmx4g -jar snpEff.jar ann NC_045512.2 -verbose -fastaProt ${work_dir}/sample.pass.snpeff.vcf.faa -csvStats ${work_dir}/sample.pass.snpeff.vcf.stats ${work_dir}/sample.pass.vcf > ${work_dir}/sample.pass.snpeff.vcf

# Extract fields of interest (position, ref allele, alt allele, allele reads, AA mutation) from the vcf with SnpSift
java -jar SnpSift.jar extractFields ${work_dir}/sample.pass.snpeff.vcf CHROM POS REF ALT SR ANN[*].HGVS_P > ${work_dir}/snpeff.tsv


#### 4. Classify lineages

cd ${work_dir}

# Assign NextClade clade
nextclade --input-fasta sample.consensus.fasta --output-tsv nextclade.tsv

# Assign Pango lineage
pangolin --update # being intensively updated
pangolin sample.consensus.fasta --outfile pangolin.csv -t ${threads}

source deactivate

source activate notebook

#### 6. Run generate_tsv.py
python /usr/local/bin/generate_tsv.py sample

### Move files around
mv sample.consensus.fasta consensus.fa
mv sample.sorted.bam covid19.bam
mv sample.sorted.bam.bai covid19.bam.bai
