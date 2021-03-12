#!/bin/bash

set -euo pipefail

: "${threads=4}"
: "${consensus_fasta=${1}}"

# inputs:
# 2. consensus.fasta

#### 4. Classify lineages

# Assign NextClade clade
conda run -n pangolin nextclade --input-fasta "${consensus_fasta}" --output-tsv nextclade.tsv

# Assign Pango lineage
conda run -n pangolin pangolin --update # being intensively updated

# TODO: copy pangolin database data somewhere.

conda run -n pangolin pangolin "${consensus_fasta}" --outfile pangolin.csv -t "${threads}"

# outputs

# 1. nexclade.tsv
# 2. pangolin.tsv
