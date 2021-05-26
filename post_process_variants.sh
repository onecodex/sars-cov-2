#!/bin/bash

# Install Nextclade
npm install --global @neherlab/nextclade

# Assign NextClade clade
echo "Assigning NextClade Clade"
nextclade --input-fasta "${1}" --output-tsv nextclade.tsv --output-json nextclade.json

# TODO: copy pangolin database data somewhere.
echo "Assigning Pango Lineage"
conda run -n pangolin pangolin "${1}" --outfile pangolin.csv
