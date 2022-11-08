#!/bin/bash

echo "Assigning Pango Lineage"
conda run -n pangolin pangolin "${1}" --outfile pangolin.csv

# Assign NextClade clade
echo "Assigning NextClade Clade"

/usr/local/bin/nextclade run "${1}" --input-dataset /usr/local/bin/data/sars-cov-2 --output-tsv nextclade.tsv --output-json nextclade.json   

