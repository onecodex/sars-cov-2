#!/bin/bash

set -euo pipefail

image=covid19

docker build --tag "${image}" .

docker \
  run \
  --rm \
  --workdir /data \
  --volume `pwd`:/data \
  --entrypoint /bin/bash \
  --env prefix=helloworld  \
  --env reference=reference/nCoV-2019.reference.fasta \
  --env input_fastq=data/twist-target-capture/RNA_control_spike_in_10_6_100k_reads.fastq.gz \
  --env primer_bed_file=reference/artic-v1/ARTIC-V1.bed \
  "${image}" \
  covid19_call_variants.sh
