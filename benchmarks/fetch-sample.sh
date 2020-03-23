#!/bin/bash

set -euo pipefail

out_path="${1}.fastq.gz"

fastq-dump \
  --split-spot \
  --read-filter \
  pass \
  --skip-technical \
  --readids \
  --stdout \
  --dumpbase \
  --clip \
  "${1}" \
  | pigz \
  > "${out_path}"

echo "${out_path}"
