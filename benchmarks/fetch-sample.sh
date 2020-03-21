#!/bin/bash

set -euo pipefail

fastq-dump \
  --split-spot \
  --read-filter \
  pass \
  --skip-technical \
  --readids \
  --stdout \
  --dumpbase \
  --clip \
  "${1}"
