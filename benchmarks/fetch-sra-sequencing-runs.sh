#!/bin/bash

set -euo pipefail

bq \
  query \
  --format json \
  --nouse_legacy_sql \
  --max_rows 10000 \
  'select * from `nih-sra-datastore`.sra.metadata where organism = "Severe acute respiratory syndrome coronavirus 2";' \
  > sra-sequence-runs.json
