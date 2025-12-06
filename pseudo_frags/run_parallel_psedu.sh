#!/usr/bin/env bash

set -euo pipefail

echo "Listing *raw.h5mu files under /data (sorted by size):"
find /data -name "*raw.h5mu" -print0 \
  | xargs -0 stat --format='%s %n' \
  | sort -n \
  | awk '{size=$1; $1=\"\"; printf \"%10s %s\n\", size, substr($0,2)}'
echo

find /data -name "*processed.h5mu" -print0 |
while IFS= read -r -d '' h5mu; do
  base=$(basename "${h5mu%.h5mu}")
  echo "Processing ${h5mu} → /data/pseudo_fragments/${base}_new_pseudo"
  python /data/pseudo_fragments/ShinyCell2_multiomic/pseudo_frags/h5mu_to_pseudofragments_new.py \
    --h5mu "$h5mu" \
    --out "/data/pseudo_fragments/${base}_new_pseudo"
done
