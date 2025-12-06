#!/usr/bin/env bash

set -euo pipefail

echo "Listing *raw.h5mu files under /data (sorted by size):"
raw_list=$(mktemp)
find /data -name "*raw.h5mu" -print0 | while IFS= read -r -d '' f; do
  stat --format='%s %n' "$f" >> "$raw_list"
done
if [[ -s "$raw_list" ]]; then
  sort -n "$raw_list" \
    | awk '{size=$1; $1=""; printf "%10s %s\n", size, substr($0,2)}'
else
  echo "No *raw.h5mu files found."
fi
rm -f "$raw_list"
echo

find /data -name "*processed.h5mu" -print0 |
while IFS= read -r -d '' h5mu; do
  base=$(basename "${h5mu%.h5mu}")
  echo "Processing ${h5mu} → /data/pseudo_fragments/${base}_new_pseudo"
  python /data/pseudo_fragments/ShinyCell2_multiomic/pseudo_frags/h5mu_to_pseudofragments_new.py \
    --h5mu "$h5mu" \
    --out "/data/pseudo_fragments/${base}_new_pseudo"
done
