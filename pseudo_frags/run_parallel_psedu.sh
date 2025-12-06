#!/usr/bin/env bash

set -euo pipefail

echo "Listing *raw.h5mu files under /data (sorted by size ascending):"
raw_files=$(find /data/ -name "*raw.h5mu" -print0 | xargs -0 ls -lhSr)

if [[ -n "$raw_files" ]]; then
  printf '%s\n' "$raw_files"

  mapfile -t raw_paths < <(
    find /data -name "*raw.h5mu" -print0 \
      | xargs -0 stat --format='%s %n' \
      | LC_ALL=C sort -n -k1,1 \
      | awk '{size=$1; $1=""; print substr($0,2)}'
  )

  for h5mu in "${raw_paths[@]}"; do
    base=$(basename "${h5mu%.h5mu}")
    echo "Processing ${h5mu} → /data/pseudo_fragments/${base}_new_pseudo"
    python /data/pseudo_fragments/ShinyCell2_multiomic/pseudo_frags/h5mu_to_pseudofragments_new.py \
      --h5mu "$h5mu" \
      --out "/data/pseudo_fragments/${base}_new_pseudo"
  done
else
  echo "No *raw.h5mu files found."
fi
