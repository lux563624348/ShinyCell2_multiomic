#!/usr/bin/env python3
"""
Convert an h5mu object with an ATAC-by-bin modality into a pseudo fragment TSV.

Usage
-----
python h5mu_to_pseudofragments.py --h5mu HT.h5mu --out pseudo_HT_fragments.tsv

Optional arguments exist to override the modality key or matrix layer when the
defaults do not match the input file.
"""
import argparse
import re
from pathlib import Path
from typing import Iterable

import muon as mu
import numpy as np
import scipy.sparse as sp


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Extract a pseudo fragment table from an h5mu file by expanding the "
            "non-zero entries of an ATAC-by-bin matrix."
        )
    )
    parser.add_argument(
        "--h5mu",
        required=True,
        help="Path to the input h5mu file.",
    )
    parser.add_argument(
        "--out",
        required=True,
        help="Path to the TSV file to write. Existing files are overwritten.",
    )
    parser.add_argument(
        "--modality",
        default=None,
        help=(
            "Modality key that stores the ATAC-by-bin matrix. When omitted the "
            "script tries common options such as 'atac_cell_by_bin' or 'atac'."
        ),
    )
    return parser.parse_args()


def pick_modality(mdata, override: str | None) -> "anndata.AnnData":
    candidate_keys = ["atac_cell_by_bin"]
    if override:
        candidate_keys = [override]
    for key in candidate_keys:
        if key in mdata.mod:
            return mdata.mod[key]
    # fallback to the first modality
    return next(iter(mdata.mod.values()))


def parse_bins(bin_labels: Iterable[str]) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    bins = np.asarray(bin_labels, dtype=object)
    split = np.char.split(bins.astype(str))

    def to_int(x: str) -> int:
        try:
            return int(float(x))
        except Exception:
            match = re.search(
                r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?",
                str(x),
            )
            return int(float(match.group(0))) if match else 0

    chrom = np.array([parts[0] if parts else "unk" for parts in split], dtype=object)
    start = np.array(
        [
            to_int(parts[2]) if len(parts) >= 3 else 0
            for parts in split
        ],
        dtype=np.int64,
    )
    end = start + 1

    return chrom, start, end


def write_fragments(
    coo: sp.coo_matrix,
    chrom: np.ndarray,
    start: np.ndarray,
    end: np.ndarray,
    cells: np.ndarray,
    out_path: Path,
    replicate_counts: bool,
) -> int:
    total = 0
    with out_path.open("w") as handle:
        for row, col, value in zip(coo.row, coo.col, coo.data):
            count = int(value)
            if count <= 0:
                continue
            if replicate_counts:
                for _ in range(count):
                    handle.write(
                        f"{chrom[row]}\t{start[row]}\t{end[row]}\t{cells[col]}\t1\n"
                    )
                    total += 1
            else:
                handle.write(
                    f"{chrom[row]}\t{start[row]}\t{end[row]}\t{cells[col]}\t{count}\n"
                )
                total += 1
    return total


def main() -> None:
    args = parse_args()

    h5mu_path = Path(args.h5mu)
    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    print(f"Loading MuData from {h5mu_path} ...")
    mdata = mu.read_h5mu(h5mu_path)

    adata = pick_modality(mdata, args.modality)
    print(f"Selected modality with {adata.n_obs} cells and {adata.n_vars} bins.")

    # Get sparse count matrix
    if "counts" in adata.layers:
        mat = adata.layers["counts"]
    else:
        mat = adata.X
    if not sp.issparse(mat):
        mat = sp.csr_matrix(mat)

    # Row/col names
    bins = np.array(adata.var_names)  # rows
    cells = np.array(adata.obs_names)  # cols

    # split into at least 3 fields; pad if shorter
    parts = np.char.split(bins.astype(str))
    parts = np.array([p if len(p) >= 3 else (p + ["0"]*(3-len(p))) for p in parts], dtype=object)
    chrom = np.array([p[0] for p in parts], dtype=object)

    def to_int(x):
        try:
            return int(float(x))
        except Exception:
            # fallback: extract first number from the string
            m = re.search(r'[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?', str(x))
            return int(float(m.group(0))) if m else 0

    start = np.array([to_int(p[2]) for p in parts], dtype=np.int64)
    #end = start + 500

    print(f"Writing pseudo fragments to {out_path} ...")

    # Iterate over nonzeros    ## final output
    coo = mat.tocoo()
    
    limit = 100000
    n_written = 0
    stop = False
    with out_path.open("a") as f:
        for c, r, v in zip(coo.row, coo.col, coo.data):
            v = int(v)
            if v <= 0:
                continue
            else:
                f.write(f"{chrom[r]}\t{start[r]}\t{start[r]+500}\t{cells[c]}\t{v}\n")
                n_written += 1
                if n_written >= limit:
                        stop = True
                        break
            if stop:
                break

    print(f"Wrote {n_written} fragments to {out_path}")

if __name__ == "__main__":
    main()
