#!/usr/bin/env python3
"""
Convert an h5mu ATAC-by-bin modality into a bgzip-compressed, tabix-indexed
pseudo fragment file. Coordinates tolerate scientific notation (e.g. 1e+05).

Usage:
  python h5mu_to_pseudofragments_new.py --h5mu HT.h5mu --out pseudo_HT_fragments.tsv.gz
"""

import argparse
import re
import shutil
import subprocess
from pathlib import Path

import muon as mu
import numpy as np
import scipy.sparse as sp


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Expand an ATAC-by-bin matrix into a pseudo fragment table."
    )
    parser.add_argument("--h5mu", required=True, help="Input .h5mu file")
    parser.add_argument(
        "--out",
        required=True,
        help="Output path (compressed with bgzip; .gz/.bgz added if missing)",
    )
    parser.add_argument(
        "--modality",
        default=None,
        help="Modality key; defaults to 'atac_cell_by_bin' or the first modality.",
    )
    return parser.parse_args()


def _require_cmd(name: str):
    if shutil.which(name) is None:
        raise RuntimeError(
            f"Required command '{name}' not found. Install htslib/samtools for bgzip/tabix."
        )


def _bgzip_and_tabix(src: Path, dest: Path):
    """Compress src -> dest with bgzip and create a tabix index."""
    _require_cmd("bgzip")
    _require_cmd("tabix")

    print(f"Compressing {src} -> {dest} with bgzip ...")
    with dest.open("wb") as out_f:
        subprocess.run(["bgzip", "-c", str(src)], stdout=out_f, check=True)

    print(f"Indexing {dest} with tabix ...")
    subprocess.run(
        ["tabix", "-s", "1", "-b", "2", "-e", "3", "-0", "-f", str(dest)],
        check=True,
    )

    src.unlink(missing_ok=True)


def _to_int_array(values):
    """Convert to int with float()+regex fallback to handle scientific notation."""
    def _to_int(x):
        try:
            return int(float(x))
        except Exception:
            match = re.search(
                r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?",
                str(x),
            )
            if match:
                return int(float(match.group(0)))
            raise

    arr = np.asarray(values, dtype=object)
    return np.array([_to_int(v) for v in arr], dtype=np.int64)


def pick_modality(mdata, override: str | None):
    candidates = ["atac_cell_by_bin"]
    if override:
        candidates = [override]
    for key in candidates:
        if key in mdata.mod:
            return mdata.mod[key]
    return next(iter(mdata.mod.values()))


def load_correct_bins(adata):
    var = adata.var

    # --- Case 1: interval column "chr start end" ---
    if "interval" in var.columns:
        intervals = var["interval"].astype(str).to_numpy()
        parts = np.array([s.split() for s in intervals])
        chrom = parts[:, 0]
        start = _to_int_array(parts[:, 1])
        end = _to_int_array(parts[:, 2])
        return chrom, start, end

    # --- Case 2: chrom + start (+ end) columns exist ---
    if "chrom" in var.columns and "start" in var.columns:
        chrom = var["chrom"].astype(str).to_numpy()
        start = _to_int_array(var["start"].to_numpy())
        if "end" in var.columns:
            end = _to_int_array(var["end"].to_numpy())
        else:
            end = start + 499
        return chrom, start, end

    # --- Case 3: fallback → parse var_names: "chr start end" ---
    names = adata.var_names.astype(str)
    parts = np.array([s.split() for s in names])
    if parts.shape[1] < 3:
        raise ValueError("Cannot parse bin coordinates from var_names.")

    chrom = parts[:, 0]
    start = _to_int_array(parts[:, 2])
    end = start + 499
    return chrom, start, end


def main():
    args = parse_args()

    h5mu_path = Path(args.h5mu)
    out_path = Path(args.out)

    # Normalize output to .gz/.bgz
    if out_path.suffix not in {".gz", ".bgz"}:
        gz_path = out_path.with_suffix(out_path.suffix + ".gz")
        print(f"Output will be bgzip-compressed to {gz_path} (requested {out_path}).")
    else:
        gz_path = out_path
    gz_path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = gz_path.with_name(gz_path.name + ".tmp")

    print(f"Loading MuData from {h5mu_path} ...")
    mdata = mu.read_h5mu(h5mu_path)
    adata = pick_modality(mdata, args.modality)
    print(f"Selected modality with {adata.n_obs} cells and {adata.n_vars} bins.")

    # counts matrix
    mat = adata.layers.get("counts", adata.X)
    if not sp.issparse(mat):
        mat = sp.csr_matrix(mat)
    coo = mat.tocoo()

    chrom, start, end = load_correct_bins(adata)
    cells = adata.obs_names.to_numpy()

    print(f"Writing fragments to temporary file: {tmp_path}")
    n_written = 0
    with tmp_path.open("w") as f:
        for r, c, v in zip(coo.row, coo.col, coo.data):
            v = int(v)
            if v <= 0:
                continue
            # matrix is cells x bins: row -> cell, col -> bin
            f.write(f"{chrom[c]}\t{start[c]}\t{end[c]}\t{cells[r]}\t{v}\n")
            n_written += 1

    print(f"Done. Wrote {n_written} fragments. Compressing and indexing ...")
    _bgzip_and_tabix(tmp_path, gz_path)
    print(f"Finished. Compressed file: {gz_path} (index: {gz_path}.tbi)")


if __name__ == "__main__":
    main()
