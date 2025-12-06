#!/usr/bin/env python3
"""
Convert an h5mu object with an ATAC-by-bin modality into a pseudo fragment TSV.

Usage:
    python h5mu_to_pseudofragments.py --h5mu HT.h5mu --out pseudo_HT_fragments.tsv
"""

import argparse
import re
from pathlib import Path

import muon as mu
import numpy as np
import scipy.sparse as sp


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--h5mu", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--modality", default=None)
    return p.parse_args()


def pick_modality(mdata, override):
    if override:
        if override in mdata.mod:
            return mdata.mod[override]
        raise ValueError(f"Modality {override} not found.")
    if "atac_cell_by_bin" in mdata.mod:
        return mdata.mod["atac_cell_by_bin"]
    return next(iter(mdata.mod.values()))


# --- chromosome cleaner -------------------------------------------------------

def clean_chr(x: str) -> str:
    """
    Normalize chromosome labels:
        'chr1' → 'chr1'
        '1' → 'chr1'
        'chrchr1' → 'chr1'
        weird values → 'chrUn'
    """
    x = str(x)

    # remove repeated chr
    x = re.sub(r"^chrchr", "chr", x)

    # keep chr prefix
    if x.startswith("chr"):
        return x

    # if pure number → chrN
    if x.isdigit():
        return "chr" + x

    # extract numeric part
    m = re.match(r"(\d+)", x)
    if m:
        return "chr" + m.group(1)

    # no match → unassigned
    return "chrUn"


def extract_start(x: str) -> int:
    """Extract integer start coordinate safely."""
    try:
        return max(0, int(float(x)))
    except Exception:
        m = re.search(r'\d+', str(x))
        return int(m.group(0)) if m else 0


# ------------------------------------------------------------------------------


def main():
    args = parse_args()
    out_path = Path(args.out)

    print(f"Loading {args.h5mu} ...")
    mdata = mu.read_h5mu(args.h5mu)
    adata = pick_modality(mdata, args.modality)

    # matrix
    mat = adata.layers["counts"] if "counts" in adata.layers else adata.X
    if not sp.issparse(mat):
        mat = sp.csr_matrix(mat)

    bins = np.array(adata.var_names)
    cells = np.array(adata.obs_names)

    # parse bin labels → chrom, start, end
    parts = np.char.split(bins.astype(str))
    parts = np.array([p if len(p) >= 3 else (p + ["0"]*(3-len(p))) for p in parts],
                     dtype=object)

    chrom = np.array([clean_chr(p[0]) for p in parts], dtype=object)
    start = np.array([extract_start(p[2]) for p in parts], dtype=np.int64)
    end = start + 500

    # Everything goes to COO format
    coo = mat.tocoo()

    # --- SORT for BGZIP/TABIX -----------------------------------------------
    print("Sorting fragments for tabix ...")
    rows = []
    for r, c, v in zip(coo.row, coo.col, coo.data):
        v = int(v)
        if v > 0:
            rows.append((chrom[r], start[r], end[r], cells[c], v))

    # Convert to array for sorting
    rows = np.array(rows, dtype=object)

    # Sort: chrom → start → end
    order = np.lexsort((rows[:,2].astype(np.int64),
                        rows[:,1].astype(np.int64),
                        rows[:,0]))
    rows = rows[order]

    # --- WRITE ---------------------------------------------------------------
    print(f"Writing sorted fragments to {out_path} ...")

    with out_path.open("w") as f:
        for chr_, s, e, cell, v in rows:
            f.write(f"{chr_}\t{s}\t{e}\t{cell}\t{v}\n")

    print(f"Done. Wrote {rows.shape[0]} fragments.")
    print("Now run:")
    print(f"  bgzip -f {out_path}")
    print(f"  tabix -s 1 -b 2 -e 3 -0 -f {out_path}.gz")


if __name__ == "__main__":
    main()
