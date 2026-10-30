#!/usr/bin/env python3

"""Merge 10x matrix directories (RNA or ATAC) into one AnnData object.

Each subdirectory under --input-dir is treated as a 10x matrix directory
(donor ID = subdirectory name). --modality controls which files are required
and whether a mask is loaded.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc


RNA_FILES = (
    "matrix.mtx.gz",
    "barcodes.tsv.gz",
    "features.tsv.gz",
    "rna_ref_presence.csv.gz",
)

ATAC_FILES = (
    "matrix.mtx.gz",
    "barcodes.tsv.gz",
    "features.tsv.gz",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Merge 10x matrix directories under --input-dir. "
            "Default is dry-run; add --run to write the merged h5ad."
        )
    )
    parser.add_argument(
        "--modality",
        required=True,
        choices=["rna", "atac"],
        help="Data modality: 'rna' or 'atac'.",
    )
    parser.add_argument(
        "--input-dir",
        default="data/processed/Atlas_20260528/rna",
        help="Directory whose subdirectories are 10x matrix directories.",
    )
    parser.add_argument(
        "--output",
        default="data/processed/Atlas_20260528/merged_rna.h5ad",
        help="Output h5ad path.",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=None,
        help="Only process the first N entries for testing.",
    )
    parser.add_argument(
        "--run",
        action="store_true",
        help="Actually read matrices, merge, and write output.",
    )
    return parser.parse_args()


def resolve_path(path: str, root_dir: Path) -> Path:
    input_path = Path(path)
    if input_path.is_absolute():
        return input_path
    return root_dir / input_path


def discover_samples(input_dir: Path) -> pd.DataFrame:
    if not input_dir.exists():
        raise FileNotFoundError(f"Input directory not found: {input_dir}")
    if not input_dir.is_dir():
        raise NotADirectoryError(f"Not a directory: {input_dir}")

    all_entries = sorted(input_dir.iterdir())
    dirs = [p for p in all_entries if p.is_dir() and not p.name.startswith(".")]
    if not dirs:
        raise ValueError(f"No subdirectories found in {input_dir}")

    n_skipped = len(all_entries) - len(dirs)
    if n_skipped:
        print(f"Skipped {n_skipped} non-directory or hidden entries.")

    return pd.DataFrame(
        {
            "DONOR_ID": [d.name for d in dirs],
            "matrix_dir": [str(d.resolve()) for d in dirs],
        }
    )


def read_features(features_path: Path) -> pd.DataFrame:
    return pd.read_csv(
        features_path,
        sep="\t",
        header=None,
        names=["id", "name", "feature_type"],
        dtype=str,
    )


def preflight(entries: pd.DataFrame, modality: str) -> pd.DataFrame:
    required = RNA_FILES if modality == "rna" else ATAC_FILES
    problems: list[dict[str, str]] = []
    first_features: pd.DataFrame | None = None
    first_donor: str | None = None

    for row in entries.itertuples(index=False):
        donor_id = row.DONOR_ID
        matrix_dir = Path(row.matrix_dir)

        for file_name in required:
            file_path = matrix_dir / file_name
            if not file_path.exists():
                problems.append(
                    {
                        "donor_id": donor_id,
                        "matrix_dir": str(matrix_dir),
                        "problem": f"missing {file_name}",
                    }
                )

        features_path = matrix_dir / "features.tsv.gz"
        if features_path.exists():
            features = read_features(features_path)
            if first_features is None:
                first_features = features
                first_donor = donor_id
            elif not features.equals(first_features):
                problems.append(
                    {
                        "donor_id": donor_id,
                        "matrix_dir": str(matrix_dir),
                        "problem": f"features.tsv.gz differs from {first_donor}",
                    }
                )

    return pd.DataFrame(problems)


def read_mask(mask_path: Path) -> np.ndarray:
    mask_table = pd.read_csv(mask_path, compression="gzip")
    if mask_table.shape[1] < 2:
        raise ValueError(f"Mask file must have at least two columns: {mask_path}")
    return mask_table.iloc[:, 1].to_numpy(dtype=np.float32)


def read_one_sample(
    matrix_dir: Path, donor_id: str, modality: str
) -> tuple[ad.AnnData, np.ndarray | None]:
    gex_only = modality == "rna"
    adata = sc.read_10x_mtx(
        matrix_dir,
        var_names="gene_symbols",
        make_unique=False,
        cache=False,
        gex_only=gex_only,
    )

    original_barcodes = adata.obs_names.astype(str)
    adata.obs["barcode"] = original_barcodes
    adata.obs["batch"] = donor_id
    adata.obs_names = [f"{donor_id}_{barcode}" for barcode in original_barcodes]

    if modality == "rna":
        mask = read_mask(matrix_dir / "rna_ref_presence.csv.gz")
        if len(mask) != adata.n_vars:
            raise ValueError(
                f"Mask length mismatch for {donor_id}: mask={len(mask)}, n_vars={adata.n_vars}"
            )
        return adata, mask

    return adata, None


def merge_samples(entries: pd.DataFrame, modality: str) -> ad.AnnData:
    adatas: list[ad.AnnData] = []
    masks: dict[str, np.ndarray] = {}
    first_var: pd.DataFrame | None = None

    for i, row in enumerate(entries.itertuples(index=False), start=1):
        donor_id = row.DONOR_ID
        matrix_dir = Path(row.matrix_dir)
        print(f"[{i}/{len(entries)}] reading {donor_id}: {matrix_dir}")

        adata, mask = read_one_sample(matrix_dir, donor_id, modality)
        if first_var is None:
            first_var = adata.var.copy()

        adatas.append(adata)
        if mask is not None:
            masks[f"mask_{donor_id}"] = mask

    combined = ad.concat(
        adatas,
        axis=0,
        join="inner",
        merge="same",
        uns_merge=None,
    )
    if first_var is not None:
        combined.var = first_var.loc[combined.var_names].copy()
    if masks:
        combined.uns.update(masks)

    return combined


def main() -> None:
    args = parse_args()
    root_dir = Path.cwd().resolve()
    input_dir = resolve_path(args.input_dir, root_dir)
    output_path = resolve_path(args.output, root_dir)

    entries = discover_samples(input_dir)
    if args.limit is not None:
        if args.limit <= 0:
            raise ValueError("--limit must be a positive integer.")
        entries = entries.head(args.limit).copy()

    print(f"Found {len(entries)} directories.")
    print(f"Modality: {args.modality}")
    print(f"Output: {output_path}")
    print(f"Mode: {'run' if args.run else 'dry-run'}")

    problems = preflight(entries, args.modality)
    if not problems.empty:
        print("\nPreflight failed:")
        print(problems.to_string(index=False))
        raise SystemExit(1)

    print("Preflight passed.")
    print(entries[["DONOR_ID", "matrix_dir"]].to_string(index=False))

    if not args.run:
        print("Dry-run only. Add --run to merge and write the h5ad file.")
        return

    output_path.parent.mkdir(parents=True, exist_ok=True)
    combined = merge_samples(entries, args.modality)
    combined.write_h5ad(output_path)

    print(f"Saved: {output_path}")
    print(f"Shape: {combined.n_obs} cells x {combined.n_vars} genes")
    print(f"obs columns: {list(combined.obs.columns)}")
    print(f"uns masks: {[key for key in combined.uns.keys() if key.startswith('mask_')]}")


if __name__ == "__main__":
    main()
