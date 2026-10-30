#!/usr/bin/env python3

"""Filter a brain atlas MuData object by enrollment and batch size.

This script:
1. Reads a MuData .h5mu file.
2. Reads brain atlas sample metadata.
3. Left-joins metadata into ``mdata.obs`` using ``obs['batch']`` -> ``DONOR_ID``.
4. Removes cells with ``Enrollment == 'N'`` by slicing the whole MuData object.
5. Removes batches with fewer than ``--min-cells-per-batch`` cells, again by
   slicing the whole MuData object.
6. Writes the filtered MuData object to a new .h5mu file.
7. Reopens the output in writable backed mode and removes stale
   ``mdata['rna'].uns['mask_*']`` entries whose suffixes no longer occur in
   ``mdata['rna'].obs['batch']``.

The filtering steps intentionally operate on ``mdata`` itself, not just
``mdata.obs``, to keep global and modality-specific observations synchronized.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import mudata as mu
import pandas as pd


DEFAULT_INPUT = Path("data/processed/Atlas_20260528/aligned.h5mu")
DEFAULT_METADATA = Path("parameters/2026-05-24_brain_atlas.csv")
DEFAULT_OUTPUT = Path("data/processed/Atlas_20260528/aligned_enrolled_min1000.h5mu")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Join brain atlas metadata to a MuData object, remove "
            "Enrollment == 'N' cells, and remove small batches."
        )
    )
    parser.add_argument(
        "--input",
        default=DEFAULT_INPUT,
        type=Path,
        help=f"Input MuData .h5mu path. Default: {DEFAULT_INPUT}",
    )
    parser.add_argument(
        "--metadata",
        default=DEFAULT_METADATA,
        type=Path,
        help=f"Brain atlas metadata CSV path. Default: {DEFAULT_METADATA}",
    )
    parser.add_argument(
        "--output",
        default=DEFAULT_OUTPUT,
        type=Path,
        help=f"Output filtered .h5mu path. Default: {DEFAULT_OUTPUT}",
    )
    parser.add_argument(
        "--min-cells-per-batch",
        default=1000,
        type=int,
        help="Remove batches with fewer than this many cells. Default: 1000",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Run checks and filtering summaries, but do not write output.",
    )
    return parser.parse_args()


def shape_summary(mdata: mu.MuData) -> dict[str, object]:
    """Return global and modality-specific shapes."""
    return {
        "shape": mdata.shape,
        "mod_shapes": {name: mod.shape for name, mod in mdata.mod.items()},
    }


def print_shape_summary(label: str, mdata: mu.MuData) -> None:
    """Print global and modality-specific shapes."""
    shapes = shape_summary(mdata)
    print(f"{label} MuData shape: {shapes['shape']}")
    for name, shape in shapes["mod_shapes"].items():
        print(f"{label} {name} shape: {shape}")


def check_required_columns(
    mdata: mu.MuData,
    metadata: pd.DataFrame,
    batch_col: str = "batch",
    donor_col: str = "DONOR_ID",
    enrollment_col: str = "Enrollment",
) -> None:
    """Validate columns needed by the filtering workflow."""
    if batch_col not in mdata.obs.columns:
        raise ValueError(f"mdata.obs is missing required column: {batch_col}")
    if donor_col not in metadata.columns:
        raise ValueError(f"metadata is missing required column: {donor_col}")
    if enrollment_col not in metadata.columns:
        raise ValueError(f"metadata is missing required column: {enrollment_col}")


def join_metadata(
    mdata: mu.MuData,
    metadata: pd.DataFrame,
    batch_col: str = "batch",
    donor_col: str = "DONOR_ID",
) -> mu.MuData:
    """Left-join sample metadata into mdata.obs while preserving obs index."""
    duplicated_donor_rows = int(metadata[donor_col].duplicated(keep=False).sum())
    if duplicated_donor_rows > 0:
        raise ValueError(
            f"metadata[{donor_col!r}] must be unique for many-to-one join; "
            f"found {duplicated_donor_rows} duplicated rows."
        )

    unmatched_batches = sorted(
        set(mdata.obs[batch_col].dropna().unique())
        - set(metadata[donor_col].dropna().unique())
    )
    if unmatched_batches:
        preview = unmatched_batches[:20]
        raise ValueError(
            f"Found {len(unmatched_batches)} {batch_col!r} values without "
            f"metadata match in {donor_col!r}. Preview: {preview}"
        )

    obs_index_name = mdata.obs.index.name
    obs_index_before = mdata.obs.index.copy()
    n_obs_before = len(mdata.obs)

    mdata.obs = (
        mdata.obs.reset_index(names="_obs_index")
        .merge(
            metadata,
            how="left",
            left_on=batch_col,
            right_on=donor_col,
            validate="many_to_one",
        )
        .set_index("_obs_index")
    )
    mdata.obs.index.name = obs_index_name

    if len(mdata.obs) != n_obs_before:
        raise RuntimeError(
            f"Join changed obs rows from {n_obs_before} to {len(mdata.obs)}."
        )
    if not mdata.obs.index.equals(obs_index_before):
        raise RuntimeError("Join did not preserve the original mdata.obs index.")

    unmatched_rows = int(mdata.obs[donor_col].isna().sum())
    if unmatched_rows > 0:
        raise RuntimeError(f"Join produced {unmatched_rows} unmatched obs rows.")

    print("Metadata join complete.")
    print(f"  obs rows: {n_obs_before}")
    print(f"  unmatched rows after join: {unmatched_rows}")
    return mdata


def filter_enrollment(
    mdata: mu.MuData,
    enrollment_col: str = "Enrollment",
    remove_value: str = "N",
) -> mu.MuData:
    """Remove cells with Enrollment == remove_value by slicing MuData."""
    counts = mdata.obs[enrollment_col].value_counts(dropna=False).to_dict()
    remove_mask = mdata.obs[enrollment_col].eq(remove_value)
    n_remove = int(remove_mask.sum())
    n_keep = int((~remove_mask).sum())

    print("Enrollment filter:")
    print(f"  counts before: {counts}")
    print(f"  cells to remove ({enrollment_col} == {remove_value!r}): {n_remove}")
    print(f"  cells to keep: {n_keep}")
    print_shape_summary("Before enrollment filter", mdata)

    keep_mask = ~remove_mask.to_numpy()
    mdata = mdata[keep_mask, :]

    remaining_removed = int(mdata.obs[enrollment_col].eq(remove_value).sum())
    if remaining_removed != 0:
        raise RuntimeError(
            f"Enrollment filter failed: {remaining_removed} cells with "
            f"{enrollment_col} == {remove_value!r} remain."
        )

    print_shape_summary("After enrollment filter", mdata)
    print(f"  remaining {enrollment_col} == {remove_value!r}: {remaining_removed}")
    return mdata


def filter_small_batches(
    mdata: mu.MuData,
    min_cells_per_batch: int,
    batch_col: str = "batch",
) -> mu.MuData:
    """Remove batches with fewer than min_cells_per_batch cells."""
    if min_cells_per_batch < 1:
        raise ValueError("--min-cells-per-batch must be >= 1")

    batch_counts = mdata.obs[batch_col].value_counts(dropna=False)
    small_batches = batch_counts[batch_counts < min_cells_per_batch]

    print("Small-batch filter:")
    print(f"  batches before: {int(batch_counts.shape[0])}")
    print(f"  min cells per batch threshold: {min_cells_per_batch}")
    print(f"  small batches to remove: {int(small_batches.shape[0])}")
    print(f"  cells to remove: {int(small_batches.sum())}")
    print(f"  cells to keep: {int(mdata.n_obs - small_batches.sum())}")
    print_shape_summary("Before small-batch filter", mdata)

    keep_mask = ~mdata.obs[batch_col].isin(small_batches.index).to_numpy()
    mdata = mdata[keep_mask, :]

    batch_counts_after = mdata.obs[batch_col].value_counts(dropna=False)
    remaining_small_batches = int((batch_counts_after < min_cells_per_batch).sum())
    remaining_small_batch_cells = int(
        batch_counts_after[batch_counts_after < min_cells_per_batch].sum()
    )

    if remaining_small_batches != 0:
        raise RuntimeError(
            f"Small-batch filter failed: {remaining_small_batches} batches "
            f"below {min_cells_per_batch} remain."
        )

    print_shape_summary("After small-batch filter", mdata)
    print(f"  batches after: {int(batch_counts_after.shape[0])}")
    print(f"  min cells per batch after: {int(batch_counts_after.min())}")
    print(f"  remaining batches below threshold: {remaining_small_batches}")
    print(f"  remaining cells in batches below threshold: {remaining_small_batch_cells}")
    return mdata


def find_stale_rna_uns_mask_keys(
    mdata: mu.MuData,
    modality: str = "rna",
    batch_col: str = "batch",
    mask_prefix: str = "mask_",
) -> list[str]:
    """Return RNA mask keys whose suffix is absent from obs batch values."""
    if modality not in mdata.mod:
        raise ValueError(f"mdata is missing required modality: {modality!r}")

    mod = mdata.mod[modality]
    if batch_col not in mod.obs.columns:
        raise ValueError(f"mdata[{modality!r}].obs is missing column: {batch_col!r}")

    mask_keys = [str(key) for key in mod.uns if str(key).startswith(mask_prefix)]
    batch_values = set(mod.obs[batch_col].dropna().astype(str))
    return sorted(
        key
        for key in mask_keys
        if key.removeprefix(mask_prefix) not in batch_values
    )


def report_stale_rna_uns_masks(mdata: mu.MuData) -> None:
    """Report stale RNA masks without modifying a filtered MuData view."""
    stale_keys = find_stale_rna_uns_mask_keys(mdata)
    print("RNA uns mask cleanup preview:")
    print(f"  stale mask entries to remove: {len(stale_keys)}")
    print(f"  preview: {stale_keys[:20]}")


def clean_output_rna_uns_masks(output: Path, modality: str = "rna") -> None:
    """Reopen an output file in writable backed mode and remove stale RNA masks."""
    mdata = mu.read_h5mu(output, backed="r+")
    try:
        stale_keys = find_stale_rna_uns_mask_keys(mdata, modality=modality)
        mod = mdata.mod[modality]
        for key in stale_keys:
            mod.uns.pop(key, None)

        mdata.write()
    finally:
        mdata.file.close()

    verified = mu.read_h5mu(output, backed="r")
    try:
        remaining_stale_keys = find_stale_rna_uns_mask_keys(
            verified,
            modality=modality,
        )
        remaining_mask_count = sum(
            str(key).startswith("mask_") for key in verified.mod[modality].uns
        )
    finally:
        verified.file.close()

    if remaining_stale_keys:
        raise RuntimeError(
            f"RNA uns mask cleanup failed: {len(remaining_stale_keys)} stale "
            f"mask entries remain. Preview: {remaining_stale_keys[:20]}"
        )

    print("RNA uns mask cleanup:")
    print(f"  stale mask entries removed: {len(stale_keys)}")
    print(f"  mask entries after: {remaining_mask_count}")
    print("  stale mask entries after cleanup: 0")


def write_mudata(mdata: mu.MuData, output: Path) -> None:
    """Write a filtered MuData object or view to a new file."""
    output.parent.mkdir(parents=True, exist_ok=True)
    mdata.write(output)

def main() -> None:
    args = parse_args()

    print(f"Loading MuData: {args.input}")
    mdata = mu.read_h5mu(args.input, backed="r")
    print(f"  backed: {getattr(mdata, 'isbacked', None)}")
    print_shape_summary("Initial", mdata)

    print(f"\nLoading metadata: {args.metadata}")
    metadata = pd.read_csv(args.metadata)
    print(f"  metadata rows: {len(metadata)}")
    print(f"  metadata columns: {list(metadata.columns)}")

    check_required_columns(mdata, metadata)

    print("\nJoining metadata...")
    mdata = join_metadata(mdata, metadata)

    print("\nFiltering Enrollment == 'N' cells...")
    mdata = filter_enrollment(mdata)

    print("\nFiltering small batches...")
    mdata = filter_small_batches(mdata, args.min_cells_per_batch)

    if args.dry_run:
        print("\nChecking stale RNA uns mask entries...")
        report_stale_rna_uns_masks(mdata)
        print("\nDry run requested; output was not written.")
        return

    if args.input.resolve() == args.output.resolve():
        raise ValueError("--output must differ from --input")

    print(f"\nWriting filtered MuData: {args.output}")
    write_mudata(mdata, args.output)
    mdata.file.close()
    print(f"Saved filtered MuData: {args.output}")

    print("\nReopening output and cleaning stale RNA uns mask entries...")
    clean_output_rna_uns_masks(args.output)
    print(f"Saved cleaned MuData: {args.output}")


if __name__ == "__main__":
    main()
