#!/usr/bin/env python3

"""Align RNA and ATAC AnnData objects by batch and create a MuData object.

For batches present in both RNA and ATAC (shared batches), keep only the
intersection of cells (by obs_names). For batches unique to one modality,
keep all cells. Output is a MuData object (.h5mu).
"""

from __future__ import annotations

import argparse
from pathlib import Path

import anndata as ad
import mudata as md
import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Align RNA and ATAC AnnData by batch. "
            "Shared batches: keep cell intersection. "
            "Non-shared batches: keep all cells."
        )
    )
    parser.add_argument(
        "--rna",
        required=True,
        help="Path to RNA h5ad file.",
    )
    parser.add_argument(
        "--atac",
        required=True,
        help="Path to ATAC h5ad file.",
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Output path for MuData object (.h5mu).",
    )
    return parser.parse_args()


def align_rna_atac(rna: ad.AnnData, atac: ad.AnnData) -> tuple[ad.AnnData, ad.AnnData]:
    """Align RNA and ATAC by batch, keeping cell intersection for shared batches."""
    rna_batches = set(rna.obs["batch"].unique())
    atac_batches = set(atac.obs["batch"].unique())
    
    shared_batches = rna_batches & atac_batches
    rna_only_batches = rna_batches - atac_batches
    atac_only_batches = atac_batches - rna_batches
    
    print(f"RNA batches: {len(rna_batches)}")
    print(f"ATAC batches: {len(atac_batches)}")
    print(f"Shared batches: {len(shared_batches)}")
    print(f"RNA-only batches: {len(rna_only_batches)}")
    print(f"ATAC-only batches: {len(atac_only_batches)}")
    
    # Process shared batches: keep intersection
    rna_indices_shared = []
    atac_indices_shared = []
    
    for batch in sorted(shared_batches):
        rna_mask = rna.obs["batch"] == batch
        atac_mask = atac.obs["batch"] == batch
        
        rna_obs = set(rna.obs_names[rna_mask])
        atac_obs = set(atac.obs_names[atac_mask])
        
        intersection = rna_obs & atac_obs
        
        rna_indices_shared.extend(
            idx for idx in rna.obs_names[rna_mask] if idx in intersection
        )
        atac_indices_shared.extend(
            idx for idx in atac.obs_names[atac_mask] if idx in intersection
        )
    
    print(f"\nShared batch cells: {len(rna_indices_shared)} RNA, {len(atac_indices_shared)} ATAC")
    
    # Process non-shared batches: keep all
    rna_mask_only = rna.obs["batch"].isin(rna_only_batches)
    atac_mask_only = atac.obs["batch"].isin(atac_only_batches)
    
    rna_indices_only = rna.obs_names[rna_mask_only].tolist()
    atac_indices_only = atac.obs_names[atac_mask_only].tolist()
    
    print(f"Non-shared batch cells: {len(rna_indices_only)} RNA, {len(atac_indices_only)} ATAC")
    
    # Combine indices
    rna_indices_final = rna_indices_shared + rna_indices_only
    atac_indices_final = atac_indices_shared + atac_indices_only
    
    print(f"\nFinal: {len(rna_indices_final)} RNA cells, {len(atac_indices_final)} ATAC cells")
    
    # Subset (handle both backed and in-memory modes)
    rna_aligned = rna[rna_indices_final]
    atac_aligned = atac[atac_indices_final]
    
    # Convert to memory if backed
    if rna_aligned.isbacked:
        rna_aligned = rna_aligned.to_memory()
    if atac_aligned.isbacked:
        atac_aligned = atac_aligned.to_memory()
    
    return rna_aligned, atac_aligned


def create_mudata(rna: ad.AnnData, atac: ad.AnnData) -> md.MuData:
    """Create MuData object with merged obs."""
    mdata = md.MuData({"rna": rna, "atac": atac})
    
    # Add unified batch column
    mdata.obs["batch"] = mdata.obs["rna:batch"].combine_first(mdata.obs["atac:batch"])
    
    return mdata


def main() -> None:
    args = parse_args()
    
    rna_path = Path(args.rna)
    atac_path = Path(args.atac)
    output_path = Path(args.output)
    
    print(f"Loading RNA: {rna_path}")
    rna = ad.read_h5ad(rna_path)
    print(f"  Shape: {rna.shape}")
    
    print(f"\nLoading ATAC: {atac_path}")
    atac = ad.read_h5ad(atac_path)
    print(f"  Shape: {atac.shape}")
    
    print("\nAligning by batch...")
    rna_aligned, atac_aligned = align_rna_atac(rna, atac)
    
    print(f"\nRNA aligned: {rna_aligned.shape}")
    print(f"ATAC aligned: {atac_aligned.shape}")
    
    print("\nCreating MuData object...")
    mdata = create_mudata(rna_aligned, atac_aligned)
    
    print(f"\nMuData shape: {mdata.shape}")
    print(f"MuData obs columns: {list(mdata.obs.columns)}")
    
    output_path.parent.mkdir(parents=True, exist_ok=True)
    mdata.write(output_path)
    
    print(f"\nSaved: {output_path}")


if __name__ == "__main__":
    main()
