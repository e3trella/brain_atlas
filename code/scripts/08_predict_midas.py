#!/usr/bin/env python
"""Stream MIDAS predictions for RNA+ATAC MuData to disk.

Expected inputs:
- A `.h5mu` file containing `rna` and `atac` modalities.
- A model directory produced by `model.save(...)`, containing `model.pt` and `setup.json`.
- A TOML feature-dimension file, typically with keys like:

    rna = [4000]
    atac = [2897, 2007, 3150, 1988]

The script writes only:
- joint latent representation (`z/joint`, later loadable as `z_c`/`z_u`)
- batch-corrected reconstructed matrices (`x_bc/rna`, `x_bc/atac`)

Outputs are streamed as mini-batch `.npy` shards by default.
"""
from __future__ import annotations

import argparse
import json
import shutil
from pathlib import Path
from typing import Any

import mudata as mu
import toml

import scmidas


REQUIRED_MODALITIES = ("rna", "atac")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Load a saved MIDAS model and stream joint latent + batch-corrected "
            "RNA/ATAC predictions for a .h5mu MuData file to disk."
        )
    )
    parser.add_argument(
        "--mdata",
        required=True,
        type=Path,
        help="Path to input .h5mu MuData file.",
    )
    parser.add_argument(
        "--model-dir",
        required=True,
        type=Path,
        help="Directory produced by model.save(...), containing model.pt and setup.json.",
    )
    parser.add_argument(
        "--feat-dims",
        required=True,
        type=Path,
        help="TOML file containing RNA/ATAC feature chunks, e.g. rna=[...], atac=[...].",
    )
    parser.add_argument(
        "--out-dir",
        required=True,
        type=Path,
        help="Directory where streamed prediction shards will be written.",
    )
    parser.add_argument(
        "--batch-key",
        default="batch",
        help="Batch column present in each modality's .obs. Default: batch.",
    )
    parser.add_argument(
        "--batch-size",
        default=256,
        type=int,
        help="Prediction DataLoader batch size used by MIDAS. Default: 256.",
    )
    parser.add_argument(
        "--save-format",
        default="npy",
        choices=("npy", "csv"),
        help="Shard format written by model.predict(...). Default: npy.",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Remove --out-dir first if it already exists.",
    )
    parser.add_argument(
        "--quiet",
        action="store_true",
        help="Disable MIDAS prediction progress bars/log messages where supported.",
    )
    return parser.parse_args()


def load_feat_dims(path: Path, mdata: mu.MuData) -> dict[str, list[int]]:
    raw: dict[str, Any] = toml.load(path)
    dims_x: dict[str, list[int]] = {}

    for modality in REQUIRED_MODALITIES:
        if modality in raw:
            value = raw[modality]
            if not isinstance(value, list) or not value:
                raise ValueError(f"{path}: {modality!r} must be a non-empty list of integers.")
            dims_x[modality] = [int(v) for v in value]
        else:
            dims_x[modality] = [int(mdata[modality].n_vars)]

    for modality, chunks in dims_x.items():
        expected = int(mdata[modality].n_vars)
        observed = int(sum(chunks))
        if observed != expected:
            raise ValueError(
                f"Feature dimension mismatch for {modality!r}: "
                f"sum(feat_dims)={observed}, but mdata[{modality!r}].n_vars={expected}."
            )

    return dims_x


def validate_inputs(args: argparse.Namespace) -> None:
    if not args.mdata.exists():
        raise FileNotFoundError(f"Input MuData not found: {args.mdata}")
    if not args.model_dir.is_dir():
        raise FileNotFoundError(f"Model directory not found: {args.model_dir}")
    if not (args.model_dir / "model.pt").exists():
        raise FileNotFoundError(f"Missing model file: {args.model_dir / 'model.pt'}")
    if not (args.model_dir / "setup.json").exists():
        raise FileNotFoundError(f"Missing setup file: {args.model_dir / 'setup.json'}")
    if not args.feat_dims.exists():
        raise FileNotFoundError(f"Feature-dimension TOML not found: {args.feat_dims}")
    if args.out_dir.exists() and any(args.out_dir.iterdir()):
        if not args.overwrite:
            raise FileExistsError(
                f"Output directory already exists and is non-empty: {args.out_dir}. "
                "Pass --overwrite to remove it first."
            )
        shutil.rmtree(args.out_dir)


def main() -> None:
    args = parse_args()
    validate_inputs(args)

    print(f"Loading MuData: {args.mdata}")
    mdata = mu.read_h5mu(str(args.mdata))

    missing_modalities = [m for m in REQUIRED_MODALITIES if m not in mdata.mod]
    if missing_modalities:
        raise ValueError(
            f"Input MuData is missing required modalities: {missing_modalities}. "
            f"Available modalities: {list(mdata.mod.keys())}"
        )

    missing_batch = [m for m in REQUIRED_MODALITIES if args.batch_key not in mdata[m].obs.columns]
    if missing_batch:
        raise ValueError(
            f"batch_key={args.batch_key!r} missing from .obs of modalities {missing_batch}."
        )

    dims_x = load_feat_dims(args.feat_dims, mdata)
    print(f"Registering MuData with dims_x={dims_x}")
    scmidas.MIDAS.setup_mudata(mdata, batch_key=args.batch_key, dims_x=dims_x)

    print(f"Loading MIDAS model: {args.model_dir}")
    model = scmidas.MIDAS.load(
        str(args.model_dir),
        mdata,
        batch_size=args.batch_size,
        save_model_path=str(args.model_dir),
    )

    print(f"Streaming predictions to: {args.out_dir}")
    manifest = model.predict(
        return_in_memory=False,
        save_dir=str(args.out_dir),
        save_format=args.save_format,
        joint_latent=True,
        mod_latent=False,
        impute=False,
        batch_correct=True,
        translate=False,
        input=False,
        verbose=not args.quiet,
    )

    args.out_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = args.out_dir / "predict_manifest.json"
    with manifest_path.open("w", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2)

    print("Done.")
    print(f"Manifest: {manifest_path}")
    print("Expected output groups include:")
    print("  - <out-dir>/<batch>/z/joint/*.npy")
    print("  - <out-dir>/<batch>/x_bc/rna/*.npy")
    print("  - <out-dir>/<batch>/x_bc/atac/*.npy")


if __name__ == "__main__":
    main()