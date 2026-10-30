#!/usr/bin/env python
"""Convert raw MIDAS `.pt` checkpoints into `model.save(...)` directories.

This script is intended for training directories that contain checkpoint files such as
`model_epoch100_20260605-153526.pt`. For each checkpoint, it instantiates a MIDAS
model against the original MuData object, loads the checkpoint with
`model.load_checkpoint(...)`, and writes a sibling subdirectory via `model.save(...)`.

By default, existing non-empty output directories are skipped. Pass `--overwrite` to
replace outputs through the MIDAS save API.
"""
from __future__ import annotations

import argparse
import logging
import re
import tomllib
from pathlib import Path

import mudata as mu

import scmidas
from scmidas.config import load_config


REQUIRED_MODALITIES = ("rna", "atac")
CHECKPOINT_RE = re.compile(r"^model_epoch(?P<epoch>\d+)_.+\.pt$")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Load MIDAS .pt checkpoints and save each one as a model.save(...) "
            "directory containing model.pt and setup.json."
        )
    )
    parser.add_argument(
        "--mdata",
        required=True,
        type=Path,
        help="MuData file used to train the checkpoints.",
    )
    parser.add_argument(
        "--checkpoint-dir",
        required=True,
        type=Path,
        help="Directory containing raw .pt checkpoints.",
    )
    parser.add_argument(
        "--feat-dims",
        required=True,
        type=Path,
        help="TOML file containing RNA/ATAC feature chunks used for training.",
    )
    parser.add_argument(
        "--checkpoint",
        type=Path,
        action="append",
        default=None,
        help=(
            "Specific .pt checkpoint to convert. May be supplied multiple times. "
            "If omitted, all model_epoch*.pt files in --checkpoint-dir are converted."
        ),
    )
    parser.add_argument(
        "--batch-key",
        default="batch",
        help="Batch column used for MIDAS setup. Default: batch.",
    )
    parser.add_argument(
        "--batch-size",
        type=int,
        default=512,
        help="Batch size used to instantiate MIDAS. Default: 512.",
    )
    parser.add_argument(
        "--num-workers",
        type=int,
        default=4,
        help="MIDAS config num_workers override. Default: 4.",
    )
    parser.add_argument(
        "--n-max",
        type=int,
        default=3500,
        help="MIDAS config n_max override. Default: 3500.",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=None,
        help="Convert only the first N checkpoints after sorting by epoch.",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite existing model.save(...) directories instead of skipping them.",
    )
    return parser.parse_args()


def checkpoint_sort_key(path: Path) -> tuple[int, str]:
    match = CHECKPOINT_RE.match(path.name)
    if match is None:
        return (10**12, path.name)
    return (int(match.group("epoch")), path.name)


def resolve_checkpoint_path(checkpoint: Path, checkpoint_dir: Path) -> Path:
    if checkpoint.is_absolute():
        return checkpoint
    if checkpoint.exists():
        return checkpoint
    return checkpoint_dir / checkpoint


def find_checkpoints(
    checkpoint_dir: Path,
    limit: int | None,
    selected_checkpoints: list[Path] | None = None,
) -> list[Path]:
    if selected_checkpoints:
        checkpoints = [
            resolve_checkpoint_path(checkpoint, checkpoint_dir)
            for checkpoint in selected_checkpoints
        ]
    else:
        checkpoints = sorted(
            [path for path in checkpoint_dir.glob("*.pt") if CHECKPOINT_RE.match(path.name)],
            key=checkpoint_sort_key,
        )

    if limit is not None:
        checkpoints = checkpoints[:limit]
    return checkpoints


def validate_inputs(args: argparse.Namespace) -> None:
    if not args.mdata.exists():
        raise FileNotFoundError(f"Input MuData not found: {args.mdata}")
    if not args.checkpoint_dir.is_dir():
        raise FileNotFoundError(f"Checkpoint directory not found: {args.checkpoint_dir}")
    if not args.feat_dims.is_file():
        raise FileNotFoundError(f"Feature-dimension TOML not found: {args.feat_dims}")
    if args.limit is not None and args.limit < 1:
        raise ValueError("--limit must be a positive integer when provided.")
    if args.checkpoint:
        for checkpoint in args.checkpoint:
            checkpoint_path = resolve_checkpoint_path(checkpoint, args.checkpoint_dir)
            if not checkpoint_path.is_file():
                raise FileNotFoundError(f"Checkpoint file not found: {checkpoint_path}")
            if checkpoint_path.suffix != ".pt":
                raise ValueError(f"Checkpoint path must end in .pt: {checkpoint_path}")


def validate_feat_dims(mdata: mu.MuData, feat_dims: dict[str, list[int]]) -> None:
    missing_modalities = [modality for modality in feat_dims if modality not in mdata.mod]
    if missing_modalities:
        raise ValueError(
            f"Input MuData is missing required modalities: {missing_modalities}. "
            f"Available modalities: {list(mdata.mod.keys())}"
        )

    for modality, chunks in feat_dims.items():
        expected = int(mdata[modality].n_vars)
        observed = int(sum(chunks))
        if observed != expected:
            raise ValueError(
                f"Feature dimension mismatch for {modality!r}: "
                f"sum(feat_dims)={observed}, but mdata[{modality!r}].n_vars={expected}."
            )


def load_feat_dims(path: Path) -> dict[str, list[int]]:
    with path.open("rb") as file:
        raw = tomllib.load(file)

    return {
        modality: [int(value) for value in raw[modality]]
        for modality in REQUIRED_MODALITIES
    }


def main() -> None:
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    args = parse_args()
    validate_inputs(args)

    checkpoints = find_checkpoints(args.checkpoint_dir, args.limit, args.checkpoint)
    if not checkpoints:
        raise FileNotFoundError(f"No model_epoch*.pt checkpoints found in {args.checkpoint_dir}")
    logging.info("Found %d checkpoint(s)", len(checkpoints))

    logging.info("Reading MuData: %s", args.mdata)
    mdata = mu.read_h5mu(str(args.mdata))

    feat_dims = load_feat_dims(args.feat_dims)
    validate_feat_dims(mdata, feat_dims)

    scmidas.MIDAS.setup_mudata(mdata, batch_key=args.batch_key, dims_x=feat_dims)

    configs = load_config()
    configs["num_workers"] = args.num_workers
    configs["n_max"] = args.n_max
    model = scmidas.MIDAS(
        mdata,
        configs=configs,
        save_model_path=str(args.checkpoint_dir),
        batch_size=args.batch_size,
    )

    total = len(checkpoints)
    for index, checkpoint_path in enumerate(checkpoints, start=1):
        output_dir = checkpoint_path.with_suffix("")
        if output_dir.exists() and any(output_dir.iterdir()) and not args.overwrite:
            logging.info(
                "[Checkpoint %d/%d] skipped existing output: %s",
                index,
                total,
                output_dir,
            )
            continue

        logging.info("[Checkpoint %d/%d] loading: %s", index, total, checkpoint_path)
        model.load_checkpoint(str(checkpoint_path))

        logging.info("[Checkpoint %d/%d] saving: %s", index, total, output_dir)
        model.save(str(output_dir), overwrite=args.overwrite)

    logging.info("All checkpoint conversions completed")


if __name__ == "__main__":
    main()
