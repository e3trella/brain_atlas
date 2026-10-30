#!/usr/bin/env python
"""Convert streamed MIDAS prediction shards into a single-modality AnnData file.

This script reads a prediction directory produced by
`model.predict(..., save_dir=...)`, such as the output of
`scripts/predict_midas_rna_atac.py`, and exports one selected modality
(`rna` or `atac`) as a standalone `.h5ad` file.

The exported AnnData contains:
- `adata.X`: batch-corrected reconstructed matrix from `x_bc[modality]`
- `adata.obsm['X_midas']`: MIDAS joint biological latent `z_c['joint']`
- `adata.obsm['X_midas_u']`: MIDAS joint technical latent `z_u['joint']`

If `--mdata` is provided, cell names and feature names are copied from the
original MuData. This is strongly recommended. The file is opened with
`backed='r'` because this script only reads metadata and does not touch `X`.
"""
from __future__ import annotations

import argparse
from datetime import datetime
from pathlib import Path
from typing import Any
import re

import anndata as ad
import matplotlib.pyplot as plt
import mudata as mu
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
import torch


def log_progress(message: str) -> None:
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{timestamp}] {message}", flush=True)


def check_dependencies() -> None:
    """Verify optional runtime dependencies before the long run starts.

    The packages imported at the top of this module are hard requirements and
    fail at import time. This function covers dependencies that are only needed
    later in the pipeline, so a missing one is reported immediately instead of
    after the shard concatenation has already run.
    """
    try:
        import igraph  # noqa: F401
    except ImportError as exc:
        raise ImportError(
            "igraph is required for sc.tl.leiden(flavor='igraph'). "
            "Install it with: uv add igraph"
        ) from exc


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Convert MIDAS streamed prediction shards for one modality to .h5ad."
    )
    parser.add_argument(
        "--pred-dir",
        required=True,
        type=Path,
        help="Prediction directory produced by model.predict(..., save_dir=...).",
    )
    parser.add_argument(
        "--modality",
        required=True,
        choices=("rna", "atac"),
        help="Which modality to export as AnnData.",
    )
    parser.add_argument(
        "--output",
        required=True,
        type=Path,
        help="Output .h5ad file path.",
    )
    parser.add_argument(
        "--mdata",
        type=Path,
        default=None,
        help="Optional original .h5mu file. Recommended to preserve obs_names/var_names.",
    )
    parser.add_argument(
        "--model-dir",
        type=Path,
        default=None,
        help="Optional MIDAS model.save(...) directory. Used to infer dim_c from model.pt configs.",
    )
    parser.add_argument(
        "--dim-c",
        type=int,
        default=None,
        help="Biological latent dimension. Required unless --model-dir can infer it.",
    )
    parser.add_argument(
        "--save-format",
        default="npy",
        choices=("npy", "csv"),
        help="Shard format used by prediction. Default: npy.",
    )
    parser.add_argument(
        "--batch-key",
        default="batch",
        help="Batch column in original mdata. Default: batch.",
    )
    return parser.parse_args()


def infer_dim_c(model_dir: Path | None, dim_c: int | None) -> int:
    if dim_c is not None:
        return dim_c
    if model_dir is not None:
        model_pt = model_dir / "model.pt"
        if not model_pt.exists():
            raise FileNotFoundError(f"Missing model file: {model_pt}")
        state: dict[str, Any] = torch.load(str(model_pt), map_location="cpu", weights_only=False)
        configs = state.get("configs") or {}
        if "dim_c" in configs:
            return int(configs["dim_c"])
    raise ValueError("Please provide --dim-c, or provide --model-dir so dim_c can be inferred.")


_NUM_SUFFIX = re.compile(r"(\d+)(?=\.[^.]+$)")


def _extract_index(path: Path) -> int | None:
    match = _NUM_SUFFIX.search(path.name)
    if match is None:
        return None
    return int(match.group(1))


def _chunk_sort_key(path: Path) -> tuple[bool, int, str]:
    index = _extract_index(path)
    return (index is None, index or 0, path.name)


def _load_array_file(path: Path, save_format: str) -> np.ndarray:
    if save_format == "npy":
        return np.load(path, allow_pickle=False)
    if save_format == "csv":
        return np.loadtxt(path, delimiter=",")
    raise ValueError(f"Unsupported save_format={save_format!r}.")


def _load_chunked_array(key_dir: Path, save_format: str) -> np.ndarray:
    if not key_dir.is_dir():
        raise FileNotFoundError(f"Missing prediction shard directory: {key_dir}")

    chunk_paths = sorted(
        [path for path in key_dir.iterdir() if path.is_file() and path.name != "meta.npy"],
        key=_chunk_sort_key,
    )
    if not chunk_paths:
        raise FileNotFoundError(f"No prediction shard files found in: {key_dir}")

    log_progress(f"Loading {len(chunk_paths)} shard(s): {key_dir}")
    chunks = []
    for chunk_index, path in enumerate(chunk_paths, start=1):
        array = _load_array_file(path, save_format)
        if array.ndim == 1:
            array = array.reshape(-1, 1)
        chunks.append(array)
        if chunk_index == 1 or chunk_index == len(chunk_paths) or chunk_index % 10 == 0:
            log_progress(
                f"Loaded shard {chunk_index}/{len(chunk_paths)} from {key_dir.name}: "
                f"shape={array.shape}"
            )

    log_progress(f"Concatenating shards for {key_dir}")
    result = np.concatenate(chunks, axis=0)
    log_progress(f"Finished {key_dir}: shape={result.shape}, dtype={result.dtype}")
    return result


def _list_batches(pred_dir: Path) -> list[str]:
    if not pred_dir.is_dir():
        raise FileNotFoundError(f"Prediction directory not found: {pred_dir}")
    return sorted(path.name for path in pred_dir.iterdir() if path.is_dir())


def _collect_batch_obs(
    batch_name: str,
    n_cells: int,
    modality: str,
    original: mu.MuData | None,
    batch_key: str,
    cell_offset: int,
) -> pd.DataFrame:
    if original is None:
        obs_names = [f"cell_{i}" for i in range(cell_offset, cell_offset + n_cells)]
        return pd.DataFrame({batch_key: [batch_name] * n_cells}, index=obs_names)

    if modality not in original.mod:
        raise KeyError(f"Modality {modality!r} not found in original MuData.")
    if batch_key not in original.obs.columns:
        raise KeyError(f"batch_key={batch_key!r} missing from original mdata.obs.")

    mask = original.obs[batch_key].astype(str) == str(batch_name)
    obs_names = original.obs_names[mask].tolist()
    if len(obs_names) != n_cells:
        raise RuntimeError(
            f"Cell-name mismatch for batch {batch_name!r}: original mdata.obs has "
            f"{len(obs_names)} cells, but predictions have {n_cells} rows."
        )

    return pd.DataFrame({batch_key: [batch_name] * n_cells}, index=obs_names)


def collect_prediction_parts(
    pred_dir: Path,
    modality: str,
    save_format: str,
    dim_c: int,
    original: mu.MuData | None,
    batch_key: str,
) -> tuple[list[sp.csr_matrix], list[np.ndarray], list[np.ndarray], pd.DataFrame]:
    matrix_parts: list[sp.csr_matrix] = []
    z_c_parts: list[np.ndarray] = []
    z_u_parts: list[np.ndarray] = []
    obs_parts: list[pd.DataFrame] = []
    cell_offset = 0

    batch_names = _list_batches(pred_dir)
    total_batches = len(batch_names)

    for batch_index, batch_name in enumerate(batch_names, start=1):
        log_progress(f"Reading batch {batch_index}/{total_batches}: {batch_name}")
        batch_dir = pred_dir / batch_name
        z = _load_chunked_array(batch_dir / "z" / "joint", save_format)
        x_bc_dense = _load_chunked_array(batch_dir / "x_bc" / modality, save_format)

        if z.shape[0] != x_bc_dense.shape[0]:
            raise RuntimeError(
                f"Row mismatch for batch {batch_name!r}: z/joint has "
                f"{z.shape[0]} rows but x_bc/{modality} has {x_bc_dense.shape[0]} rows."
            )
        if z.shape[1] < dim_c:
            raise RuntimeError(
                f"Latent dimension mismatch for batch {batch_name!r}: z/joint has "
                f"{z.shape[1]} columns but dim_c={dim_c}."
            )

        n_cells = x_bc_dense.shape[0]
        log_progress(
            f"Collecting obs for batch {batch_index}/{total_batches}: "
            f"{batch_name}, n_cells={n_cells}"
        )
        obs_parts.append(
            _collect_batch_obs(
                batch_name=batch_name,
                n_cells=n_cells,
                modality=modality,
                original=original,
                batch_key=batch_key,
                cell_offset=cell_offset,
            )
        )
        log_progress(
            f"Converting batch {batch_index}/{total_batches} modality matrix to CSR: "
            f"{batch_name}, shape={x_bc_dense.shape}, dtype={x_bc_dense.dtype}"
        )
        x_bc = sp.csr_matrix(x_bc_dense)
        del x_bc_dense
        log_progress(
            f"Finished batch CSR conversion: {batch_name}, "
            f"shape={x_bc.shape}, nnz={x_bc.nnz}, dtype={x_bc.dtype}"
        )
        matrix_parts.append(x_bc)
        z_c_parts.append(z[:, :dim_c])
        z_u_parts.append(z[:, dim_c:])
        cell_offset += n_cells
        log_progress(
            f"Finished batch {batch_index}/{total_batches}: {batch_name}; "
            f"cumulative_cells={cell_offset}"
        )

    if not matrix_parts:
        raise FileNotFoundError(f"No batch prediction directories found in: {pred_dir}")

    log_progress(f"Concatenating obs from {len(obs_parts)} batch(es)")
    obs = pd.concat(obs_parts, axis=0)
    obs[batch_key] = pd.Categorical(obs[batch_key])
    log_progress(f"Finished obs: shape={obs.shape}")
    return matrix_parts, z_c_parts, z_u_parts, obs


def build_adata(
    matrix_parts: list[sp.csr_matrix],
    z_c_parts: list[np.ndarray],
    z_u_parts: list[np.ndarray],
    modality: str,
    obs: pd.DataFrame,
    original: mu.MuData | None,
) -> ad.AnnData:
    log_progress(
        f"Combining CSR modality matrix from {len(matrix_parts)} batch(es); "
        f"part_shapes={[part.shape for part in matrix_parts]}"
    )
    matrix = sp.vstack(matrix_parts, format="csr")
    log_progress(
        f"Finished combining CSR modality matrix: "
        f"shape={matrix.shape}, nnz={matrix.nnz}, dtype={matrix.dtype}"
    )

    log_progress(f"Stacking biological latent matrix from {len(z_c_parts)} batch(es)")
    z_c = np.vstack(z_c_parts)
    log_progress(f"Finished stacking biological latent matrix: shape={z_c.shape}, dtype={z_c.dtype}")

    log_progress(f"Stacking technical latent matrix from {len(z_u_parts)} batch(es)")
    z_u = np.vstack(z_u_parts)
    log_progress(f"Finished stacking technical latent matrix: shape={z_u.shape}, dtype={z_u.dtype}")

    if matrix.shape[0] != obs.shape[0]:
        raise RuntimeError(
            f"Row mismatch: matrix has {matrix.shape[0]} cells but obs has {obs.shape[0]} rows."
        )
    if z_c.shape[0] != matrix.shape[0] or z_u.shape[0] != matrix.shape[0]:
        raise RuntimeError("Latent rows do not match modality matrix rows.")

    var = None
    if original is not None:
        log_progress(f"Copying var metadata for modality={modality}")
        var = original[modality].var.copy()
        log_progress(f"Finished copying var metadata: shape={var.shape}")

    log_progress("Constructing AnnData object")
    adata = ad.AnnData(X=matrix, obs=obs.copy(), var=var)
    adata.obsm["X_midas"] = z_c
    adata.obsm["X_midas_u"] = z_u
    adata.uns["midas_modality"] = modality
    adata.uns["midas_source"] = "x_bc"
    log_progress(f"Finished AnnData object: shape={adata.shape}")
    return adata


def compute_midas_graph_and_umap(adata: ad.AnnData, plot_output: Path) -> None:
    log_progress("[Neighbors] started")
    sc.pp.neighbors(
        adata,
        n_neighbors=30,
        use_rep="X_midas",
        key_added="midas_neighbors",
        random_state=8731,
    )
    log_progress("[Neighbors] completed")

    log_progress("[Leiden clustering] started")
    sc.tl.leiden(
        adata,
        neighbors_key="midas_neighbors",
        key_added="midas_leiden",
        resolution=1.0,
        random_state=8731,
        flavor="igraph",
        n_iterations=2,
        directed=False
    )
    log_progress("[Leiden clustering] completed")

    log_progress("[UMAP] started")
    sc.tl.umap(
        adata,
        neighbors_key="midas_neighbors",
        random_state=8731,
    )
    adata.obsm["X_midas_umap"] = adata.obsm["X_umap"].copy()
    log_progress("[UMAP] completed")

    log_progress(f"[UMAP plot] saving to: {plot_output}")
    sc.pl.embedding(
        adata,
        basis="midas_umap",
        color=None,
        show=False,
    )
    plt.savefig(plot_output, dpi=300, bbox_inches="tight")
    plt.close()
    log_progress(f"[UMAP plot] saved: {plot_output}")


def main() -> None:
    args = parse_args()
    check_dependencies()
    if args.output.suffix.lower() != ".h5ad":
        raise ValueError("--output must end with .h5ad")

    dim_c = infer_dim_c(args.model_dir, args.dim_c)
    log_progress(f"Using dim_c={dim_c}")
    if args.mdata is not None:
        log_progress(f"Opening original MuData in backed read mode: {args.mdata}")
        original = mu.read_h5mu(str(args.mdata), backed="r")
        log_progress(f"Opened original MuData: n_obs={original.n_obs}, modalities={list(original.mod.keys())}")
    else:
        original = None

    matrix_parts, z_c_parts, z_u_parts, obs = collect_prediction_parts(
        args.pred_dir,
        modality=args.modality,
        save_format=args.save_format,
        dim_c=dim_c,
        original=original,
        batch_key=args.batch_key,
    )
    adata = build_adata(
        matrix_parts,
        z_c_parts,
        z_u_parts,
        args.modality,
        obs,
        original,
    )

    plot_output = args.output.with_name(f"{args.output.stem}_midas_umap.png")
    compute_midas_graph_and_umap(adata, plot_output)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    log_progress(f"Writing AnnData to: {args.output}")
    adata.write_h5ad(args.output)
    log_progress(f"Wrote: {args.output}")


if __name__ == "__main__":
    main()
