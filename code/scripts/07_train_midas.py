import warnings
warnings.filterwarnings("ignore")

import logging
logging.basicConfig(level=logging.INFO)

import os
import subprocess
import tomllib
from pathlib import Path

import lightning as L
from lightning.pytorch import loggers as pl_loggers

import numpy as np
import pandas as pd
import scanpy as sc
import mudata as mu

import scmidas
from scmidas.config import load_config
from scmidas.data import download_data, download_models, download_script


# =============================================================================
# Adjustable parameters
# =============================================================================

# Reproducibility and plotting
SEED = 42
FIGURE_SIZE = (4, 4)

# Task and paths
TASK = "Atlas_20260721_imm"
DATA_PATH = "data/processed/Atlas_20260721/merged_imm.h5mu"
SAVE_MODEL_PATH = f"data/processed/Atlas_20260721/saved_models/{TASK}"
LOG_DIR = "./logs/"
PRETRAINED_CHECKPOINT_PATH = f"data/processed/Atlas_20260721/saved_models/{TASK}.pt"

# GPU / Lightning configuration
# Single-GPU default:
GPUS = "0"
STRATEGY = "auto"
ACCELERATOR = "gpu"

# Multi-GPU options:
# - In a notebook/script-like interactive run, DDP startup can be slower:
#       GPUS = "0,1"
#       STRATEGY = "ddp_notebook"
# - As a normal Python script, recommended for production multi-GPU:
#       GPUS = "0,1"
#       STRATEGY = "ddp"

# Data loading and training
BATCH_KEY = "batch"
BATCH_SIZE = 256
NUM_WORKERS = 4
MAX_EPOCHS = 4000
N_SAVE = 100
N_MAX = 500 # max number of cells per batch
SAMPLER_TYPE = "auto"
VIZ_UMAP_TB = True
USE_PRETRAINED = False  # To train from scratch, set to False.

# Feature dimensions used by MIDAS setup.
# RNA has one feature block; ATAC is split into chromosome/chunk-sized blocks.
FEAT_DIMS_PATH = Path(__file__).resolve().parents[2] / "parameters/20260721-feat_dim.toml"
with FEAT_DIMS_PATH.open("rb") as file:
    FEAT_DIMS = tomllib.load(file)


# =============================================================================
# Runtime setup
# =============================================================================

if GPUS is not None:
    os.environ["CUDA_VISIBLE_DEVICES"] = GPUS

# Lightning's `devices="auto"` can default to 1 in notebooks even when more GPUs
# are visible, so derive an explicit device count from `GPUS`.
DEVICES = len(GPUS.split(",")) if GPUS else "auto"

sc.set_figure_params(figsize=FIGURE_SIZE)
L.seed_everything(SEED)

run_params = {
    "task": TASK,
    "data_path": DATA_PATH,
    "feat_dims_path": FEAT_DIMS_PATH,
    "save_model_path": SAVE_MODEL_PATH,
    "log_dir": LOG_DIR,
    "pretrained_checkpoint_path": PRETRAINED_CHECKPOINT_PATH,
    "gpus": GPUS,
    "devices": DEVICES,
    "strategy": STRATEGY,
    "accelerator": ACCELERATOR,
    "batch_key": BATCH_KEY,
    "num_workers": NUM_WORKERS,
    "batch_size": BATCH_SIZE,
    "n_save": N_SAVE,
    "sampler_type": SAMPLER_TYPE,
    "viz_umap_tb": VIZ_UMAP_TB,
    "max_epochs": MAX_EPOCHS,
    "use_pretrained": USE_PRETRAINED,
    "n_max": N_MAX,
    "seed": SEED,
}
logging.info("MIDAS run parameters:")
for key, value in run_params.items():
    logging.info("  %s: %s", key, value)


# =============================================================================
# Load data and configure model
# =============================================================================

mdata = mu.read_h5mu(DATA_PATH)
mdata.uns["feat_dims"] = FEAT_DIMS

configs = load_config()
configs["num_workers"] = NUM_WORKERS
configs["n_max"] = N_MAX
logging.info("Loaded MIDAS config overrides:")
logging.info("  num_workers: %s", configs["num_workers"])
logging.info("  n_max: %s", configs["n_max"])

# `dims_x` carries the ATAC chromosome chunk sizes that `from_dir` read from
# feat_dims.toml. ATAC counts are binarized by default.
scmidas.MIDAS.setup_mudata(mdata, batch_key=BATCH_KEY, dims_x=mdata.uns["feat_dims"])

model = scmidas.MIDAS(
    mdata,
    configs=configs,
    save_model_path=SAVE_MODEL_PATH,
    batch_size=BATCH_SIZE,
    n_save=N_SAVE,
    sampler_type=SAMPLER_TYPE,
    viz_umap_tb=VIZ_UMAP_TB,
)

tb_logger = pl_loggers.TensorBoardLogger(
    save_dir=LOG_DIR,
    version=TASK,
)

if USE_PRETRAINED:
    model.load_checkpoint(PRETRAINED_CHECKPOINT_PATH)
else:
    model.train(
        max_epochs=MAX_EPOCHS,
        accelerator=ACCELERATOR,
        devices=DEVICES,
        strategy=STRATEGY,
        logger=tb_logger,
    )
