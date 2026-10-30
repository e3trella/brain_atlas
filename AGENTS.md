# Agent Instructions

## Scope and Style

This repository supports single-cell brain atlas data integration and processing workflows.

Keep code direct, concise, and readable. Add defensive checks, wrappers, abstractions, or configuration layers only when they solve a specific problem in the current task.

## Environment

- The repository lives on a remote compute cluster. The local working directory is an sshfs mount of it.
- Use local file tools only to read and edit code and text files, and to view image files in `results/`, such as PNG or PDF plots. Run every command on the cluster through `./remote.sh`. Never call `ssh`, `scp`, or `rsync` directly, and never run `python`, `Rscript`, `git`, or other programs locally inside the mount.
  - Short commands: `./remote.sh run '<command>'`
  - Long jobs: `./remote.sh bg <name> '<command>'`, logged to `results/logs/<name>.log`
  - Job progress: `./remote.sh status <name>`. Do not restart a job that is still running.
- Commands on the cluster run as root inside a Docker container. Keep all operations inside the repository and never modify the container environment or files outside the repository.
- Set thread counts explicitly in every job, for example `OMP_NUM_THREADS`, `n_jobs`, or `future::plan()`. Use at most 16 cores.
- The user also edits code directly on the cluster, especially plotting code. Re-read a file immediately before editing it.
- Use `./remote.sh bg` for any command expected to run longer than about two minutes, including anything that loads a full dataset.

## Repository Layout

- `code/scripts/`: runnable, reproducible workflows.
- `code/src/r/` and `code/src/py/`: reusable R and Python helpers, free of project-specific hard-coded paths.
- `code/notebooks/`: exploratory analysis, QC checks, plot drafts, and temporary scripts. Move stable or reusable logic into `code/scripts/` or `code/src/`.
- `code/deprecated/`: historical code. Do not extend it unless explicitly restoring an old workflow.
- `parameters/`: path lists, reference files, and run configuration. Scripts read input paths from here instead of hard-coding them.
- `data/raw/`: read-only source data.
- `data/processed/`: reusable processed outputs.
- `results/`: final tables, plots, and logs.

## Data

- Never open files under `data/` or large binary files anywhere in the repository with local tools. This includes `*.h5ad`, `*.h5`, `*.rds`, `*.qs`, `*.loom`, `*.mtx`, `*.tsv.gz`, and `*fragments*`.
- To learn the structure of a dataset, run a command on the cluster that prints only a summary, such as dimensions, metadata columns, or label values. Never print expression matrices or full tables.
- Never write to `data/raw/`. Scripts may write and overwrite their own outputs in `data/processed/` and `results/`. Do not delete or move existing files there, or overwrite files produced by other scripts, without explicit approval.

## R

- Run R scripts with `Rscript`, for example `./remote.sh run 'Rscript code/scripts/example.R'`.
- Use Seurat v5 APIs, for example `Seurat::GetAssayData(layer = "counts")` rather than the deprecated `slot = "counts"`.
- When working with Seurat objects, complete the task from `obj@meta.data` where possible, for example filtering samples there, instead of subsetting or rebuilding the object.

## Python

- Run Python only through `uv run --no-sync`, for example `uv run --no-sync python code/scripts/example.py`. The `--no-sync` flag prevents uv from modifying or rebuilding the environment.
- The project `.venv` is a symbolic link to an environment stored elsewhere. Never delete, recreate, replace, or relink `.venv` or its target.

## Dependencies

- Never install anything yourself, including R packages, Python packages, system packages, and command-line tools.
- If a dependency is missing, report it and give the exact installation command for the user to run, using `pak` for R and `uv` for Python. Continue only after the user has installed it.

## Git

- Run only read-only git commands through `./remote.sh`, such as `git status`, `git diff`, and `git log`. Never run git commands that change the working tree, the index, or history.
