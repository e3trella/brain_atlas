#!/usr/bin/env Rscript

# Convert all RNA Seurat RDS files listed in parameters/rna_path.txt to
# reference-aligned 10x-style RNA matrices using align_rna_features().
# Default is dry-run. Add --run to write output files.

options(stringsAsFactors = FALSE)

usage <- function() {
  cat(
    paste(
      "用法:",
      "  Rscript code/scripts/02_align_all_rna_features.R [--input-list PATH]",
      "                                                   [--ref-genes PATH]",
      "                                                   [--outputdir PATH]",
      "                                                   [--assay NAME]",
      "                                                   [--limit N]",
      "                                                   [--run]",
      "",
      "说明:",
      "  1. 只处理路径列表中以 rna.rds 结尾的文件，忽略 atac.h5seurat。",
      "  2. 输出目录为 --outputdir/DONOR_ID/（默认 data/processed/Atlas_20260721/rna）。",
      "  3. 默认 dry-run，只打印计划，不读取 RDS、不写文件；加 --run 才执行转换。",
      "  4. align_rna_features() 会写出 matrix.mtx、barcodes.tsv.gz、",
      "     features.tsv.gz、rna_ref_presence.csv.gz，随后自动将 matrix.mtx 压缩为 matrix.mtx.gz。",
      "",
      "示例:",
      "  Rscript code/scripts/align_all_rna_features.R --limit 2",
      "  Rscript code/scripts/align_all_rna_features.R --limit 2 --run",
      "  Rscript code/scripts/align_all_rna_features.R --run",
      sep = "\n"
    )
  )
}

parse_args <- function(args) {
  config <- list(
    input_list = "parameters/rna_path.txt",
    ref_genes_path = "parameters/genes_20260527.csv",
    outputdir = "data/processed/Atlas_20260721/rna",
    assay = "rna",
    limit = NA_integer_,
    run = FALSE
  )

  i <- 1L
  while (i <= length(args)) {
    arg <- args[[i]]

    if (arg %in% c("-h", "--help")) {
      usage()
      quit(save = "no", status = 0)
    } else if (arg == "--input-list") {
      i <- i + 1L
      if (i > length(args)) stop("--input-list 需要一个文件路径。")
      config$input_list <- args[[i]]
    } else if (arg == "--ref-genes") {
      i <- i + 1L
      if (i > length(args)) stop("--ref-genes 需要一个文件路径。")
      config$ref_genes_path <- args[[i]]
    } else if (arg == "--outputdir") {
      i <- i + 1L
      if (i > length(args)) stop("--outputdir 需要一个目录路径。")
      config$outputdir <- args[[i]]
    } else if (arg == "--assay") {
      i <- i + 1L
      if (i > length(args)) stop("--assay 需要一个 assay 名称。")
      config$assay <- args[[i]]
    } else if (arg == "--limit") {
      i <- i + 1L
      if (i > length(args)) stop("--limit 需要一个正整数。")
      config$limit <- as.integer(args[[i]])
      if (is.na(config$limit) || config$limit <= 0L) {
        stop("--limit 必须是正整数。")
      }
    } else if (arg == "--run") {
      config$run <- TRUE
    } else {
      stop(sprintf("未知参数: %s\n可使用 --help 查看用法。", arg))
    }

    i <- i + 1L
  }

  config
}

check_packages <- function() {
  required_pkgs <- c("Matrix", "data.table", "Seurat", "R.utils")
  missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]

  if (length(missing_pkgs) > 0L) {
    stop(sprintf("缺少依赖包: %s", paste(missing_pkgs, collapse = ", ")))
  }
}

resolve_path <- function(path, root_dir) {
  if (grepl("^(/|[A-Za-z]:[/\\\\])", path)) {
    return(normalizePath(path, winslash = "/", mustWork = FALSE))
  }

  normalizePath(file.path(root_dir, path), winslash = "/", mustWork = FALSE)
}

read_rna_paths <- function(input_list, root_dir) {
  if (!file.exists(input_list)) {
    stop(sprintf("找不到路径文件: %s", input_list))
  }

  tbl <- read.delim(
    input_list,
    sep = "\t",
    header = TRUE,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  if (!("DONOR_ID" %in% names(tbl))) {
    stop("路径文件必须包含 DONOR_ID 列。")
  }
  if (!("PATH" %in% names(tbl))) {
    stop("路径文件必须包含 PATH 列。")
  }

  keep <- nzchar(trimws(tbl[["PATH"]])) & grepl("rna\\.rds$", tbl[["PATH"]], ignore.case = TRUE)
  out <- data.frame(
    donor_id = tbl[["DONOR_ID"]][keep],
    raw_path = tbl[["PATH"]][keep],
    stringsAsFactors = FALSE
  )

  if (nrow(out) == 0L) {
    stop("没有筛选到任何 rna.rds 路径。")
  }

  out$resolved_path <- vapply(out$raw_path, resolve_path, character(1), root_dir = root_dir)

  duplicated_path <- duplicated(out$resolved_path)
  if (any(duplicated_path)) {
    warning(sprintf("发现 %d 个重复 rna.rds 路径，仅保留第一次出现的记录。", sum(duplicated_path)), call. = FALSE)
    out <- out[!duplicated_path, , drop = FALSE]
  }

  rownames(out) <- NULL
  out
}

add_output_dirs <- function(rna_paths, outputdir) {
  rna_paths$output_dir <- file.path(outputdir, rna_paths$donor_id)
  rna_paths
}

main <- function() {
  config <- parse_args(commandArgs(trailingOnly = TRUE))
  root_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  utils_path <- file.path(root_dir, "code", "src", "r", "utils.r")

  if (!file.exists(utils_path)) {
    stop(sprintf("找不到 utils.r: %s", utils_path))
  }

  config$input_list <- resolve_path(config$input_list, root_dir)
  config$ref_genes_path <- resolve_path(config$ref_genes_path, root_dir)
  config$outputdir <- resolve_path(config$outputdir, root_dir)

  if (!file.exists(config$ref_genes_path)) {
    stop(sprintf("找不到参考基因文件: %s", config$ref_genes_path))
  }

  rna_paths <- read_rna_paths(config$input_list, root_dir)
  if (!is.na(config$limit)) {
    rna_paths <- head(rna_paths, config$limit)
  }
  rna_paths <- add_output_dirs(rna_paths, config$outputdir)

  message(sprintf("Found %d rna.rds files to align.", nrow(rna_paths)))
  message(sprintf("Reference genes: %s", config$ref_genes_path))
  message(sprintf("Output directory: %s", config$outputdir))
  message(sprintf("Mode: %s", if (config$run) "run" else "dry-run"))

  if (!config$run) {
    print(rna_paths[, c("donor_id", "resolved_path", "output_dir"), drop = FALSE], row.names = FALSE)
    message("Dry-run only. Add --run to read RDS files and write output.")
    return(invisible(NULL))
  }

  check_packages()
  source(utils_path)

  failures <- list()
  success_count <- 0L

  for (i in seq_len(nrow(rna_paths))) {
    sample_row <- rna_paths[i, , drop = FALSE]
    donor_id <- as.character(sample_row$donor_id)

    message(sprintf("[%d/%d] %s", i, nrow(rna_paths), donor_id))

    result <- tryCatch({
      if (!file.exists(sample_row$resolved_path)) {
        stop(sprintf("文件不存在: %s", sample_row$resolved_path))
      }

      obj <- readRDS(sample_row$resolved_path)
      align_rna_features(
        se_obj = obj,
        ref_genes_path = config$ref_genes_path,
        output_dir = sample_row$output_dir,
        assay = config$assay
      )

      mtx_path <- file.path(sample_row$output_dir, "matrix.mtx")
      mtx_gz_path <- file.path(sample_row$output_dir, "matrix.mtx.gz")
      R.utils::gzip(
        filename = mtx_path,
        destname = mtx_gz_path,
        overwrite = TRUE,
        remove = TRUE
      )

      TRUE
    }, error = function(e) e)

    if (inherits(result, "error")) {
      failures[[length(failures) + 1L]] <- data.frame(
        donor_id = donor_id,
        path = sample_row$resolved_path,
        output_dir = sample_row$output_dir,
        error = conditionMessage(result),
        stringsAsFactors = FALSE
      )
      message(sprintf("  FAILED: %s", conditionMessage(result)))
      next
    }

    success_count <- success_count + 1L
    message(sprintf("  OK: %s", sample_row$output_dir))
  }

  message(sprintf("Completed: success=%d; failed=%d", success_count, length(failures)))

  if (length(failures) > 0L) {
    failed_df <- do.call(rbind, failures)
    failed_path <- file.path("logs", "align_all_rna_features_failed.tsv")
    dir.create(dirname(failed_path), recursive = TRUE, showWarnings = FALSE)
    utils::write.table(
      failed_df,
      file = failed_path,
      sep = "\t",
      quote = TRUE,
      row.names = FALSE,
      col.names = TRUE,
      na = ""
    )
    message(sprintf("Failed file report: %s", failed_path))
  }

  invisible(NULL)
}

main()
