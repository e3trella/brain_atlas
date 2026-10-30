#!/usr/bin/env Rscript

# Count file-level gene presence across RNA Seurat RDS files listed in rna_path.txt.

options(stringsAsFactors = FALSE)

usage <- function() {
  cat(
    paste(
      "用法:",
      "  Rscript code/scripts/count_rna_gene_frequency.R [--input-list PATH]",
      "                                                 [--output-dir DIR] [--limit N] [--quiet]",
      "",
      "说明:",
      "  1. 只处理路径列表中以 rna.rds 结尾的文件，忽略 atac.h5seurat。",
      "  2. 直接读取 Seurat 对象的 rownames 作为基因名。",
      "  3. 输出文件级出现频率: frequency = n_files / 成功读取的 rna.rds 文件数。",
      "  4. 单个文件失败不会中断全流程，失败记录写入 failed_files.tsv。",
      "",
      "示例:",
      "  Rscript code/scripts/count_rna_gene_frequency.R --limit 2",
      "  Rscript code/scripts/count_rna_gene_frequency.R",
      "  Rscript code/scripts/count_rna_gene_frequency.R --output-dir logs/rna_gene_frequency",
      sep = "\n"
    )
  )
}

parse_args <- function(args) {
  config <- list(
    input_list = "parameters/rna_path.txt",
    output_dir = file.path("logs", "rna_gene_frequency"),
    limit = NA_integer_,
    quiet = FALSE
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
    } else if (arg == "--output-dir") {
      i <- i + 1L
      if (i > length(args)) stop("--output-dir 需要一个目录路径。")
      config$output_dir <- args[[i]]
    } else if (arg == "--limit") {
      i <- i + 1L
      if (i > length(args)) stop("--limit 需要一个正整数。")
      config$limit <- as.integer(args[[i]])
      if (is.na(config$limit) || config$limit <= 0L) {
        stop("--limit 必须是正整数。")
      }
    } else if (arg == "--quiet") {
      config$quiet <- TRUE
    } else {
      stop(sprintf("未知参数: %s\n可使用 --help 查看用法。", arg))
    }

    i <- i + 1L
  }

  config
}

resolve_path <- function(path, root_dir) {
  if (grepl("^(/|[A-Za-z]:[/\\\\])", path)) {
    return(normalizePath(path, winslash = "/", mustWork = FALSE))
  }

  normalizePath(file.path(root_dir, path), winslash = "/", mustWork = FALSE)
}

write_tsv <- function(df, path) {
  utils::write.table(
    df,
    file = path,
    sep = "\t",
    quote = TRUE,
    row.names = FALSE,
    col.names = TRUE,
    na = ""
  )
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

  if ("PATH" %in% names(tbl)) {
    path_col <- tbl[["PATH"]]
  } else if (ncol(tbl) >= 2L) {
    path_col <- tbl[[2L]]
  } else {
    stop("路径文件至少需要两列，且第二列应为 PATH。")
  }

  donor_col <- if (ncol(tbl) >= 1L) tbl[[1L]] else seq_len(nrow(tbl))
  keep <- nzchar(trimws(path_col)) & grepl("rna\\.rds$", path_col, ignore.case = TRUE)

  out <- data.frame(
    donor_id = donor_col[keep],
    raw_path = path_col[keep],
    stringsAsFactors = FALSE
  )

  if (nrow(out) == 0L) {
    stop("没有筛选到任何 rna.rds 路径。")
  }

  out$resolved_path <- vapply(out$raw_path, resolve_path, character(1), root_dir = root_dir)

  duplicated_path <- duplicated(out$resolved_path)
  if (any(duplicated_path)) {
    warning(
      sprintf("发现 %d 个重复 rna.rds 路径，仅保留第一次出现的记录。", sum(duplicated_path)),
      call. = FALSE
    )
    out <- out[!duplicated_path, , drop = FALSE]
  }

  rownames(out) <- NULL
  out
}

read_sample_genes <- function(sample_row) {
  if (!file.exists(sample_row$resolved_path)) {
    stop(sprintf("文件不存在: %s", sample_row$resolved_path))
  }

  obj <- readRDS(sample_row$resolved_path)
  genes <- rownames(obj)

  if (is.null(genes) || length(genes) == 0L) {
    stop("RDS 对象没有 rownames。")
  }

  unique(genes)
}

make_frequency_table <- function(gene_to_files, denominator) {
  genes <- names(gene_to_files)

  out <- data.frame(
    gene = genes,
    n_files = lengths(gene_to_files),
    frequency = lengths(gene_to_files) / denominator,
    stringsAsFactors = FALSE
  )

  out[order(-out$n_files, out$gene), , drop = FALSE]
}

main <- function() {
  config <- parse_args(commandArgs(trailingOnly = TRUE))
  root_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

  rna_paths <- read_rna_paths(config$input_list, root_dir)

  if (!is.na(config$limit)) {
    rna_paths <- head(rna_paths, config$limit)
  }

  dir.create(config$output_dir, recursive = TRUE, showWarnings = FALSE)

  gene_to_files <- list()
  failures <- list()
  success_count <- 0L

  if (!config$quiet) {
    message(sprintf("Found %d rna.rds files to scan.", nrow(rna_paths)))
  }

  for (i in seq_len(nrow(rna_paths))) {
    sample_row <- rna_paths[i, , drop = FALSE]
    donor_id <- as.character(sample_row$donor_id)

    if (!config$quiet) {
      message(sprintf("[%d/%d] %s", i, nrow(rna_paths), donor_id))
    }

    result <- tryCatch(
      read_sample_genes(sample_row),
      error = function(e) e
    )

    if (inherits(result, "error")) {
      failures[[length(failures) + 1L]] <- data.frame(
        donor_id = donor_id,
        path = sample_row$resolved_path,
        error = conditionMessage(result),
        stringsAsFactors = FALSE
      )
      next
    }

    success_count <- success_count + 1L
    for (gene in result) {
      gene_to_files[[gene]] <- c(gene_to_files[[gene]], donor_id)
    }
  }

  if (success_count == 0L) {
    stop("没有成功读取任何 rna.rds，无法计算基因出现频率。")
  }

  freq_tbl <- make_frequency_table(gene_to_files, success_count)
  freq_path <- file.path(config$output_dir, "gene_frequency.tsv")
  write_tsv(freq_tbl, freq_path)

  if (length(failures) > 0L) {
    failed_tbl <- do.call(rbind, failures)
  } else {
    failed_tbl <- data.frame(
      donor_id = character(),
      path = character(),
      error = character(),
      stringsAsFactors = FALSE
    )
  }

  failed_path <- file.path(config$output_dir, "failed_files.tsv")
  write_tsv(failed_tbl, failed_path)

  summary_tbl <- data.frame(
    input_list = config$input_list,
    total_rna_files = nrow(rna_paths),
    successful_files = success_count,
    failed_files = length(failures),
    unique_genes = nrow(freq_tbl),
    gene_frequency_path = normalizePath(freq_path, winslash = "/", mustWork = FALSE),
    failed_files_path = normalizePath(failed_path, winslash = "/", mustWork = FALSE),
    stringsAsFactors = FALSE
  )
  write_tsv(summary_tbl, file.path(config$output_dir, "summary.tsv"))

  message(sprintf("Done. successful=%d; failed=%d; unique_genes=%d", success_count, length(failures), nrow(freq_tbl)))
  message(sprintf("Wrote: %s", freq_path))
  message(sprintf("Wrote: %s", failed_path))
}

main()
