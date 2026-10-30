#!/usr/bin/env Rscript

# 批量更新 parameters/rna_path.txt 中 RNA RDS 对象的基因名。
# 默认 dry-run，只做读取、映射和统计，不覆盖保存。
# 加 --run 后，按原路径覆盖保存更新后的 rna.rds。

options(stringsAsFactors = FALSE)

usage <- function() {
  cat(
    paste(
      "用法:",
      "  Rscript code/scripts/update_rna_hgnc_symbols.R [--input-list PATH] [--limit N]",
      "                                                 [--assay NAME] [--log-dir DIR] [--run]",
      "",
      "说明:",
      "  1. 只处理路径列表中以 rna.rds 结尾的文件，忽略 atac.h5seurat。",
      "  2. 使用 scCustomize::Updated_HGNC_Symbols() 统一基因名。",
      "  3. 如果一个旧基因名对应多个新基因名，只取第一条映射。",
      "  4. 更新后仅保留 counts；data 和 scale.data 会被丢弃。",
      "  5. 默认会输出详细日志到 logs/update_rna_hgnc_symbols/<timestamp>/。",
      "  6. 默认是 dry-run，不覆盖 rna.rds；加 --run 才会原路径覆盖保存。",
      "",
      "示例:",
      "  Rscript code/scripts/update_rna_hgnc_symbols.R",
      "  Rscript code/scripts/update_rna_hgnc_symbols.R --limit 2",
      "  Rscript code/scripts/update_rna_hgnc_symbols.R --run",
      "  Rscript code/scripts/update_rna_hgnc_symbols.R --log-dir logs/my_hgnc_logs",
      sep = "\n"
    )
  )
}

parse_args <- function(args) {
  config <- list(
    input_list = "parameters/rna_path.txt",
    assay = "rna",
    limit = NA_integer_,
    log_dir = file.path("logs", "update_rna_hgnc_symbols"),
    run = FALSE,
    case_check_as_warn = TRUE,
    update_symbol_data = NULL,
    verbose = TRUE
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
    } else if (arg == "--assay") {
      i <- i + 1L
      if (i > length(args)) stop("--assay 需要一个 assay 名称。")
      config$assay <- args[[i]]
    } else if (arg == "--limit") {
      i <- i + 1L
      if (i > length(args)) stop("--limit 需要一个整数。")
      config$limit <- as.integer(args[[i]])
      if (is.na(config$limit) || config$limit <= 0L) {
        stop("--limit 必须是正整数。")
      }
    } else if (arg == "--log-dir") {
      i <- i + 1L
      if (i > length(args)) stop("--log-dir 需要一个目录路径。")
      config$log_dir <- args[[i]]
    } else if (arg == "--run") {
      config$run <- TRUE
    } else if (arg == "--refresh-hgnc") {
      config$update_symbol_data <- TRUE
    } else if (arg == "--quiet") {
      config$verbose <- FALSE
    } else {
      stop(sprintf("未知参数: %s\n可使用 --help 查看用法。", arg))
    }

    i <- i + 1L
  }

  config
}

check_packages <- function() {
  required_pkgs <- c("BiocFileCache", "Matrix", "Seurat", "SeuratObject", "scCustomize")
  missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]

  if (length(missing_pkgs) > 0L) {
    stop(
      sprintf(
        "缺少依赖包: %s",
        paste(missing_pkgs, collapse = ", ")
      )
    )
  }
}

resolve_path <- function(path, root_dir) {
  if (grepl("^(/|[A-Za-z]:[/\\\\])", path)) {
    return(normalizePath(path, winslash = "/", mustWork = FALSE))
  }

  normalizePath(file.path(root_dir, path), winslash = "/", mustWork = FALSE)
}

safe_path_component <- function(x) {
  x <- gsub("[^A-Za-z0-9._-]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  if (!nzchar(x)) {
    x <- "unknown"
  }
  x
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
  tbl_out <- data.frame(
    donor_id = donor_col[keep],
    raw_path = path_col[keep],
    stringsAsFactors = FALSE
  )

  if (nrow(tbl_out) == 0L) {
    stop("没有筛选到任何 rna.rds 路径。")
  }

  tbl_out$resolved_path <- vapply(tbl_out$raw_path, resolve_path, character(1), root_dir = root_dir)

  duplicated_path <- duplicated(tbl_out$resolved_path)
  if (any(duplicated_path)) {
    warning(
      sprintf(
        "发现 %d 个重复 rna.rds 路径，仅保留第一次出现的记录。",
        sum(duplicated_path)
      ),
      call. = FALSE
    )
    tbl_out <- tbl_out[!duplicated_path, , drop = FALSE]
  }

  rownames(tbl_out) <- NULL
  tbl_out
}

collapse_sparse_rows <- function(mat, new_names) {
  if (nrow(mat) != length(new_names)) {
    stop("矩阵行数与 new_names 长度不一致。")
  }

  if (!anyDuplicated(new_names)) {
    rownames(mat) <- new_names
    return(as(mat, "dgCMatrix"))
  }

  unique_names <- unique(new_names)
  group_idx <- match(new_names, unique_names)
  triplet <- Matrix::summary(mat)

  collapsed <- Matrix::sparseMatrix(
    i = group_idx[triplet$i],
    j = triplet$j,
    x = triplet$x,
    dims = c(length(unique_names), ncol(mat)),
    dimnames = list(unique_names, colnames(mat))
  )

  as(collapsed, "dgCMatrix")
}

prepare_symbol_mapping <- function(feature_names, update_symbol_data = NULL, case_check_as_warn = TRUE, verbose = TRUE) {
  mapping_raw <- scCustomize::Updated_HGNC_Symbols(
    input_data = feature_names,
    update_symbol_data = update_symbol_data,
    case_check_as_warn = case_check_as_warn,
    verbose = verbose
  )

  mapping_raw$candidate_rank <- ave(
    rep.int(1L, nrow(mapping_raw)),
    mapping_raw$input_features,
    FUN = seq_along
  )

  ambiguous_features <- unique(mapping_raw$input_features[duplicated(mapping_raw$input_features)])
  mapping_raw$is_ambiguous_input <- mapping_raw$input_features %in% ambiguous_features
  mapping_raw$is_selected <- mapping_raw$candidate_rank == 1L

  mapping_first <- mapping_raw[!duplicated(mapping_raw$input_features), , drop = FALSE]
  mapping_first <- mapping_first[match(feature_names, mapping_first$input_features), , drop = FALSE]
  mapping_first$selected_candidate_rank <- mapping_first$candidate_rank

  output_features <- mapping_first$Output_Features
  invalid_output <- is.na(output_features) | !nzchar(output_features)
  output_features[invalid_output] <- feature_names[invalid_output]

  list(
    mapping_raw = mapping_raw,
    mapping = mapping_first,
    output_features = output_features,
    ambiguous_features = ambiguous_features
  )
}

update_seurat_rna_assay <- function(seurat_obj, assay_name, update_symbol_data = NULL, case_check_as_warn = TRUE, verbose = TRUE) {
  if (!inherits(seurat_obj, "Seurat")) {
    stop("readRDS 读入对象不是 Seurat。")
  }

  assay_names <- names(seurat_obj@assays)
  if (!assay_name %in% assay_names) {
    stop(sprintf("对象中不存在 assay '%s'。可用 assays: %s", assay_name, paste(assay_names, collapse = ", ")))
  }

  counts <- SeuratObject::GetAssayData(seurat_obj, assay = assay_name, layer = "counts")
  counts <- as(counts, "dgCMatrix")
  original_features <- rownames(counts)

  duplicated_input_features <- unique(original_features[duplicated(original_features)])
  if (length(duplicated_input_features) > 0L) {
    counts <- collapse_sparse_rows(counts, original_features)
    original_features <- rownames(counts)
  }

  mapping_info <- prepare_symbol_mapping(
    feature_names = original_features,
    update_symbol_data = update_symbol_data,
    case_check_as_warn = case_check_as_warn,
    verbose = verbose
  )

  new_counts <- collapse_sparse_rows(counts, mapping_info$output_features)
  merge_sizes <- table(mapping_info$output_features)
  merged_output_features <- names(merge_sizes[merge_sizes > 1L])

  mapping_selected <- mapping_info$mapping
  mapping_selected$mapping_status <- ifelse(
    !is.na(mapping_selected$Updated_Symbol),
    "updated",
    ifelse(
      !is.na(mapping_selected$Approved_Symbol),
      "approved",
      ifelse(
        !is.na(mapping_selected$Not_Found_Symbol),
        "not_found",
        "unknown"
      )
    )
  )
  mapping_selected$final_feature <- mapping_info$output_features
  mapping_selected$is_ambiguous_input <- mapping_selected$input_features %in% mapping_info$ambiguous_features
  mapping_selected$is_merged_output <- mapping_selected$final_feature %in% merged_output_features
  mapping_selected$merge_group_size <- unname(as.integer(merge_sizes[mapping_selected$final_feature]))

  ambiguous_log <- mapping_info$mapping_raw[mapping_info$mapping_raw$is_ambiguous_input, , drop = FALSE]

  if (length(merged_output_features) > 0L) {
    merged_log <- do.call(
      rbind,
      lapply(merged_output_features, function(feature_name) {
        subset_df <- mapping_selected[mapping_selected$final_feature == feature_name, , drop = FALSE]
        data.frame(
          final_feature = feature_name,
          n_input_features = nrow(subset_df),
          input_features = paste(subset_df$input_features, collapse = ";"),
          mapping_status = paste(subset_df$mapping_status, collapse = ";"),
          updated_symbol = paste(unique(stats::na.omit(subset_df$Updated_Symbol)), collapse = ";"),
          approved_symbol = paste(unique(stats::na.omit(subset_df$Approved_Symbol)), collapse = ";"),
          not_found_symbol = paste(unique(stats::na.omit(subset_df$Not_Found_Symbol)), collapse = ";"),
          stringsAsFactors = FALSE
        )
      })
    )
  } else {
    merged_log <- data.frame(
      final_feature = character(0),
      n_input_features = integer(0),
      input_features = character(0),
      mapping_status = character(0),
      updated_symbol = character(0),
      approved_symbol = character(0),
      not_found_symbol = character(0),
      stringsAsFactors = FALSE
    )
  }

  old_assay <- seurat_obj[[assay_name]]
  old_key <- tryCatch(SeuratObject::Key(old_assay), error = function(e) NULL)
  old_var_features <- tryCatch(SeuratObject::VariableFeatures(old_assay), error = function(e) character(0))

  new_assay <- SeuratObject::CreateAssayObject(counts = new_counts, min.cells = 0, min.features = 0)

  if (!is.null(old_key) && nzchar(old_key)) {
    SeuratObject::Key(new_assay) <- old_key
  }

  if (length(old_var_features) > 0L) {
    mapped_var_features <- mapping_info$output_features[match(old_var_features, original_features)]
    mapped_var_features <- unique(mapped_var_features[!is.na(mapped_var_features)])
    mapped_var_features <- mapped_var_features[mapped_var_features %in% rownames(new_counts)]
    SeuratObject::VariableFeatures(new_assay) <- mapped_var_features
  }

  suppressWarnings({
    seurat_obj[[assay_name]] <- new_assay
  })
  SeuratObject::DefaultAssay(seurat_obj) <- assay_name

  stats <- list(
    input_features = length(original_features),
    output_features = nrow(new_counts),
    updated_features = sum(!is.na(mapping_info$mapping$Updated_Symbol)),
    approved_features = sum(!is.na(mapping_info$mapping$Approved_Symbol)),
    not_found_features = sum(!is.na(mapping_info$mapping$Not_Found_Symbol)),
    ambiguous_input_features = length(mapping_info$ambiguous_features),
    duplicated_input_features = length(duplicated_input_features),
    merged_features = length(original_features) - nrow(new_counts)
  )

  list(
    object = seurat_obj,
    stats = stats,
    logs = list(
      selected_mapping = mapping_selected,
      updated_genes = mapping_selected[!is.na(mapping_selected$Updated_Symbol), , drop = FALSE],
      ambiguous_candidates = ambiguous_log,
      merged_features = merged_log
    )
  )
}

format_stats <- function(stats) {
  paste(
    sprintf("features %d -> %d", stats$input_features, stats$output_features),
    sprintf("updated=%d", stats$updated_features),
    sprintf("approved=%d", stats$approved_features),
    sprintf("not_found=%d", stats$not_found_features),
    sprintf("ambiguous_first_pick=%d", stats$ambiguous_input_features),
    sprintf("merged=%d", stats$merged_features),
    sep = "; "
  )
}

create_run_log_dir <- function(base_dir, run_mode) {
  run_label <- sprintf(
    "%s_%s",
    tolower(run_mode),
    format(Sys.time(), "%Y%m%d_%H%M%S")
  )
  run_log_dir <- file.path(base_dir, run_label)
  dir.create(run_log_dir, recursive = TRUE, showWarnings = FALSE)
  run_log_dir
}

write_sample_logs <- function(log_root, sample_idx, donor_id, run_result, raw_path, resolved_path) {
  sample_dir <- file.path(
    log_root,
    sprintf("%03d_%s", sample_idx, safe_path_component(donor_id))
  )
  dir.create(sample_dir, recursive = TRUE, showWarnings = FALSE)

  selected_mapping_path <- file.path(sample_dir, "selected_mapping.tsv")
  updated_genes_path <- file.path(sample_dir, "updated_genes.tsv")
  ambiguous_candidates_path <- file.path(sample_dir, "ambiguous_candidates.tsv")
  merged_features_path <- file.path(sample_dir, "merged_features.tsv")
  sample_summary_path <- file.path(sample_dir, "sample_summary.tsv")

  write_tsv(run_result$logs$selected_mapping, selected_mapping_path)
  write_tsv(run_result$logs$updated_genes, updated_genes_path)
  write_tsv(run_result$logs$ambiguous_candidates, ambiguous_candidates_path)
  write_tsv(run_result$logs$merged_features, merged_features_path)

  sample_summary <- data.frame(
    donor_id = donor_id,
    raw_path = raw_path,
    resolved_path = resolved_path,
    input_features = run_result$stats$input_features,
    output_features = run_result$stats$output_features,
    updated_features = run_result$stats$updated_features,
    approved_features = run_result$stats$approved_features,
    not_found_features = run_result$stats$not_found_features,
    ambiguous_first_pick = run_result$stats$ambiguous_input_features,
    merged_features = run_result$stats$merged_features,
    stringsAsFactors = FALSE
  )
  write_tsv(sample_summary, sample_summary_path)

  list(
    sample_dir = sample_dir,
    selected_mapping = selected_mapping_path,
    updated_genes = updated_genes_path,
    ambiguous_candidates = ambiguous_candidates_path,
    merged_features = merged_features_path,
    sample_summary = sample_summary_path
  )
}

main <- function() {
  config <- parse_args(commandArgs(trailingOnly = TRUE))
  check_packages()

  root_dir <- getwd()
  path_table <- read_rna_paths(config$input_list, root_dir = root_dir)

  if (!is.na(config$limit)) {
    path_table <- head(path_table, config$limit)
  }

  mode_label <- if (config$run) "RUN" else "DRY-RUN"
  run_log_dir <- create_run_log_dir(config$log_dir, mode_label)
  cat(sprintf("模式: %s\n", mode_label))
  cat(sprintf("输入列表: %s\n", config$input_list))
  cat(sprintf("目标 assay: %s\n", config$assay))
  cat(sprintf("日志目录: %s\n", run_log_dir))
  cat(sprintf("待处理 rna.rds 文件数: %d\n\n", nrow(path_table)))

  results <- vector("list", length = nrow(path_table))

  for (idx in seq_len(nrow(path_table))) {
    donor_id <- path_table$donor_id[[idx]]
    raw_path <- path_table$raw_path[[idx]]
    resolved_path <- path_table$resolved_path[[idx]]

    cat(sprintf("[%d/%d] %s\n", idx, nrow(path_table), donor_id))
    cat(sprintf("  path: %s\n", raw_path))

    if (!file.exists(resolved_path)) {
      msg <- "文件不存在"
      cat(sprintf("  status: FAIL - %s\n\n", msg))
      results[[idx]] <- data.frame(
        donor_id = donor_id,
        path = raw_path,
        status = "fail",
        message = msg,
        log_dir = NA_character_,
        stringsAsFactors = FALSE
      )
      next
    }

    run_result <- tryCatch({
      obj <- readRDS(resolved_path)
      updated <- update_seurat_rna_assay(
        seurat_obj = obj,
        assay_name = config$assay,
        update_symbol_data = config$update_symbol_data,
        case_check_as_warn = config$case_check_as_warn,
        verbose = config$verbose
      )

      if (config$run) {
        saveRDS(updated$object, file = resolved_path)
      }

      updated
    }, error = function(e) {
      e
    })

    if (inherits(run_result, "error")) {
      cat(sprintf("  status: FAIL - %s\n\n", conditionMessage(run_result)))
      results[[idx]] <- data.frame(
        donor_id = donor_id,
        path = raw_path,
        status = "fail",
        message = conditionMessage(run_result),
        log_dir = NA_character_,
        stringsAsFactors = FALSE
      )
      next
    }

    log_paths <- write_sample_logs(
      log_root = run_log_dir,
      sample_idx = idx,
      donor_id = donor_id,
      run_result = run_result,
      raw_path = raw_path,
      resolved_path = resolved_path
    )

    cat(sprintf("  status: OK - %s\n", format_stats(run_result$stats)))
    cat(sprintf("  logs: %s\n", log_paths$sample_dir))
    if (!config$run) {
      cat("  action: dry-run，未写回文件\n\n")
    } else {
      cat("  action: 已按原路径覆盖保存\n\n")
    }

    results[[idx]] <- data.frame(
      donor_id = donor_id,
      path = raw_path,
      status = "ok",
      message = format_stats(run_result$stats),
      log_dir = log_paths$sample_dir,
      stringsAsFactors = FALSE
    )

    invisible(gc())
  }

  summary_df <- do.call(rbind, results)
  n_ok <- sum(summary_df$status == "ok")
  n_fail <- sum(summary_df$status == "fail")
  run_summary_path <- file.path(run_log_dir, "run_summary.tsv")
  write_tsv(summary_df, run_summary_path)

  cat("处理完成。\n")
  cat(sprintf("成功: %d\n", n_ok))
  cat(sprintf("失败: %d\n", n_fail))
  cat(sprintf("运行汇总: %s\n", run_summary_path))

  if (n_fail > 0L) {
    cat("\n失败列表:\n")
    fail_df <- summary_df[summary_df$status == "fail", , drop = FALSE]
    for (i in seq_len(nrow(fail_df))) {
      cat(sprintf("  - %s | %s | %s\n", fail_df$donor_id[[i]], fail_df$path[[i]], fail_df$message[[i]]))
    }
    quit(save = "no", status = 1)
  }
}

main()
