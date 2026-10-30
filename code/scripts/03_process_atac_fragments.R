#!/usr/bin/env Rscript

# Process ATAC fragment files into 10x matrix directories.
# Default is dry-run. Add --run to actually process samples.

options(stringsAsFactors = FALSE)

library(Signac)
library(Seurat)
library(GenomicRanges)
library(Matrix)
library(R.utils)
library(EnsDb.Hsapiens.v86)

# ---------------------------------------------------------------------------
# Annotation (hard-coded hg38)
# ---------------------------------------------------------------------------

annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
lv <- seqlevels(annotation)
lv[lv == "MT"] <- "M"
seqlevels(annotation) <- paste0("chr", lv)
genome(annotation) <- "hg38"

# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

usage <- function() {
  cat(paste(
    "用法:",
    "  Rscript code/scripts/03_process_atac_fragments.R [OPTIONS]",
    "",
    "选项:",
    "  --input-list PATH       Fragment 路径表 (TSV, 列: DONOR_ID, PATH)",
    "                          [默认: parameters/fragments_path.txt]",
    "  --peaks-bed PATH        参考 peaks BED 文件",
    "                          [默认: parameters/reference_peaks.bed]",
    "  --output-dir PATH       输出根目录",
    "                          [默认: data/processed/Atlas_20260528/atac]",
    "  --min-counts N          nCount_peaks 下限 [默认: 1000]",
    "  --max-counts N          nCount_peaks 上限 [默认: 50000]",
    "  --min-tss N             TSS.enrichment 下限 [默认: 2]",
    "  --max-nucleosome N      nucleosome_signal 上限 [默认: 2]",
    "  --min-frip N            FRiP 下限 [默认: 0]",
    "  --save-qc               保存 QC 指标 (qc_metrics.tsv.gz)",
    "  --limit N               只处理前 N 个样本 (测试用)",
    "  --run                   实际执行 (默认 dry-run)",
    "  -h, --help              显示帮助",
    "",
    "示例:",
    "  Rscript code/scripts/03_process_atac_fragments.R --limit 2",
    "  Rscript code/scripts/03_process_atac_fragments.R --limit 2 --run",
    "  Rscript code/scripts/03_process_atac_fragments.R --run --save-qc",
    sep = "\n"
  ))
}

parse_args <- function(args) {
  config <- list(
    input_list     = "parameters/fragments_path.txt",
    peaks_bed      = "parameters/reference_peaks.bed",
    output_dir     = "data/processed/Atlas_20260528/atac",
    min_counts     = 1000,
    max_counts     = 50000,
    min_tss        = 2,
    max_nucleosome = 2,
    min_frip       = 0,
    save_qc        = FALSE,
    limit          = NA_integer_,
    run            = FALSE
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
    } else if (arg == "--peaks-bed") {
      i <- i + 1L
      if (i > length(args)) stop("--peaks-bed 需要一个文件路径。")
      config$peaks_bed <- args[[i]]
    } else if (arg == "--output-dir") {
      i <- i + 1L
      if (i > length(args)) stop("--output-dir 需要一个路径。")
      config$output_dir <- args[[i]]
    } else if (arg == "--min-counts") {
      i <- i + 1L
      if (i > length(args)) stop("--min-counts 需要一个数值。")
      config$min_counts <- as.numeric(args[[i]])
    } else if (arg == "--max-counts") {
      i <- i + 1L
      if (i > length(args)) stop("--max-counts 需要一个数值。")
      config$max_counts <- as.numeric(args[[i]])
    } else if (arg == "--min-tss") {
      i <- i + 1L
      if (i > length(args)) stop("--min-tss 需要一个数值。")
      config$min_tss <- as.numeric(args[[i]])
    } else if (arg == "--max-nucleosome") {
      i <- i + 1L
      if (i > length(args)) stop("--max-nucleosome 需要一个数值。")
      config$max_nucleosome <- as.numeric(args[[i]])
    } else if (arg == "--min-frip") {
      i <- i + 1L
      if (i > length(args)) stop("--min-frip 需要一个数值。")
      config$min_frip <- as.numeric(args[[i]])
    } else if (arg == "--save-qc") {
      config$save_qc <- TRUE
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

# ---------------------------------------------------------------------------
# Path resolution
# ---------------------------------------------------------------------------

resolve_path <- function(path) {
  if (grepl("^(/|[A-Za-z]:[/\\\\])", path)) {
    return(normalizePath(path, winslash = "/", mustWork = FALSE))
  }
  normalizePath(file.path(getwd(), path), winslash = "/", mustWork = FALSE)
}

# ---------------------------------------------------------------------------
# Fragment path table
# ---------------------------------------------------------------------------

read_frag_paths <- function(input_list) {
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

  if (!("DONOR_ID" %in% names(tbl))) stop("路径文件必须包含 DONOR_ID 列。")
  if (!("PATH" %in% names(tbl)))   stop("路径文件必须包含 PATH 列。")

  keep <- nzchar(trimws(tbl[["PATH"]]))
  out <- data.frame(
    donor_id       = tbl[["DONOR_ID"]][keep],
    fragments_path = tbl[["PATH"]][keep],
    stringsAsFactors = FALSE
  )

  if (nrow(out) == 0L) stop("没有筛选到任何有效 fragment 路径。")
  rownames(out) <- NULL
  out
}

# ---------------------------------------------------------------------------
# Core processing function (unchanged except annotation is a closure)
# ---------------------------------------------------------------------------

process_atac_frag <- function(fragments_path, peaks_bed_path, output_dir,
                              annotation,
                              min_counts     = 1000,
                              max_counts     = 50000,
                              min_tss        = 2,
                              max_nucleosome = 2,
                              min_frip       = 0,
                              save_qc        = FALSE) {

  ## 1. peaks -> GRanges
  bed <- read.table(peaks_bed_path, header = FALSE, sep = "\t",
                    col.names = c("chr", "start", "end"),
                    colClasses = c("character", "integer", "integer"))
  peaks_gr <- GRanges(seqnames = bed$chr,
                      ranges = IRanges(start = bed$start + 1L, end = bed$end))
  peak_names <- paste0(bed$chr, ":", bed$start, "-", bed$end)
  names(peaks_gr) <- peak_names

  ## 2. fragment 对象 + 全量计数
  frags <- CreateFragmentObject(path = fragments_path, validate.fragments = FALSE)
  count_mat <- FeatureMatrix(fragments = frags, features = peaks_gr)
  rownames(count_mat) <- peak_names

  ## 3. 预过滤
  cps <- Matrix::colSums(count_mat)
  count_mat <- count_mat[, cps >= min_counts, drop = FALSE]
  if (ncol(count_mat) == 0L)
    stop("预过滤后无 barcode 剩余, 检查 min_counts 是否过高")

  ## 4. 构建对象并计算 fragment-based QC
  chrom <- CreateChromatinAssay(counts = count_mat, sep = c(":", "-"),
                                fragments = fragments_path,
                                annotation = annotation)
  obj <- CreateSeuratObject(chrom, assay = "peaks")

  obj <- TSSEnrichment(obj, fast = TRUE)
  obj <- NucleosomeSignal(obj)

  fragcounts <- CountFragments(fragments = fragments_path)
  rownames(fragcounts) <- fragcounts$CB
  obj$total_frag <- fragcounts[colnames(obj), "frequency_count"]
  obj <- FRiP(obj, assay = "peaks", total.fragments = "total_frag")

  ## 5. 硬阈值过滤
  keep <- obj$nCount_peaks      >= min_counts &
          obj$TSS.enrichment    >= min_tss &
          obj$nucleosome_signal <= max_nucleosome &
          obj$FRiP              >= min_frip
  if (!is.null(max_counts)) keep <- keep & obj$nCount_peaks <= max_counts
  keep[is.na(keep)] <- FALSE

  message(sprintf("  QC: %d / %d barcodes 通过 (%.1f%%)",
                  sum(keep), length(keep), 100 * mean(keep)))

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  if (save_qc) {
    qc <- data.frame(barcode = colnames(obj),
                     nCount_peaks = obj$nCount_peaks,
                     TSS.enrichment = obj$TSS.enrichment,
                     nucleosome_signal = obj$nucleosome_signal,
                     FRiP = obj$FRiP, keep = keep)
    gzq <- gzfile(file.path(output_dir, "qc_metrics.tsv.gz"), "w")
    write.table(qc, gzq, sep = "\t", quote = FALSE, row.names = FALSE)
    close(gzq)
  }

  pass_cells <- colnames(obj)[keep]
  count_mat  <- count_mat[, pass_cells, drop = FALSE]

  ## 6. 写出 10x 格式
  gz <- gzfile(file.path(output_dir, "barcodes.tsv.gz"), "w")
  writeLines(pass_cells, gz); close(gz)

  features_df <- data.frame(id = peak_names, name = peak_names,
                            feature_type = "Peaks", stringsAsFactors = FALSE)
  gz <- gzfile(file.path(output_dir, "features.tsv.gz"), "w")
  write.table(features_df, gz, sep = "\t", quote = FALSE,
              row.names = FALSE, col.names = FALSE); close(gz)

  mtx_file <- file.path(output_dir, "matrix.mtx")
  Matrix::writeMM(count_mat, mtx_file)
  R.utils::gzip(mtx_file, destname = paste0(mtx_file, ".gz"),
                overwrite = TRUE, remove = TRUE)

  message("Done. Output written to: ", output_dir)
  invisible(count_mat)
}

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

main <- function() {
  config <- parse_args(commandArgs(trailingOnly = TRUE))

  input_list  <- resolve_path(config$input_list)
  peaks_bed   <- resolve_path(config$peaks_bed)
  output_root <- config$output_dir  # relative to cwd; resolved per-sample below

  frag_paths <- read_frag_paths(input_list)
  if (!is.na(config$limit)) {
    frag_paths <- head(frag_paths, config$limit)
  }

  message(sprintf("Found %d fragment entries.", nrow(frag_paths)))
  message(sprintf("Peaks BED: %s", peaks_bed))
  message(sprintf("Output root: %s", output_root))
  message(sprintf("QC thresholds: min_counts=%d, max_counts=%d, min_tss=%.1f, max_nucleosome=%.1f, min_frip=%.1f",
                  config$min_counts, config$max_counts,
                  config$min_tss, config$max_nucleosome, config$min_frip))
  message(sprintf("Save QC: %s", config$save_qc))
  message(sprintf("Mode: %s", if (config$run) "run" else "dry-run"))

  if (!config$run) {
    for (i in seq_len(nrow(frag_paths))) {
      row <- frag_paths[i, , drop = FALSE]
      out_dir <- file.path(output_root, row$donor_id)
      message(sprintf("  [%d/%d] %s -> %s", i, nrow(frag_paths), row$donor_id, out_dir))
    }
    message("Dry-run only. Add --run to process fragment files.")
    return(invisible(NULL))
  }

  failures <- list()
  status   <- character(nrow(frag_paths))

  for (i in seq_len(nrow(frag_paths))) {
    row <- frag_paths[i, , drop = FALSE]
    donor_id       <- row$donor_id
    fragments_path <- row$fragments_path
    out_dir        <- file.path(output_root, donor_id)

    message(sprintf("[%d/%d] Processing %s", i, nrow(frag_paths), donor_id))

    result <- tryCatch(
      {
        process_atac_frag(
          fragments_path = fragments_path,
          peaks_bed_path = peaks_bed,
          output_dir     = out_dir,
          annotation     = annotation,
          min_counts     = config$min_counts,
          max_counts     = config$max_counts,
          min_tss        = config$min_tss,
          max_nucleosome = config$max_nucleosome,
          min_frip       = config$min_frip,
          save_qc        = config$save_qc
        )
        "OK"
      },
      error = function(e) e
    )

    if (inherits(result, "error")) {
      status[[i]] <- "FAILED"
      failures[[length(failures) + 1L]] <- data.frame(
        donor_id       = donor_id,
        fragments_path = fragments_path,
        error          = conditionMessage(result),
        stringsAsFactors = FALSE
      )
      message(sprintf("  FAILED: %s", conditionMessage(result)))
    } else {
      status[[i]] <- "OK"
    }
  }

  message(sprintf(
    "Completed: OK=%d; failed=%d",
    sum(status == "OK"), length(failures)
  ))

  if (length(failures) > 0L) {
    failed_df <- do.call(rbind, failures)
    print(failed_df, row.names = FALSE)
    quit(save = "no", status = 1)
  }

  invisible(NULL)
}

main()
