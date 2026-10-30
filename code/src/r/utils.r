#' rename file.path function
pj <- file.path

# rename print function
prt <- function(...) {
  cat(paste0(..., "\n"))
}

# 定义一个mkdir函数，用于创建目录，如果需要，可以先删除旧目录
#' @param directory: 要创建的目录路径
#' @param remove_old: 是否删除旧目录
#' @return : 无返回值
mkdir <- function(directory, remove_old = FALSE) {
  # 如果remove_old为真，则尝试删除旧的目录
  if (remove_old) {
    # 检查目录是否存在
    if (dir.exists(directory)) {
      # 打印正在删除的目录信息
      print("Removing directory ", directory)
      # 删除目录及其内容
      unlink(directory, recursive = TRUE)
    }
  }
  # 检查目录是否存在，如果不存在，则创建新目录
  if (!dir.exists(directory)) {
    # 创建目录，如果目录结构有多层，会递归创建
    dir.create(directory, recursive = TRUE)
  }
}

# 定义一个名为mkdirs的函数，用于批量创建目录
mkdirs <- function(directories, remove_old = FALSE) {
  # 遍历目录列表，为每个目录创建文件夹
  for (directory in directories) {
    # 调用mkdir函数创建目录，如果remove_old为TRUE，则删除旧目录后再创建
    mkdir(directory, remove_old = remove_old)
  }
}

random_round <- function(mat) {
  mat_floor <- floor(mat)
  res <- mat - mat_floor
  res[] <- rbinom(n = nrow(res) * ncol(res), size = 1, prob = res)
  mat <- mat_floor + res
  mode(mat) <- "integer"
  return(mat)
}

# 获取路径中的文件名
get_file_name <- function(path) {
  basename(path)
}

# 获取路径中最后一个文件夹的名称
get_root_folder_name <- function(path) {
  # 如果路径指向的是文件，使用dirname去除文件部分
  dir_path <- if (file.info(path)$isdir) path else dirname(path)

  # 使用basename获取最后一个文件夹名
  basename(normalizePath(dir_path, mustWork = FALSE))
}

# 获取路径中的文件后缀
get_file_suffix <- function(path) {
  # 使用strsplit函数将路径按"."分割，得到一个包含多个字符串的向量
  splits <- strsplit(path, "\\.")[[1]]
  # 如果向量长度大于1，则最后一个字符串即为文件后缀
  if (length(splits) > 1) {
    return(splits[length(splits)])
  } else {
    return("")
  }
}

# 计算多个列表之间的交集
# 函数接受一个列表作为输入，返回所有列表的交集
list_intersection <- function(lst) {
  if (length(lst) == 0) {
    return(NULL)  # 如果列表为空，返回 NULL
  }

  # 使用 Reduce 和 intersect 函数计算交集
  result <- Reduce(intersect, lst)
  return(result)
}

# 函数用于根据给定的样本路径列表拆分出不同类型的样本路径
# @param sample_path_list: 样本路径列表
# @return : 返回一个包含 RNA、ATAC 和 ADT 路径的列表
split_sample_paths <- function(sample_path_list) {

  # 筛选 RNA 列表（包含 10x_multiome 和 CITE_seq 的路径）
  rna_list <- sample_path_list[grepl("10x_multiome|CITE_seq", sample_path_list)]
  # 筛选 ATAC 列表（仅包含 10x_multiome 的路径）
  atac_list <- sample_path_list[grepl("10x_multiome", sample_path_list)]
  # 筛选 ADT 列表（仅包含 CITE_seq 的路径）
  adt_list <- sample_path_list[grepl("CITE_seq", sample_path_list)]

  return(list(rna = rna_list, atac = atac_list, adt = adt_list))
}

# check a vector has duplicated elements,
# if yes, print the duplicated elements and return TRUE, otherwise return FALSE
check_duplicated_elements <- function(vec) {
  if (length(vec) != length(unique(vec))) {
    duplicated_elements <- vec[duplicated(vec)]
    print(paste0("Duplicated elements: ", duplicated_elements))
    return(TRUE)
  } else {
    return(FALSE)
  }
}

#' Align RNA features in a Seurat object to a reference gene set
#'
#' This function aligns the RNA count matrix in one Seurat object to a common
#' reference gene set and writes standard 10x Genomics v3 files.
#'
#' @param se_obj Seurat object containing an RNA assay.
#' @param ref_genes_path Path to the reference gene list file. Supports .rds or
#'        .csv formats. For .csv files, genes are read from the second column.
#' @param output_dir Directory where output files will be written. Created if it
#'        does not exist.
#' @param assay Assay name used to extract raw counts. Default is \code{"rna"}.
#'
#' @return Invisibly returns NULL. Writes four files to \code{output_dir}:
#'         \itemize{
#'           \item \code{matrix.mtx}: Sparse count matrix in genes x cells order
#'           \item \code{barcodes.tsv.gz}: Cell barcodes
#'           \item \code{features.tsv.gz}: 10x v3 feature table with gene_id,
#'                 gene_name, and feature_type columns
#'           \item \code{rna_ref_presence.csv.gz}: Binary indicator of reference
#'                 gene presence in the input object
#'         }
#'
#' @details
#' Genes not present in the input object have all-zero counts in the output
#' matrix. The matrix is written in standard 10x orientation: genes x cells.
#'
#' @importFrom Matrix Matrix writeMM
#' @importFrom data.table fwrite
#' @importFrom Seurat GetAssayData
#'
#' @examples
#' \dontrun{
#' align_rna_features(
#'   se_obj = rna_obj,
#'   ref_genes_path = "reference_genes.rds",
#'   output_dir = "sample1/filtered_feature_bc_matrix"
#' )
#' }
#'
#' @export
align_rna_features <- function(se_obj,
                               ref_genes_path,
                               output_dir,
                               assay = "rna") {
  if (missing(ref_genes_path) || is.null(ref_genes_path)) {
    stop("ref_genes_path must be provided")
  }
  if (!file.exists(ref_genes_path)) stop("ref_genes_path does not exist")
  if (missing(output_dir) || is.null(output_dir)) {
    stop("output_dir must be provided")
  }

  library(Matrix)
  library(data.table)
  library(Seurat)

  # Auto-detect file format and read accordingly
  if (grepl("\\.rds$", ref_genes_path, ignore.case = TRUE)) {
    ref_genes <- readRDS(ref_genes_path)
  } else if (grepl("\\.csv$", ref_genes_path, ignore.case = TRUE)) {
    ref_genes <- read.csv(ref_genes_path)
    if (!("Output_Features" %in% colnames(ref_genes))) {
      stop("CSV reference gene file must contain an Output_Features column")
    }
    ref_genes <- ref_genes[["Output_Features"]]
  } else {
    stop("Unsupported file format. Expected .rds or .csv")
  }
  ref_genes <- as.character(ref_genes)
  if (length(ref_genes) == 0) stop("ref_genes is empty")

  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  counts <- Seurat::GetAssayData(se_obj, layer = "counts", assay = assay)
  counts <- as(counts, "dgCMatrix")

  n_cells <- ncol(counts)
  out_counts <- Matrix::Matrix(0, nrow = length(ref_genes), ncol = n_cells, sparse = TRUE)
  rownames(out_counts) <- ref_genes
  colnames(out_counts) <- colnames(counts)

  present_genes <- intersect(rownames(counts), ref_genes)
  if (length(present_genes) > 0) {
    rows_in_ref <- match(present_genes, ref_genes)
    out_counts[rows_in_ref, , drop = FALSE] <- counts[present_genes, , drop = FALSE]
  }

  presence_vec <- as.integer(ref_genes %in% rownames(counts))

  Matrix::writeMM(out_counts, file = pj(output_dir, "matrix.mtx"))

  data.table::fwrite(data.frame(barcode = colnames(out_counts)),
                     file = pj(output_dir, "barcodes.tsv.gz"),
                     sep = "\t",
                     col.names = FALSE)

  features <- data.frame(
    gene_id = ref_genes,
    gene_name = ref_genes,
    feature_type = "Gene Expression"
  )
  data.table::fwrite(features,
                     file = pj(output_dir, "features.tsv.gz"),
                     sep = "\t",
                     col.names = FALSE)

  presence <- data.frame(gene = ref_genes, presence = presence_vec)
  data.table::fwrite(presence,
                     file = pj(output_dir, "rna_ref_presence.csv.gz"))

  invisible(NULL)
}

#' Move and reorganize RNA matrix files from MIDAS output to processed directory
#'
#' This function moves RNA count matrices (rna.mtx) and gene mask 
#' (rna_ref_presence.csv) from raw MIDAS output folders to a standardized processed
#' directory structure. Files are organized into 'mat' and 'mask' subdirectories.
#'
#' @param source_dir Character string specifying the source directory containing
#'        subset folders with MIDAS outputs. Default is "data/raw/20260221_midas".
#' @param dest_base Character string specifying the base destination directory
#'        for processed data. Default is "data/processed/Atlas_20260221/data".
#'
#' @return Invisibly returns a list with:
#'         \itemize{
#'           \item \code{processed}: Number of folders successfully processed
#'           \item \code{skipped}: Number of folders skipped (missing rna.mtx)
#'         }
#'
#' @details
#' The function performs the following operations for each subset folder:
#' \enumerate{
#'   \item Scans \code{source_dir} for folders matching pattern "^subset"
#'   \item Checks for existence of \code{rna.mtx} file (required)
#'   \item Creates destination directory structure: \code{dest_base/<subset>/mat/} 
#'         and \code{dest_base/<subset>/mask/}
#'   \item Moves \code{rna.mtx} to the \code{mat/} subdirectory
#'   \item Moves and renames \code{rna_ref_presence.csv} to \code{mask/rna.csv}
#' }
#'
#' Folder naming convention:
#' \itemize{
#'   \item Source: \code{<source_dir>/<subset_N>/rna.mtx}
#'   \item Destination mat: \code{<dest_base>/<subset_N>/mat/rna.mtx}
#'   \item Destination mask: \code{<dest_base>/<subset_N>/mask/rna.csv}
#' }
#'
#' @note
#' If \code{rna_ref_presence.csv} is missing, a warning is printed but processing
#' continues. Folders without \code{rna.mtx} are skipped entirely.
#'
#' @examples
#' \dontrun{
#' # Use default paths
#' move_rna_files()
#'
#' # Custom paths
#' move_rna_files(
#'   source_dir = "data/raw/midas_output",
#'   dest_base = "data/processed/my_atlas"
#' )
#' }
#'
#' @seealso \code{\link{align_rna_features}} for generating the input files
#'
#' @export
move_rna_files <- function(
  source_dir = "data/raw/20260221_midas",
  dest_base = "data/processed/Atlas_20260221/data"
) {
  # 遍历源目录下的所有 subset* 文件夹
  subset_folders <- list.files(
    source_dir,
    pattern = "^subset",
    full.names = FALSE
  )
  
  subset_folders <- sort(subset_folders)
  
  processed_count <- 0
  skipped_count <- 0
  
  for (subset_folder in subset_folders) {
    source_path <- pj(source_dir, subset_folder)
    
    # 检查是否是目录
    if (!dir.exists(source_path)) {
      next
    }
    
    # 检查是否有 rna.mtx 文件
    rna_mtx_path <- pj(source_path, "rna.mtx")
    if (!file.exists(rna_mtx_path)) {
      prt("跳过", subset_folder, ": 没有 rna.mtx 文件")
      skipped_count <- skipped_count + 1
      next
    }
    
    # 创建目标目录结构
    mat_dir <- pj(dest_base, subset_folder, "mat")
    mask_dir <- pj(dest_base, subset_folder, "mask")
    
    mkdir(directory = mat_dir, remove_old = FALSE)
    mkdir(directory = mask_dir, remove_old = FALSE)
    
    # 移动 rna.mtx 到 mat 文件夹
    dest_mtx <- pj(mat_dir, "rna.mtx")
    file.rename(rna_mtx_path, dest_mtx)
    prt("移动: ", rna_mtx_path, " -> ", dest_mtx)
    
    # 移动并重命名 rna_ref_presence.csv 为 rna.csv 到 mask 文件夹
    rna_csv_source <- pj(source_path, "rna_ref_presence.csv")
    if (file.exists(rna_csv_source)) {
      dest_csv <- pj(mask_dir, "rna.csv")
      file.rename(rna_csv_source, dest_csv)
      prt("移动并重命名: ", rna_csv_source, " -> ", dest_csv)
    } else {
      prt("警告: ", subset_folder, " 中没有找到 rna_ref_presence.csv")
    }
    
    processed_count <- processed_count + 1
  }
  
  prt("\n处理完成！")
  prt("处理文件夹数: ", processed_count)
  prt("跳过文件夹数: ", skipped_count)
  
  invisible(list(
    processed = processed_count,
    skipped = skipped_count
  ))
}


library(data.table)
library(Matrix)

convert_csv_to_mtx <- function(base_path) {
  # 查找所有 mat 文件夹
  mat_dirs <- list.dirs(base_path, recursive = TRUE, full.names = TRUE)
  mat_dirs <- mat_dirs[grepl("/mat$", mat_dirs)]
  
  if (length(mat_dirs) == 0) {
    message("未找到任何 mat 文件夹")
    return(invisible(NULL))
  }
  
  # 处理每个 mat 文件夹
  for (mat_dir in mat_dirs) {
    # 查找该文件夹内所有 csv 文件
    csv_files <- list.files(mat_dir, pattern = "\\.csv$", full.names = TRUE)
    
    for (csv_file in csv_files) {
      message("处理: ", csv_file)
      
      # 读取 CSV
      data <- fread(csv_file, data.table = FALSE)
      
      # 提取行名
      rownames(data) <- data[[1]]
      data <- data[, -1]
      
      # 转换为稀疏矩阵
      sparse_matrix <- as(as.matrix(data), "dgCMatrix")
      
      # 生成输出文件名（同名 .mtx）
      mtx_file <- sub("\\.csv$", ".mtx", csv_file)
      
      # 写入 MTX
      writeMM(sparse_matrix, file = mtx_file)
      message("已保存: ", mtx_file)
    }
  }
  
  message("转换完成")
}


