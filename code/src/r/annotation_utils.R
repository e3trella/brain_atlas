# annotation_utils.R
# 基于注释表的细胞标注工具，适用于任意细胞类型的子聚类结果。
#
# 注释表要求：
#   必需列  cluster，以及至少一个 level 列，例如 level1、level2、level3
#   可选列  order、confidence、key_markers、literature_match、note，以及其他任意列
# 所有列按字符读入，cluster 必须与对象中聚类列的取值一致。

read_annotation <- function(file) {
  anno <- utils::read.csv(file, colClasses = "character", check.names = FALSE,
                          na.strings = "", strip.white = TRUE)
  if (!"cluster" %in% names(anno)) stop("注释表缺少 cluster 列: ", file)
  if (anyDuplicated(anno$cluster)) {
    stop("注释表中 cluster 有重复: ",
         paste(unique(anno$cluster[duplicated(anno$cluster)]), collapse = ", "))
  }
  if (length(grep("^level[0-9]+$", names(anno))) == 0) {
    stop("注释表至少需要一个 level 列，例如 level1: ", file)
  }
  if ("order" %in% names(anno)) {
    anno$order <- as.integer(anno$order)
    if (anyNA(anno$order)) stop("order 列存在空值或非整数")
    anno <- anno[order(anno$order), , drop = FALSE]
  }
  rownames(anno) <- NULL
  anno
}

# 将注释表写入 Seurat 对象的 metadata
#   cluster_col  对象中聚类结果所在的 metadata 列
#   prefix       新列的前缀，不同细胞类型请使用不同前缀，例如 mg_、ast_
#   cols         要写入的列，默认写入全部 level 列和 confidence 列
#   name         在 obj@misc$annotations 中记录本次注释时使用的名称
#   overwrite    同名列已存在时是否覆盖
apply_annotation <- function(obj, file, cluster_col, prefix,
                             cols = NULL, name = prefix, overwrite = FALSE) {
  anno <- read_annotation(file)
  md <- obj@meta.data
  if (!cluster_col %in% names(md)) stop("对象中不存在聚类列: ", cluster_col)

  if (is.null(cols)) {
    cols <- c(grep("^level[0-9]+$", names(anno), value = TRUE),
              intersect("confidence", names(anno)))
  }
  bad <- setdiff(cols, names(anno))
  if (length(bad) > 0) stop("注释表中没有这些列: ", paste(bad, collapse = ", "))

  cl <- as.character(md[[cluster_col]])
  missing_in_table <- setdiff(unique(cl), anno$cluster)
  missing_in_obj   <- setdiff(anno$cluster, unique(cl))
  if (length(missing_in_table) > 0) {
    stop("注释表缺少以下 Cluster: ", paste(sort(missing_in_table), collapse = ", "))
  }
  if (length(missing_in_obj) > 0) {
    warning("注释表中的以下 Cluster 不在对象中: ", paste(missing_in_obj, collapse = ", "))
  }

  new_cols <- paste0(prefix, cols)
  exists <- intersect(new_cols, names(md))
  if (length(exists) > 0 && !overwrite) {
    stop("以下列已存在，如需覆盖请设置 overwrite = TRUE: ", paste(exists, collapse = ", "))
  }

  idx <- match(cl, anno$cluster)
  for (i in seq_along(cols)) {
    values <- anno[[cols[i]]][idx]
    md[[new_cols[i]]] <- factor(values, levels = unique(anno[[cols[i]]]))
  }
  md[[cluster_col]] <- factor(cl, levels = intersect(anno$cluster, unique(cl)))
  obj@meta.data <- md

  if (is.null(obj@misc$annotations)) obj@misc$annotations <- list()
  obj@misc$annotations[[name]] <- list(file = normalizePath(file, mustWork = FALSE),
                                       cluster_col = cluster_col,
                                       columns = new_cols,
                                       date = Sys.time(),
                                       table = anno)
  message("已写入 ", length(new_cols), " 列: ", paste(new_cols, collapse = ", "))
  obj
}

# 从注释表中提取 key_markers，可用于 DotPlot
annotation_markers <- function(x) {
  anno <- if (is.character(x)) read_annotation(x) else x
  if (!"key_markers" %in% names(anno)) stop("注释表中没有 key_markers 列")
  m <- anno$key_markers[!is.na(anno$key_markers) & anno$key_markers != "none"]
  unique(trimws(unlist(strsplit(m, ";"))))
}

# 将子集对象的注释按细胞条码合并回全图谱对象
#   from_cols  子集对象中的列，例如 c("mg_level2", "mg_level3")
#   to_cols    全图谱对象中的目标列，例如 c("atlas_level2", "atlas_level3")
# 不在子集中的细胞保持原值，因此可以依次合并多个细胞类型。
merge_subset_annotation <- function(full, subset_obj, from_cols, to_cols = from_cols) {
  if (length(from_cols) != length(to_cols)) stop("from_cols 与 to_cols 长度不一致")
  sub_md <- subset_obj@meta.data
  bad <- setdiff(from_cols, names(sub_md))
  if (length(bad) > 0) stop("子集对象中没有这些列: ", paste(bad, collapse = ", "))

  idx <- match(colnames(full), colnames(subset_obj))
  hit <- !is.na(idx)
  not_found <- setdiff(colnames(subset_obj), colnames(full))
  if (length(not_found) > 0) {
    warning(length(not_found), " 个子集细胞不在全图谱对象中，已忽略")
  }

  md <- full@meta.data
  for (i in seq_along(from_cols)) {
    src <- sub_md[[from_cols[i]]]
    old <- if (to_cols[i] %in% names(md)) md[[to_cols[i]]] else factor(rep(NA_character_, nrow(md)))
    old_levels <- if (is.factor(old)) levels(old) else unique(stats::na.omit(as.character(old)))
    new_levels <- if (is.factor(src)) levels(src) else unique(stats::na.omit(as.character(src)))
    vals <- as.character(old)
    vals[hit] <- as.character(src)[idx[hit]]
    md[[to_cols[i]]] <- factor(vals, levels = unique(c(old_levels, new_levels)))
  }
  full@meta.data <- md
  message("已合并 ", sum(hit), " 个细胞的注释")
  full
}
