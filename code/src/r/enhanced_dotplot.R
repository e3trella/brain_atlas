# =============================================================================
# EnhancedDotPlot — 优化版 Seurat DotPlot
# 支持：基因聚类 / 细胞亚群聚类 / 自定义顺序 / scale & normalize
# =============================================================================
# 依赖：Seurat, ggplot2, dplyr, tidyr
# 安装：install.packages(c("ggplot2","dplyr","tidyr"))
#        BiocManager::install("Seurat")
# =============================================================================

EnhancedDotPlot <- function(
    seurat_obj,
    features,                        # character vector，marker 基因列表
    group_by        = "seurat_clusters", # 细胞分组 meta.data 列名

    # --- 排序控制 ---
    cluster_genes   = TRUE,          # 基因按 hierarchical clustering 排序
    cluster_cells   = TRUE,          # 细胞亚群按基因表达 hierarchical clustering 排序
    gene_order      = NULL,          # 手动指定基因顺序（覆盖 cluster_genes）
    cell_order      = NULL,          # 手动指定细胞亚群顺序（覆盖 cluster_cells）

    # --- 表达值处理 ---
    scale_data      = TRUE,          # 是否对 avg_exp 做 z-score（跨亚群）
    normalize_method = "log1p",      # "log1p" | "zscore" | "none"
                                     #   log1p : log1p(avg_exp)
                                     #   zscore: scale across groups per gene
                                     #   none  : raw avg_exp

    # --- 视觉参数 ---
    dot_scale       = 6,             # 最大点半径（mm）
    col_low         = "#313695",
    col_mid         = "#ffffbf",
    col_high        = "#d73027",
    dot_min         = 0,             # pct.exp 下限过滤（0~1）
    theme_base      = theme_bw()
) {

  # ── 0. 包检查 ──────────────────────────────────────────────────────────────
  for (pkg in c("Seurat", "ggplot2", "dplyr", "tidyr")) {
    if (!requireNamespace(pkg, quietly = TRUE))
      stop(sprintf("请先安装包：%s", pkg))
  }
  library(ggplot2); library(dplyr); library(tidyr)

  # ── 1. 过滤不存在的基因 ────────────────────────────────────────────────────
  all_genes <- rownames(seurat_obj)
  missing   <- setdiff(features, all_genes)
  if (length(missing) > 0)
    warning("以下基因在 SeuratObject 中不存在，已跳过：",
            paste(missing, collapse = ", "))
  features <- intersect(features, all_genes)
  if (length(features) == 0) stop("没有有效的基因，请检查 features 参数。")

  # ── 2. 提取点图数据（avg_exp + pct.exp）──────────────────────────────────
  dp_data <- Seurat::DotPlot(seurat_obj,
                              features = features,
                              group.by = group_by)$data
  # 列名标准化
  colnames(dp_data) <- gsub("\\.", "_", colnames(dp_data))
  # dp_data 列：avg_exp, pct_exp, avg_exp_scaled, features_plot, id

  # ── 3. Normalize avg_exp ──────────────────────────────────────────────────
  dp_data <- dp_data %>%
    group_by(features_plot) %>%
    mutate(
      avg_exp_norm = switch(
        normalize_method,
        "log1p"  = log1p(avg_exp),
        "zscore" = {
          s <- sd(avg_exp, na.rm = TRUE)
          if (is.na(s) || s == 0) rep(0, n()) else (avg_exp - mean(avg_exp)) / s
        },
        "none"   = avg_exp,
        stop("normalize_method 须为 'log1p'、'zscore' 或 'none'")
      )
    ) %>%
    ungroup()

  # ── 4. Scale across groups（可选，在 normalize 之后再做）──────────────────
  if (scale_data) {
    dp_data <- dp_data %>%
      group_by(features_plot) %>%
      mutate(avg_exp_plot = {
        s <- sd(avg_exp_norm, na.rm = TRUE)
        if (is.na(s) || s == 0) rep(0, n()) else
          (avg_exp_norm - mean(avg_exp_norm, na.rm = TRUE)) / s
      }) %>%
      ungroup()
  } else {
    dp_data$avg_exp_plot <- dp_data$avg_exp_norm
  }

  # ── 5. 过滤低 pct.exp 的点 ────────────────────────────────────────────────
  dp_data <- dp_data %>% filter(pct_exp >= dot_min * 100)

  # ── 6. 构建宽矩阵（用于聚类）─────────────────────────────────────────────
  mat_wide <- dp_data %>%
    select(features_plot, id, avg_exp_plot) %>%
    pivot_wider(names_from = id, values_from = avg_exp_plot, values_fill = 0) %>%
    tibble::column_to_rownames("features_plot") %>%
    as.matrix()

  # ── 7. 基因顺序 ───────────────────────────────────────────────────────────
  if (!is.null(gene_order)) {
    gene_order <- intersect(gene_order, rownames(mat_wide))
    gene_levels <- gene_order
  } else if (cluster_genes && nrow(mat_wide) > 1) {
    hc_genes   <- hclust(dist(mat_wide), method = "ward.D2")
    gene_levels <- rownames(mat_wide)[hc_genes$order]
  } else {
    gene_levels <- rownames(mat_wide)
  }

  # ── 8. 细胞亚群顺序 ───────────────────────────────────────────────────────
  if (!is.null(cell_order)) {
    cell_order  <- intersect(cell_order, colnames(mat_wide))
    cell_levels <- cell_order
  } else if (cluster_cells && ncol(mat_wide) > 1) {
    hc_cells    <- hclust(dist(t(mat_wide)), method = "ward.D2")
    cell_levels <- colnames(mat_wide)[hc_cells$order]
  } else {
    cell_levels <- colnames(mat_wide)
  }

  # ── 9. 设置 factor 顺序 ───────────────────────────────────────────────────
  dp_data <- dp_data %>%
    mutate(
      features_plot = factor(features_plot, levels = gene_levels),
      id            = factor(id,            levels = cell_levels)
    )

  # ── 10. 颜色范围 ──────────────────────────────────────────────────────────
  exp_vals   <- dp_data$avg_exp_plot
  color_mid  <- mean(range(exp_vals, na.rm = TRUE))

  # ── 11. 绘图 ──────────────────────────────────────────────────────────────
  p <- ggplot(dp_data,
              aes(x     = features_plot,
                  y     = id,
                  size  = pct_exp,
                  color = avg_exp_plot)) +
    geom_point() +
    scale_size_continuous(
      name   = "% Expressed",
      range  = c(0, dot_scale),
      limits = c(0, 100)
    ) +
    scale_color_gradient2(
      name     = ifelse(scale_data, "Scaled Expr.", "Avg. Expr."),
      low      = col_low,
      mid      = col_mid,
      high     = col_high,
      midpoint = color_mid
    ) +
    theme_base +
    theme(
      axis.text.x  = element_text(angle = 45, hjust = 1, vjust = 1),
      axis.title   = element_blank(),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.3)
    ) +
    guides(
      size  = guide_legend(override.aes = list(color = "grey40")),
      color = guide_colorbar(barwidth = 1, barheight = 6)
    )

  return(p)
}


# =============================================================================
# 使用示例
# =============================================================================

# library(Seurat)
# library(ggplot2)
#
# markers <- c("CD3D","CD3E","CD8A","CD4",
#              "MS4A1","CD79A","GNLY","NKG7",
#              "LYZ","CST3","PPBP")
#
# # 基本用法
# p <- EnhancedDotPlot(pbmc, features = markers)
# p
#
# # 自定义细胞亚群顺序 + 不聚类基因
# p2 <- EnhancedDotPlot(
#   pbmc,
#   features      = markers,
#   cell_order    = c("0","1","2","3","4","5","6","7","8"),
#   cluster_genes = FALSE
# )
#
# # 不 scale，用 log1p normalize，叠加标题
# p3 <- EnhancedDotPlot(
#   pbmc,
#   features          = markers,
#   scale_data        = FALSE,
#   normalize_method  = "log1p"
# ) + ggtitle("My Enhanced DotPlot")