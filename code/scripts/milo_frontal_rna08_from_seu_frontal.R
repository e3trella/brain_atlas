# 从已经存在的 seu_frontal 对象开始，进行 Milo 细胞组成 / differential abundance 分析。
# 本脚本不负责创建 seu_frontal；运行前请先在 R 环境中准备好：
# seu_frontal <- subset(seu, subset = RL3_1 == "Frontal lobe")
# # KNN 构建邻接矩阵
# seu_frontal <- FindNeighbors(seu_frontal, reduction = "Xmidas_", dims = 1:32)

# # 聚类
# seu_frontal <- FindClusters(seu_frontal, resolution = 0.8)

# # 从 z_c 计算 UMAP
# seu_frontal <- RunUMAP(seu_frontal, reduction = "Xmidas_", dims = 1:32)

# 分析目标：
# 1. 将 seu_frontal 转换为 SingleCellExperiment / Milo 对象。
# 2. 使用 Xmidas_ 降维构建 Milo KNN graph 和 neighborhoods。
# 3. 以 batch 作为样本单位，对每个 neighborhood 计数。
# 4. 使用 DISEASE_NEW 做疾病组 vs Healthy 的 Milo DA 检验。
# 5. 使用 RNA_snn_res.0.8 注释每个 Milo neighborhood。
# 6. 生成按 RNA_snn_res.0.8 cluster 汇总的最终结果表和图。

library(Seurat)
library(SingleCellExperiment)
library(miloR)
library(dplyr)
library(ggplot2)

# Step 1: 检查输入对象和后续分析所需的 metadata / reduction 是否存在。
# 如果缺少关键列或降维结果，脚本会在这里停止，避免后面产生难以解释的错误。
print("[Input validation] started")

stopifnot(exists("seu_frontal"))
stopifnot(inherits(seu_frontal, "Seurat"))
stopifnot("batch" %in% colnames(seu_frontal@meta.data))
stopifnot("DISEASE_NEW" %in% colnames(seu_frontal@meta.data))
stopifnot("PUBLICATION" %in% colnames(seu_frontal@meta.data))
stopifnot("RNA_snn_res.0.8" %in% colnames(seu_frontal@meta.data))
stopifnot("Xmidas_" %in% Seurat::Reductions(seu_frontal))

print("[Input validation] completed")

# Step 2: 设定 Milo graph 使用的低维空间维度。
# 这里沿用前面分析中使用的 Xmidas_ 前 32 维。
selected_dims <- 1:32

# Step 3: 将 Seurat 对象转换为 SingleCellExperiment。
# Milo 以 SingleCellExperiment 为输入，因此需要先转换对象格式。
print(paste("[SingleCellExperiment conversion] started, total", ncol(seu_frontal), "cells"))

sce_frontal <- as.SingleCellExperiment(seu_frontal)

print("[SingleCellExperiment conversion] completed")

# Step 4: 将 Seurat 里的 Xmidas_ embedding 复制到 SingleCellExperiment 的 reducedDim 中。
# Seurat 转换后 reduction 名称可能会改变；这里显式创建一个名为 Xmidas 的 reducedDim，
# 后续所有 Milo graph、neighborhood distance 都使用这个相同的 reducedDim。
print("[Xmidas reducedDim transfer] started")

reducedDim(sce_frontal, "Xmidas") <- Embeddings(seu_frontal, "Xmidas_")

print("[Xmidas reducedDim transfer] completed")

# Step 5: 创建 Milo 对象。
# 这个对象会保存 graph、neighborhoods、样本计数矩阵和 DA 检验需要的信息。
print("[Milo object preparation] started")

milo_frontal <- Milo(sce_frontal)

print("[Milo object preparation] completed")

# Step 6: 固定随机种子。
# makeNhoods() 中 refined sampling 涉及随机抽样；固定 seed 可以让结果更容易复现。
set.seed(7186)

# Step 7: 使用 Xmidas 低维空间构建 KNN graph。
# k = 30 表示每个细胞连接到 30 个近邻；d = 32 表示使用 Xmidas 的前 32 维。
print(paste("[Milo graph construction] started, total", ncol(seu_frontal), "cells"))

milo_frontal <- buildGraph(
  milo_frontal,
  k = 30,
  d = length(selected_dims),
  reduced.dim = "Xmidas"
)

print("[Milo graph construction] completed")

# Step 8: 在 KNN graph 上抽样并构建 Milo neighborhoods。
# prop = 0.1 表示大约以 10% 的细胞作为候选 neighborhood index。
# refined = TRUE 会在低维空间中 refine neighborhood sampling，减少冗余 neighborhoods。
print("[Milo neighborhood construction] started")

milo_frontal <- makeNhoods(
  milo_frontal,
  prop = 0.1,
  k = 30,
  d = length(selected_dims),
  refined = TRUE,
  reduced_dims = "Xmidas"
)

print("[Milo neighborhood construction] completed")

# Step 9: 对每个 Milo neighborhood 按 batch 统计细胞数。
# 这里 batch 被当作样本 / biological replicate 单位。
# 后续 DA 模型不是按单个细胞建模，而是按 neighborhood x batch 的计数矩阵建模。
print("[Milo cell counting] started")

milo_frontal <- countCells(
  milo_frontal,
  meta.data = as.data.frame(colData(milo_frontal)),
  samples = "batch"
)

print("[Milo cell counting] completed")

# Step 10: 构建样本级 design table。
# design_frontal 每一行对应一个 batch，包含该 batch 的疾病状态 DISEASE_NEW 和数据来源 PUBLICATION。
# Milo 的 design.df 必须是样本级 metadata，而不是细胞级 metadata。
print("[Milo sample design] started")

design_frontal <- as.data.frame(colData(milo_frontal)) |>
  distinct(batch, DISEASE_NEW, PUBLICATION) |>
  arrange(batch)

rownames(design_frontal) <- design_frontal$batch

# Step 11: 检查每个 batch 是否只对应一个疾病状态。
# 如果同一个 batch 同时属于多个 DISEASE_NEW，说明 batch 不能直接作为样本单位。
batch_disease_check <- as.data.frame(colData(milo_frontal)) |>
  distinct(batch, DISEASE_NEW) |>
  count(batch, name = "n_disease") |>
  filter(n_disease > 1)

stopifnot(nrow(batch_disease_check) == 0)

# Step 12: 检查每个 batch 是否只对应一个 PUBLICATION。
# 这里不把 PUBLICATION 放入最终模型，只做一致性检查，避免样本 metadata 混乱。
batch_publication_check <- as.data.frame(colData(milo_frontal)) |>
  distinct(batch, PUBLICATION) |>
  count(batch, name = "n_publication") |>
  filter(n_publication > 1)

stopifnot(nrow(batch_publication_check) == 0)

# Step 13: 检查 Healthy 是否存在。
# 后续所有 contrast 都以 Healthy 为参考组；没有 Healthy 就无法运行疾病 vs Healthy 比较。
stopifnot("Healthy" %in% design_frontal$DISEASE_NEW)

# Step 14: 设置疾病因子的顺序。
# 把 Healthy 放在第一位，使模型默认以 Healthy 作为 reference level。
# 其他疾病按 batch 数量从多到少排列。
disease_levels <- c(
  "Healthy",
  setdiff(
    design_frontal |>
      count(DISEASE_NEW, name = "n_batches") |>
      arrange(desc(n_batches)) |>
      pull(DISEASE_NEW),
    "Healthy"
  )
)

design_frontal$DISEASE_NEW <- factor(design_frontal$DISEASE_NEW, levels = disease_levels)

# Step 15: 创建用于 Milo / limma contrast 的安全疾病名称。
# 比如 COVID-19 含有连字符，不是合法的 R 模型系数名；make.names() 会将其变成 COVID.19。
# DISEASE_NEW 保留原始标签，DISEASE_MILO 只用于建模。
design_frontal$DISEASE_MILO <- factor(
  make.names(as.character(design_frontal$DISEASE_NEW)),
  levels = make.names(levels(design_frontal$DISEASE_NEW))
)

# Step 16: 建立安全疾病名称和原始疾病名称之间的对应表。
# 后面把结果中的 disease_safe 转回原始疾病名称，方便阅读。
disease_milo_lookup <- design_frontal |>
  distinct(DISEASE_NEW, DISEASE_MILO) |>
  arrange(DISEASE_NEW)

# Step 17: 创建所有疾病 vs Healthy 的 contrast 列表。
# 每个疾病单独和 Healthy 比较，不做一个总体多组检验。
disease_contrasts_milo <- setdiff(levels(design_frontal$DISEASE_MILO), "Healthy")

print("[Milo sample design] completed")

# Step 18: 计算 neighborhood 之间的距离。
# 这些距离用于 Milo 的 spatial FDR correction。
# reduced.dim 必须和前面构图时使用的 reducedDim 保持一致。
print("[Milo neighborhood distance] started")

milo_frontal <- calcNhoodDistance(
  milo_frontal,
  d = length(selected_dims),
  reduced.dim = "Xmidas"
)

print("[Milo neighborhood distance] completed")

# Step 19: 对每一个疾病 vs Healthy contrast 运行 Milo DA 检验。
# testNhoods() 使用 neighborhood x batch 的计数矩阵，并用 DISEASE_MILO 作为解释变量。
# 每个 result 会额外添加 contrast 和 disease 两列，方便后续合并和汇总。
print(paste("[Milo DA test] started, total", length(disease_contrasts_milo), "contrasts"))

milo_da_results <- lapply(disease_contrasts_milo, function(disease_safe) {
  contrast <- paste0("DISEASE_MILO", disease_safe)

  result <- testNhoods(
    milo_frontal,
    design = ~ DISEASE_MILO,
    design.df = design_frontal,
    model.contrasts = contrast,
    reduced.dim = "Xmidas"
  )

  disease_original <- as.character(
    disease_milo_lookup$DISEASE_NEW[disease_milo_lookup$DISEASE_MILO == disease_safe]
  )

  result$contrast <- paste(disease_original, "vs Healthy")
  result$disease <- disease_original

  result
})

names(milo_da_results) <- disease_contrasts_milo

print("[Milo DA test] completed")

# Step 20: 汇总每个疾病比较中显著 DA neighborhoods 的数量。
# SpatialFDR < 0.05 被视为显著。
# logFC > 0 表示该 neighborhood 在疾病组相对 Healthy 增加；logFC < 0 表示减少。
print("[Milo DA summary] started")

milo_da_summary <- bind_rows(milo_da_results) |>
  group_by(contrast, disease) |>
  summarise(
    n_neighborhoods = n(),
    n_spatial_fdr_0_05 = sum(SpatialFDR < 0.05, na.rm = TRUE),
    n_increased = sum(SpatialFDR < 0.05 & logFC > 0, na.rm = TRUE),
    n_decreased = sum(SpatialFDR < 0.05 & logFC < 0, na.rm = TRUE),
    min_spatial_fdr = min(SpatialFDR, na.rm = TRUE),
    .groups = "drop"
  ) |>
  arrange(desc(n_spatial_fdr_0_05))

print("[Milo DA summary] completed")

# Step 21: 使用 RNA_snn_res.0.8 注释每个 Milo neighborhood。
# annotateNhoods() 会给每个 neighborhood 赋予其中占比最高的 RNA_snn_res.0.8 cluster。
# 同时生成 RNA_snn_res.0.8_fraction，表示主导 cluster 在该 neighborhood 内的细胞比例。
print("[RNA_snn_res.0.8 annotation] started")

milo_da_results_rna08 <- lapply(milo_da_results, function(result) {
  annotateNhoods(
    milo_frontal,
    result,
    coldata_col = "RNA_snn_res.0.8"
  )
})

# Step 22: 合并所有疾病 contrast 的注释后 Milo DA 结果。
# milo_da_rna08 是后续按 cluster 汇总的基础表。
milo_da_rna08 <- bind_rows(milo_da_results_rna08)

print("[RNA_snn_res.0.8 annotation] completed")

# Step 23: 生成最终的 RNA_snn_res.0.8 cluster 层级汇总表。
# milo_rna08_summary：每个 contrast 只保留显著 neighborhood 数量最多的前 10 个 cluster-direction 组合。
# 这个表适合快速查看每个疾病比较中变化最明显的 cluster。
print("[RNA_snn_res.0.8 summary] started")

milo_rna08_summary <- milo_da_rna08 |>
  filter(SpatialFDR < 0.05) |>
  mutate(direction = ifelse(logFC > 0, "increased", "decreased")) |>
  count(contrast, disease, RNA_snn_res.0.8, direction, name = "n_sig_nhoods") |>
  group_by(contrast) |>
  slice_max(n_sig_nhoods, n = 10, with_ties = FALSE) |>
  ungroup() |>
  arrange(contrast, desc(n_sig_nhoods))

# Step 24: 生成完整的 RNA_snn_res.0.8 cluster 层级汇总表。
# milo_rna08_totals 保留所有显著 cluster-direction 组合，不只限前 10。
# 如果后续要导出或进一步筛选，优先使用这个对象。
milo_rna08_totals <- milo_da_rna08 |>
  filter(SpatialFDR < 0.05) |>
  mutate(direction = ifelse(logFC > 0, "increased", "decreased")) |>
  count(contrast, disease, RNA_snn_res.0.8, direction, name = "n_sig_nhoods") |>
  arrange(desc(n_sig_nhoods))

print("[RNA_snn_res.0.8 summary] completed")

# Step 25: 准备用于可视化的数据。
# 每个 contrast 取显著 neighborhood 数量最多的前 8 个 cluster-direction 组合。
# increased 画成正值，decreased 画成负值，方便在同一张图中比较方向。
print("[RNA_snn_res.0.8 plot] started")

milo_rna08_plot_data <- milo_rna08_totals |>
  group_by(contrast) |>
  slice_max(n_sig_nhoods, n = 8, with_ties = FALSE) |>
  ungroup() |>
  mutate(
    signed_n = ifelse(direction == "increased", n_sig_nhoods, -n_sig_nhoods),
    cluster = paste0("cluster ", RNA_snn_res.0.8)
  )

# Step 26: 绘制最终 summary plot。
# x 轴为显著 neighborhoods 数量；负数表示疾病组相对 Healthy decreased，正数表示 increased。
# 每个 facet 是一个疾病 vs Healthy contrast。
milo_rna08_plot <- ggplot(
  milo_rna08_plot_data,
  aes(x = signed_n, y = reorder(cluster, signed_n), fill = direction)
) +
  geom_col() +
  facet_wrap(~ contrast, scales = "free_y") +
  labs(
    x = "Significant neighborhoods (negative = decreased)",
    y = "RNA_snn_res.0.8 cluster",
    fill = "Direction"
  )

print(milo_rna08_plot)

print("[RNA_snn_res.0.8 plot] completed")

# Step 27: 输出关键结果对象。
# 脚本运行结束后，以下对象会留在当前 R 环境：
# - sce_frontal: 从 seu_frontal 转换得到的 SingleCellExperiment。
# - milo_frontal: Milo object，包含 graph、neighborhoods、counts 和 distance。
# - design_frontal: 样本级 design table。
# - milo_da_results: 每个疾病 vs Healthy 的 Milo DA 结果 list。
# - milo_da_summary: 每个 contrast 的显著 DA neighborhood 总数汇总。
# - milo_da_results_rna08: 加入 RNA_snn_res.0.8 注释后的 DA 结果 list。
# - milo_da_rna08: 合并后的所有注释 DA 结果。
# - milo_rna08_summary: 每个 contrast 前 10 个显著 cluster-direction 汇总。
# - milo_rna08_totals: 完整显著 cluster-direction 汇总。
# - milo_rna08_plot: 最终可视化图对象。
print("[Milo frontal RNA_snn_res.0.8 analysis] completed")

milo_da_summary
milo_rna08_summary
milo_rna08_totals
