# 对 pure_imm WNN 对象中的小胶质细胞做 Slingshot 拟时序分析。
# 分析目标：
# 1. 在 wnn.umap 的 2 维坐标上，以 mg_level3 作为簇标签推断谱系。
# 2. 起点簇设为 C1_Homeo_CX3CR1，不指定终点簇，由数据自己推出分叉结构。
# 3. 输出每个细胞在每条谱系上的拟时序和权重，并画出谱系曲线叠加的 UMAP。

library(Seurat)
library(slingshot)
library(ggplot2)

in_rds <- "data/processed/Atlas_20260721/output/wnn/pure_imm/candidate_pure_imm_wnn_seurat_formonocle.rds"
out_dir <- "data/processed/Atlas_20260721/output/wnn/pure_imm/slingshot"
table_dir <- "results/tables/slingshot_pure_imm_wnn"
plot_dir <- "results/plots/slingshot_pure_imm_wnn"
for (d in c(out_dir, table_dir, plot_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

cluster_col <- "mg_level3"
reduction <- "wnn.umap"
start_clus <- "C1_Homeo_CX3CR1"

set.seed(1)

# Step 1: 读取对象，只取降维坐标和 meta.data，然后释放对象。
print("[Input] started")

seu <- readRDS(in_rds)
stopifnot(cluster_col %in% colnames(seu@meta.data))
stopifnot(reduction %in% Seurat::Reductions(seu))

emb <- Embeddings(seu, reduction)
meta <- seu@meta.data
rm(seu)
gc()

cluster_labels <- droplevels(factor(meta[[cluster_col]]))
stopifnot(start_clus %in% levels(cluster_labels))
print(table(cluster_labels))

print("[Input] completed")

# Step 2: 运行 Slingshot。
# 传入矩阵和簇标签时返回 SlingshotDataSet，不需要构建 SingleCellExperiment。
print("[Slingshot] started")

sds <- slingshot(emb, clusterLabels = cluster_labels, start.clus = start_clus)

print(slingLineages(sds))

print("[Slingshot] completed")

# Step 3: 保存拟时序和曲线权重。NA 表示该细胞不在这条谱系上。
pt <- slingPseudotime(sds)
wt <- slingCurveWeights(sds)
colnames(wt) <- paste0(colnames(pt), "_weight")

pt_table <- data.frame(
  cell = rownames(emb),
  cluster = cluster_labels,
  pt,
  wt,
  check.names = FALSE
)
write.csv(pt_table, file.path(table_dir, "pseudotime.csv"), row.names = FALSE)
saveRDS(sds, file.path(out_dir, "slingshot_sds.rds"))

# Step 4: 画图。每条谱系一张拟时序着色的 UMAP，另有一张按簇着色并叠加所有曲线。
print("[Plot] started")

df <- data.frame(UMAP1 = emb[, 1], UMAP2 = emb[, 2], cluster = cluster_labels)
curves <- do.call(rbind, lapply(seq_along(slingCurves(sds)), function(i) {
  cv <- slingCurves(sds)[[i]]
  data.frame(
    UMAP1 = cv$s[cv$ord, 1],
    UMAP2 = cv$s[cv$ord, 2],
    lineage = colnames(pt)[i]
  )
}))

p_cluster <- ggplot(df, aes(UMAP1, UMAP2, color = cluster)) +
  geom_point(size = 0.1, alpha = 0.5) +
  geom_path(data = curves, aes(group = lineage), color = "black", linewidth = 0.8) +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  theme_classic()
ggsave(file.path(plot_dir, "umap_cluster_curves.png"), p_cluster, width = 7, height = 6, dpi = 200)

for (i in seq_len(ncol(pt))) {
  lin <- colnames(pt)[i]
  df$pseudotime <- pt[, i]
  p <- ggplot(df[order(!is.na(df$pseudotime)), ], aes(UMAP1, UMAP2, color = pseudotime)) +
    geom_point(size = 0.1) +
    geom_path(data = curves[curves$lineage == lin, ], color = "black", linewidth = 0.8, inherit.aes = FALSE,
              aes(UMAP1, UMAP2)) +
    scale_color_viridis_c(na.value = "grey85") +
    ggtitle(lin) +
    theme_classic()
  ggsave(file.path(plot_dir, paste0("umap_pseudotime_", lin, ".png")), p, width = 6, height = 5, dpi = 200)
}

print("[Plot] completed")
