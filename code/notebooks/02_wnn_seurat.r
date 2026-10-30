## =====================================================================
## 用 Seurat 与 Signac 对同一组细胞的 RNA 与 ATAC 矩阵做 WNN 联合降维
## 前提：两个矩阵已经没有批次效应，因此不做跨批次整合
## =====================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(Matrix)
  library(ggplot2)
  library(GenomicRanges)
  library(GenomeInfoDb)
  library(EnsDb.Hsapiens.v86)   # Ensembl 86，基于 GRCh38
})
options(future.globals.maxSize = 50 * 1024^3)

## ---------------- 参数 ----------------
rna_file  <- "data/processed/Atlas_20260721/output/wnn/pure_astro/pure_astro_rna.h5ad"
atac_file <- "data/processed/Atlas_20260721/output/wnn/pure_astro/pure_astro_atac.h5ad"
out_dir   <- "data/processed/Atlas_20260721/output/wnn/pure_astro"
out_rds_name <- "pure_astro_wnn_seurat.rds"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## 根据 subset_h5ad.py 的诊断输出设置矩阵的数值性质
##   "counts"     整数原始计数
##   "normalized" 已归一化或经 MIDAS 校正的连续值
RNA_INPUT  <- "counts"
ATAC_INPUT <- "counts"

n_pcs    <- 50
rna_dims <- 1:30
lsi_dims <- 3:30      # 依据 DepthCor 结果调整，剔除与深度高度相关的成分

nfeatures              <- 3000
top_feature_min_cutoff <- "q5"
prune_snn              <- 1/20
k_nn                   <- 20
cluster_resolution     <- 0.5
cluster_algorithm      <- 3
cluster_n_start        <- 10

param_names <- c(
  "rna_file", "atac_file", "out_dir", "out_rds_name",
  "RNA_INPUT", "ATAC_INPUT", "n_pcs", "rna_dims", "lsi_dims",
  "nfeatures", "top_feature_min_cutoff", "prune_snn", "k_nn",
  "cluster_resolution", "cluster_algorithm", "cluster_n_start"
)
params <- mget(param_names, envir = environment())
cat("可调参数：\n")
for (nm in names(params)) {
  cat(sprintf("  %-24s = %s\n", nm, paste(params[[nm]], collapse = ", ")))
}

## ---------------- 读取 ----------------
cat("[1/6] 读取 h5ad ...\n")
rna_tmp  <- schard::h5ad2seurat(rna_file,  layers = c(counts = "counts"))
atac_tmp <- schard::h5ad2seurat(atac_file, layers = c(counts = "counts"))

rna_mat  <- LayerData(rna_tmp,  assay = DefaultAssay(rna_tmp),  layer = "counts")
atac_mat <- LayerData(atac_tmp, assay = DefaultAssay(atac_tmp), layer = "counts")
if (!inherits(rna_mat,  "dgCMatrix")) rna_mat  <- as(rna_mat,  "dgCMatrix")
if (!inherits(atac_mat, "dgCMatrix")) atac_mat <- as(atac_mat, "dgCMatrix")

stopifnot(identical(colnames(rna_mat), colnames(atac_mat)))
cat("细胞数：", ncol(rna_mat),
    "  基因数：", nrow(rna_mat),
    "  peak 数：", nrow(atac_mat), "\n")

## ---------------- 构建对象 ----------------
cat("[2/6] 构建 Seurat 对象与 ATAC assay ...\n")
obj <- CreateSeuratObject(counts = rna_mat, assay = "RNA",
                          meta.data = rna_tmp@meta.data[colnames(rna_mat), , drop = FALSE])

## ---- peak 坐标解析与过滤 ----
## 行名格式为 chr1:9925-10490，按冒号与连字符拆分
cat("  - 解析 peak 坐标 ...\n")
peak_gr <- StringToGRanges(rownames(atac_mat), sep = c(":", "-"))
cat("peak 所在染色体：\n"); print(table(seqnames(peak_gr)))

## 仅保留常规染色体，去除 KI270、GL000 一类未定位片段与线粒体
keep_chr <- paste0("chr", c(1:22, "X", "Y"))
cat("  keep_chr =", paste(keep_chr, collapse = ", "), "\n")
keep <- as.character(seqnames(peak_gr)) %in% keep_chr
cat("保留 peak：", sum(keep), " / ", length(keep), "\n")
atac_mat <- atac_mat[keep, ]
peak_gr  <- peak_gr[keep]

## ---- 基因注释 ----
cat("  - 构建基因注释 ...\n")
anno_file <- "parameters/annotation_hg38_EnsDb.v86.rds"
if (file.exists(anno_file)) {
  annotation <- readRDS(anno_file)
} else {
  annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
  seqlevelsStyle(annotation) <- "UCSC"
  genome(annotation) <- "hg38"
}

## ---- ChromatinAssay，不提供 fragments ----
cat("  - 创建 ChromatinAssay ...\n")
atac_assay <- CreateChromatinAssay(
  counts     = atac_mat,
  ranges     = peak_gr,
  genome     = "hg38",
  annotation = annotation,
  min.cells  = 0,
  min.features = 0
)
obj[["ATAC"]] <- atac_assay
obj$nCount_ATAC   <- colSums(atac_mat)
obj$nFeature_ATAC <- colSums(atac_mat > 0)

## 检查坐标是否正确写入
print(granges(obj[["ATAC"]])[1:5])
cat("与基因体或启动子重叠的 peak 比例：",
    mean(overlapsAny(granges(obj[["ATAC"]]),
                     promoters(annotation, upstream = 2000, downstream = 0)) |
         overlapsAny(granges(obj[["ATAC"]]), annotation)), "\n")

## 保留 MIDAS 的嵌入，便于与 WNN 结果对照
for (r in Reductions(rna_tmp)) obj[[paste0("midas_", r)]] <- rna_tmp[[r]]

rm(rna_tmp, atac_tmp, rna_mat, atac_mat); gc()

## ---------------- RNA 降维 ----------------
cat("[3/6] RNA 降维：NormalizeData / FindVariableFeatures / ScaleData / RunPCA ...\n")
DefaultAssay(obj) <- "RNA"
if (RNA_INPUT == "counts") {
  obj <- NormalizeData(obj)
} else {
  LayerData(obj, assay = "RNA", layer = "data") <- LayerData(obj, assay = "RNA", layer = "counts")
}
obj <- FindVariableFeatures(obj, nfeatures = nfeatures)
obj <- ScaleData(obj, features = VariableFeatures(obj))
obj <- RunPCA(obj, npcs = n_pcs, reduction.name = "pca", verbose = FALSE)
ElbowPlot(obj, ndims = n_pcs, reduction = "pca")
ggsave(file.path(out_dir, "rna_elbow.pdf"), width = 5, height = 4)

## ---------------- ATAC 降维 ----------------
cat("[4/6] ATAC 降维：RunTFIDF / FindTopFeatures / RunSVD ...\n")
DefaultAssay(obj) <- "ATAC"
if (ATAC_INPUT == "counts") {
  obj <- RunTFIDF(obj)
} else {
  ## 已归一化的连续值不再做 TF-IDF，直接作为 data 层进入 SVD
  obj <- SetAssayData(obj, assay = "ATAC", slot = "data",
                      new.data = GetAssayData(obj, assay = "ATAC", slot = "counts"))
}
obj <- FindTopFeatures(obj, min.cutoff = top_feature_min_cutoff)
obj <- RunSVD(obj, n = n_pcs, reduction.name = "lsi")

DepthCor(obj, reduction = "lsi", n = n_pcs)
ggsave(file.path(out_dir, "lsi_depthcor.pdf"), width = 6, height = 4)

## ---------------- 单模态 UMAP，用于对照 ----------------
cat("[5/6] 单模态 UMAP (RNA 与 ATAC) ...\n")
obj <- RunUMAP(obj, reduction = "pca", dims = rna_dims,
               reduction.name = "umap.rna", reduction.key = "rnaUMAP_")
obj <- RunUMAP(obj, reduction = "lsi", dims = lsi_dims,
               reduction.name = "umap.atac", reduction.key = "atacUMAP_")

## ---------------- WNN 联合降维 ----------------
cat("[6/6] WNN 联合降维：FindMultiModalNeighbors / RunUMAP / FindClusters ...\n")
cat("  - FindMultiModalNeighbors ...\n")
obj <- FindMultiModalNeighbors(
  obj,
  reduction.list = list("pca", "lsi"),
  prune.SNN = prune_snn, 
  k.nn = k_nn,
  dims.list = list(rna_dims, lsi_dims),
  modality.weight.name = c("RNA.weight", "ATAC.weight")
)
cat("  - RunUMAP (wnn.umap) ...\n")
obj <- RunUMAP(obj, nn.name = "weighted.nn",
               reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
cat("  - FindClusters (wsnn, resolution = 0.5) ...\n")
obj <- FindClusters(obj, graph.name = "wsnn", algorithm = cluster_algorithm, n.start = cluster_n_start,
                    resolution = cluster_resolution, verbose = TRUE)

## ---------------- 可视化与检查 ----------------

cat("WNN 完成，导出结果 ...\n")
if ("midas_leiden" %in% colnames(obj@meta.data)) {
  write.csv(as.data.frame.matrix(table(WNN = obj$seurat_clusters, MIDAS = obj$midas_leiden)),
            file.path(out_dir, "wnn_vs_midas_clusters.csv"))
}

saveRDS(obj, file.path(out_dir, out_rds_name))
cat("已保存：", file.path(out_dir, out_rds_name), "\n")