run_midas_clustering <- function(
    seu,
    reduction = "Xmidas_",
    dims = 1:32,
    resolution = 0.8) {
  seu <- Seurat::FindNeighbors(seu, reduction = reduction, dims = dims)
  seu <- Seurat::FindClusters(seu, resolution = resolution)
  seu <- Seurat::RunUMAP(seu, reduction = reduction, dims = dims)
  seu
}
