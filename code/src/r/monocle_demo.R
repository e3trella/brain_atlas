library(Seurat)
library(SeuratWrappers)
library(monocle3)
library(igraph)
library(ggplot2)

## --- Seurat → Monocle3 ---
DefaultAssay(seu_f) <- "RNA"
if (!"UMAP" %in% Reductions(seu_f)) seu_f[["UMAP"]] <- seu_f[["umap"]]

cds <- as.cell_data_set(seu_f, assay = "RNA", reductions = "umap")
cds <- cluster_cells(cds, reduction_method = "UMAP")
cds <- learn_graph(cds, use_partition = TRUE, close_loop = FALSE)

cell_ids <- colnames(cds)[colData(cds)$L2_CELL_TYPE == "Homeo-Glu-SLC1A2"]
closest <- as.matrix(
  cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex[colnames(cds), ]
)
root <- igraph::V(principal_graph(cds)[["UMAP"]])$name[
  as.numeric(names(which.max(table(closest[cell_ids, 1]))))
]
cds <- order_cells(cds, reduction_method = "UMAP", root_pr_nodes = root)

## --- 指定两条谱系 ---
leaf_keep <- c("Y_592", "Y_177")
g <- principal_graph(cds)[["UMAP"]]
stopifnot(all(c(root, leaf_keep) %in% names(V(g))))

Y <- cds@principal_graph_aux[["UMAP"]]$dp_mst
if (nrow(Y) == 2 && ncol(Y) != 2) {
  if (is.null(colnames(Y))) colnames(Y) <- names(V(g))
  Y <- t(Y)
}
if (is.null(rownames(Y))) rownames(Y) <- names(V(g))
Y <- Y[names(V(g)), 1:2, drop = FALSE]

umap <- as.data.frame(reducedDims(cds)[["UMAP"]])
colnames(umap) <- c("UMAP1", "UMAP2")
umap$pseudotime <- as.numeric(pseudotime(cds)[colnames(cds)])

cell_node <- paste0("Y_", closest[, 1])
path_nodes <- function(to) names(shortest_paths(g, from = root, to = to, output = "vpath")$vpath[[1]])
nodes_592 <- path_nodes("Y_592")
nodes_177 <- path_nodes("Y_177")

umap$lineage <- "other"
umap$lineage[cell_node %in% nodes_592 & !cell_node %in% nodes_177] <- "to_Y_592"
umap$lineage[!cell_node %in% nodes_592 & cell_node %in% nodes_177] <- "to_Y_177"
umap$lineage[cell_node %in% nodes_592 & cell_node %in% nodes_177]  <- "shared_trunk"

smooth_path <- function(nodes, n = 300, spar = 0.55) {
  xy <- Y[nodes, 1:2, drop = FALSE]
  if (nrow(xy) < 4) return(data.frame(UMAP1 = xy[, 1], UMAP2 = xy[, 2]))
  d <- c(0, cumsum(sqrt(rowSums(diff(xy)^2))))
  tt <- seq(min(d), max(d), length.out = n)
  data.frame(
    UMAP1 = predict(smooth.spline(d, xy[, 1], spar = spar), tt)$y,
    UMAP2 = predict(smooth.spline(d, xy[, 2], spar = spar), tt)$y
  )
}
curve_df <- rbind(
  transform(smooth_path(nodes_592), lineage = "to_Y_592"),
  transform(smooth_path(nodes_177), lineage = "to_Y_177")
)

ggplot(umap, aes(UMAP1, UMAP2)) +
  geom_point(aes(color = pseudotime), size = 0.25) +
  scale_color_viridis_c(option = "plasma", na.value = "grey80") +
  geom_path(
    data = curve_df,
    aes(UMAP1, UMAP2, group = lineage),
    color = "black", linewidth = 1.1, inherit.aes = FALSE
  ) +
  theme_classic() +
  coord_fixed()