library(Seurat)
library(ggplot2)
library(schard)
library(tidyr)
library(dplyr)
library(paletteer)
library(scCustomize)
# source("code/src/r/run_midas_clustering.R")
env <- new.env()
sys.source("code/src/r/annotation_utils.R", envir = env, keep.source = FALSE)

cols <- paletteer::paletteer_d(
  "pals::glasbey",
  n = 21,
  type = "discrete"
)
metadata = readxl::read_excel("data/metadata/2026-08-10_brain_atlas.xlsx")

obj <- readRDS("data/processed/Atlas_20260721/output/wnn/pure_imm/candidate_pure_imm_wnn_seurat.rds")

age_breaks <- c(2, 12, 20, 40, 60, Inf)
age_labels <- c("Child", "Adolescent", "YoungAdult", "MiddleAdult", "OlderAdult")

obj$AGE_GROUP <- cut(
  as.numeric(obj$AGE),
  breaks = age_breaks,
  right = FALSE,
  labels = age_labels
)

obj$GROUP <- factor(
  ifelse(
    obj$DISEASE == "Healthy",
    as.character(obj$AGE_GROUP),
    as.character(obj$DISEASE)
  ),
  levels = c(
    age_labels,
    "AD", "TDP", "FTLD", "PD", "Huntington", "ALS", "MS",
    "COVID-19", "Alcohol", "Drug", "MDD", "PTSD", "ASD"
  )
)

# obj@meta.data <- obj@meta.data |>
#   tibble::rownames_to_column("cell_id") |>
#   dplyr::left_join(metadata, by = c("batch" = "DONOR_ID")) |>
#   tibble::column_to_rownames("cell_id")

DefaultAssay(obj) <- "RNA"
Idents(obj) <- "wsnn_res.0.5"

obj <- env$apply_annotation(obj, 
  "parameters/annotation/microglia_wsnn_res0.5_v2.csv", 
  cluster_col = "wsnn_res.0.5", prefix = "mg_",overwrite = TRUE)
DotPlot(obj,
  features = env$annotation_markers("parameters/annotation/microglia_wsnn_res0.5_v2.csv"),
  group.by = "mg_level3"
) + RotatedAxis()+
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5))

DimPlot(obj,
  group.by = 'wsnn_res.0.5',
  reduction = 'wnn.umap',
  label = T)

plot_umap <- function(x, group.by, label = TRUE, cols = NULL, legend = TRUE) {
  p <- Seurat::DimPlot(
    x,
    group.by = group.by,
    reduction = "wnn.umap",
    cols = cols,
    label = label
  )
  if (!legend) p <- p + Seurat::NoLegend()
  p
}

plot_umap(obj, "wsnn_res.0.5", cols = cols)
plot_umap(obj, "RL3_1", cols = cols, label = FALSE)
plot_umap(obj, "DISEASE", cols = cols)
plot_umap(obj, "TITLE", label = FALSE, legend = FALSE)

L1_marker <- list(
  Neurons_Pan = c("SNAP25", "SYT1", "RBFOX3", "MAP2"),
  Excitatory = c("SLC17A7", "SLC17A6", "CAMK2A", "CUX2", "RORB"),
  Inhibitory = c("SLC32A1", "GAD1", "GAD2", "PVALB", "SST", "VIP", "LAMP5"),
  Dopaminergic = c("TH", "SLC6A3", "SLC18A2", "NR4A2", "LMX1A"),
  Astrocytes = c("SLC1A2", "AQP4", "GFAP", "ALDH1L1", "SOX9", "GLUL"),
  OPCs = c("PDGFRA", "CSPG4", "VCAN"),
  Oligodendrocytes = c("PLP1", "ST18", "MOBP", "MOG", "MBP"),
  Bergmann = c("ADAMTSL1", "GDF10", "TNC"),
  Microglia = c("C1QA", "C1QB", "CSF1R", "HEXB", "P2RY12", "TMEM119", "CX3CR1"),
  Macrophage_CAMs = c("MRC1", "CD163", "MS4A7","LYVE1","F13A1"),
  Immune_T_NK = c("PTPRC", "CD3D", "CD3E", "CD8A", "CD4", "SKAP1"),
  Endothelial = c("CLDN5", "FLT1", "PECAM1", "ABCB1"),
  Pericytes = c("PDGFRB", "VTN", "HIGD1B"),
  SMC = c("ACTA2", "TAGLN", "MYH11"),
  VLMC_Fibroblasts = c("DCN", "LUM", "COL1A1", "COL1A2"),
  Ependymal = c("FOXJ1", "PIFO"),
  Choroid_plexus = c("TTR", "FOLR1", "PRLR")
)

Inflammation_marker <- list(
  mg_cytokines = c(
    "IL1A", "IL1B", "IL1RN", "TNF", "IL6", "IL10", "IL15", "IL18",
    "TGFB1", "CSF1", "OSM"
  ),
  mg_chemokines = c(
    "CCL2", "CCL3", "CCL4", "CCL5", "CCL8", "CXCL1", "CXCL2",
    "CXCL8", "CXCL9", "CXCL10", "CXCL16"
  ),
  mg_inflam_proxy = c(
    "NFKBIA", "TNFAIP3", "PTGS2", "NLRP3", "CASP1", "C1QA", "C1QB", "C1QC"
  ),
  bam_identity = c(
    "MRC1", "LYVE1", "CD163", "F13A1", "MS4A7", "STAB1",
    "FOLR2", "PF4", "COLEC12", "SIGLEC1"
  ),
  bam_chemokines = c(
    "CXCL13", "CXCL12"
  ),
  bam_mhcII_inflam = c(
    "HLA-DRA", "HLA-DRB1", "HLA-DQA1", "CD74", "CIITA"
  ),
  t_cytokines = c(
    "IFNG", "CSF2", "IL2", "IL17A", "IL17F", "IL21", "IL22", "LTA"
  ),
  t_chemokines = c(
    "XCL1"
  ),
  t_trm_context = c(
    "CD69", "CXCR6", "ITGAE", "ZNF683", "PRDM1", "PDCD1"
  )
)

L1_marker <- lapply(
  L1_marker,
  \(genes) intersect(genes, rownames(obj))
)

DotPlot(
  obj,
  features = c(L1_marker,"PCDH9"),
  group.by = "wsnn_res.0.5",
  assay = "RNA",
  cluster.idents = TRUE
) +
  RotatedAxis() +
  labs(x = NULL, y = NULL)

DimPlot(
  obj, 
  cols = cols,
  split.by = "RL3_1",
  group.by = 'RL3_1',
  label = F,
  reduction = "wnn.umap")

obj@meta.data |>
  ggplot(aes(x = wsnn_res.0.5, fill = DISEASE)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = cols) +
  labs(x = "wsnn_res.0.5", y = "Proportion", fill = "Disease") +
  theme_minimal()

obj@meta.data |>
  ggplot(aes(x = DISEASE, fill = mg_level3)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = cols) +
  theme_minimal() +
  theme(
    panel.background = ggplot2::element_rect(fill = "white", colour = NA),
    plot.background = ggplot2::element_rect(fill = "white", colour = NA),
    panel.grid = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_text(size = 14, colour = "black", angle = 90, hjust = 1, vjust = 0.5),
    axis.title.y = ggplot2::element_blank()
  )

markers_wsnn_res.0.5 <- FindAllMarkers(
  obj,
  only.pos = TRUE,
  min.pct = 0.2,
  logfc.threshold = 0.25,
  test.use = "wilcox",
  max.cells.per.ident = 3000,
  random.seed = 3187
)

markers_wsnn_res.0.5 <- scCustomize::Add_Pct_Diff(markers_wsnn_res.0.5)

markers_filtered_ranked <- markers_wsnn_res.0.5 |>
  filter(
    p_val_adj < 0.05,
    !grepl("^(MT-|MT\\.|RPL|RPS|RNA)", gene, ignore.case = FALSE)
  ) |>
  arrange(cluster, desc(avg_log2FC), desc(pct_diff))


FeaturePlot(
  seu_f,
  features = c("PCDH9"),
  reduction = "wnn.umap",
  cols = c("gray95", "#911d1c"),
  pt.size = 0.5,
  min.cutoff = "q5",  # 去除底部10%的值
  max.cutoff = "q95"   # 去除顶部10%的值
)+patchwork::plot_annotation("All Regions")


