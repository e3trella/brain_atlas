library(Seurat)
library(ggplot2)
library(schard)
library(tidyr)
library(dplyr)
library(paletteer)
library(scCustomize)
source("code/src/r/run_midas_clustering.R")

obj = schard::h5ad2seurat(
  "data/processed/Atlas_20260721/output/astro_reMIDAS.h5ad",
  layers = c(counts = "counts", data = "lognorm")
)

metadata = readxl::read_excel("data/metadata/2026-08-10_brain_atlas.xlsx")

obj@meta.data <- obj@meta.data |>
  tibble::rownames_to_column("cell_id") |>
  dplyr::left_join(metadata, by = c("batch" = "DONOR_ID")) |>
  tibble::column_to_rownames("cell_id")

# KNN 构建邻接矩阵
obj <- FindNeighbors(obj, reduction = "Xmidas_", dims = 1:32)

# 聚类
obj <- FindClusters(obj, resolution = 0.5)

# 从 z_c 计算 UMAP
obj <- RunUMAP(
  object      = obj,
  reduction   = "Xmidas_",
  dims        = 1:32,
  n.neighbors = 30,
  min.dist    = 0.3,
  seed.use    = 42
)

cols <- paletteer::paletteer_d(
  "pals::glasbey",
  n = 23,
  type = "discrete"
)

obj <- readr::read_rds("data/processed/Atlas_20260721/output/astro_reMIDAS_subFrontal.rds")

DimPlot(seu_f,
  group.by = 'RNA_snn_res.0.5',
  reduction = 'umap',
  label = T)

L1_marker <- list(
  Neurons_Pan = c("SNAP25", "SYT1", "RBFOX3", "MAP2"),
  Excitatory = c("SLC17A7", "SLC17A6", "CAMK2A", "CUX2", "RORB"),
  Inhibitory = c("SLC32A1", "GAD1", "GAD2", "PVALB", "SST", "VIP", "LAMP5"),
  Dopaminergic = c("TH", "SLC6A3", "SLC18A2", "NR4A2", "LMX1A"),
  Astrocytes = c("SLC1A2", "AQP4", "GFAP", "ALDH1L1", "SOX9", "GLUL"),
  OPCs = c("PDGFRA", "CSPG4", "VCAN"),
  Oligodendrocytes = c("PLP1", "ST18", "MOBP", "MOG", "MBP"),
  Bergmann = c("ADAMTSL1", "GDF10", "TNC"),
  Microglia = c("C1QA", "C1QB", "CSF1R", "HEXB", "LRRC25"),
  Macrophage_CAMs = c("MRC1", "CD163", "MS4A7"),
  Immune_T_NK = c("PTPRC", "CD3D", "CD3E", "CD8A", "SKAP1"),
  Endothelial = c("CLDN5", "FLT1", "PECAM1", "ABCB1"),
  Pericytes = c("PDGFRB", "VTN", "HIGD1B"),
  SMC = c("ACTA2", "TAGLN", "MYH11"),
  VLMC_Fibroblasts = c("DCN", "LUM", "COL1A1", "COL1A2"),
  Ependymal = c("FOXJ1", "PIFO"),
  Choroid_plexus = c("TTR", "FOLR1", "PRLR")
)

L1_marker <- lapply(
  L1_marker,
  \(genes) intersect(genes, rownames(seu_f))
)

DotPlot(
  seu_f,
  features = L1_marker,
  group.by = "RNA_snn_res.0.5",
  assay = DefaultAssay(seu_f),
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
  reduction = "umap")

DimPlot(
  obj, 
  cols = cols,
  group.by = 'DISEASE',
  label = TRUE,
  reduction = "umap")

obj@meta.data |>
  ggplot(aes(x = RNA_snn_res.0.8, fill = DISEASE_NEW)) +
  geom_bar(position = "fill") +
  labs(x = "RNA_snn_res.0.8", y = "Proportion", fill = "Disease") +
  theme_minimal()

obj@meta.data |>
  ggplot(aes(x = DISEASE, fill = RL_3)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = cols) +
  labs(x = "DISEASE", y = "Proportion", fill = "RL_3") +
  theme_minimal()

DefaultAssay(obj) <- "RNA"
obj <- NormalizeData(obj, assay = "RNA")
Idents(obj) <- "RNA_snn_res.0.8"
markers_rna_snn_res_0_8 <- FindAllMarkers(
  obj,
  only.pos = TRUE,
  min.pct = 0.2,
  logfc.threshold = 0.25,
  test.use = "wilcox",
  max.cells.per.ident = 3000,
  random.seed = 3187
)

markers_filtered_ranked <- markers_rna_snn_res_0_8 |>
  filter(
    p_val_adj < 0.05,
    !grepl("^(MT-|MT\\.|RPL|RPS|RNA)", gene, ignore.case = FALSE)
  ) |>
  mutate(pct_diff = pct.1 - pct.2) |>
  arrange(cluster, desc(avg_log2FC), desc(pct_diff))

metadata <- readxl::read_xlsx("data/metadata/astro_mapper.xlsx")
obj$final_dec <- metadata$final_dec[match(as.character(Idents(obj)), as.character(metadata$cluster))]
obj$morphological <- metadata$morphological[match(as.character(Idents(obj)), as.character(metadata$cluster))]
obj$Function <- metadata$Function[match(as.character(Idents(obj)), as.character(metadata$cluster))]
obj$is_contain_wm <- metadata$is_contain_wm[match(as.character(Idents(obj)), as.character(metadata$cluster))]

# Add module score ----
gene_sets <- list(
  homeo = c(
    'SLC1A2', 'SLC1A3', 'GLUL', 'SLC6A11', 'KCNJ10', 'AQP4', 
    'SLC4A4', 'ATP1A2', 'SPARCL1', 'MEGF10', 'GPC6', 'ALDH1L1', 
    'FGFR3', 'GJB6', 'HEPACAM', 'ALDOC', 'NRXN1',
    'S100B', 'SOX9', 'GJA1', 'PRODH', 'PHGDH'
  ),

  reactive = c(
    'GFAP', 'TNC', 'CD44', 'CHI3L1', 'VCAN', 'HSPB1', 'CRYAB', 
    'HSPA1A', 'FKBP5', 'IL6ST', 'ID1', 'ID3', 'FOS', 'JUN', 
    'TIMP3', 'SERPINI2',
    'LCN2','CXCL10','C3','STAT3','VIM','NES'
  ),

  Neurotransmitter_Homeostasis = c(
    "SLC1A2", "SLC1A3", "GLUL", "SLC38A3", "SLC38A5",
    "SLC6A11", "SLC6A1", "GRM3", "GRM5", "GLUD1"
  ),
  Ion_Water_Homeostasis = c(
    "KCNJ10", "KCNJ16", "AQP4", "GJA1", "GJB6",
    "SLC12A2", "ATP1A2", "ATP1B2", "SLC4A4"
  ),
  Synaptogenesis = c(
    "THBS1", "THBS2", "THBS4", "SPARCL1", "SPARC",
    "GPC4", "GPC6", "CHRDL1", "TGFB1"
  ),
  Synaptic_Elimination = c(
    "MEGF10", "MERTK", "MFGE8", "C3", "C1S", "C4A"
  ),
  Metabolic_Support = c(
    "SLC16A1", "SLC16A4", "LDHA", "LDHB", "SLC2A1",
    "GYS1", "PYGB", "PCK2", "HK1", "PFKFB3"
  ),
  BBB_Neurovascular = c(
    "AQP4", "GJA1", "LAMA2", "DAG1", "AGRN",
    "PTGS2", "CYP2C8", "SLC7A10", "NOS3", "VEGFA"
  ),
  Immune_Inflammatory = c(
    "C3", "C1R", "C1S", "C4B", "SERPING1", "IL6", "IL1B",
    "TNF", "CXCL10", "CCL2", "TLR4", "NFKBIA", "STAT3", "SOCS3", "LCN2"
  ),
  Gliotransmission_Calcium = c(
    "GJA1", "GJB6", "PANX1", "SRR", "SLC7A11",
    "ITPR2", "S100B", "ENTPD1", "NT5E", "P2RY1", "VAMP2"
  ),
  Antioxidant_Neuroprotection = c(
    "NFE2L2", "GCLC", "GCLM", "GSS", "GPX1", "SOD1", "SOD2",
    "HMOX1", "NQO1", "MT1E", "MT2A", "SLC7A11", "BDNF", "GDNF", "CNTF"
  ),
  ECM_Glial_Scar = c(
    "GFAP", "VIM", "NES", "CSPG4", "ACAN", "VCAN",
    "TNC", "NCAN", "SERPINA3", "TIMP1", "MMP2"
  ),
  Lipid_Cholesterol = c(
    "APOE", "CLU", "ABCA1", "HMGCR", "HMGCS1", "FDFT1",
    "SQLE", "DHCR24", "CYP46A1", "LDLR", "ABCG1", "PLIN2"
  )
)

# 打印每个 gene set 中不在数据里的基因
lapply(gene_sets, setdiff, y = rownames(seu_f))
# 筛选不在 seuobj 中的基因
gene_sets <- lapply(gene_sets, intersect, y = rownames(seu_f))

## addModuleScore 计算得分
seu_f <- AddModuleScore(
  seu_f,
  features = gene_sets,
  name = names(gene_sets)
)

FeaturePlot(
  seu_f,
  features = names(gene_sets),
  reduction = "umap",
  cols = c("gray95", "#911d1c"),
  pt.size = 0.5,
  min.cutoff = "q5",  # 去除底部10%的值
  max.cutoff = "q95"   # 去除顶部10%的值
)+patchwork::plot_annotation("All Regions")


# Cerebrum 部分 ----
seu_cere <- subset(
  obj,
  subset = RL_3 == "Cerebral cortex" &
    !is.na(Function)
)

seu_cere <- run_midas_clustering(
  seu_cere,
  reduction = "Xmidas_",
  dims = 1:32,
  resolution = 0.8
)

DefaultAssay(seu_cere)
Idents(seu_cere) <- seu_cere$DISEASE_NEW
disease_levels <- setdiff(unique(seu_cere$DISEASE_NEW), "Healthy")

deg_by_disease <- lapply(
  disease_levels,
  function(disease) {
    FindMarkers(
      seu_cere,
      ident.1 = disease,
      ident.2 = "Healthy",
      test.use = "wilcox",
      logfc.threshold = 0,
      min.pct = 0.1,
      max.cells.per.ident = 5000,
      only.pos = FALSE
    ) |>
      tibble::rownames_to_column("gene") |>
      mutate(
        disease = disease,
        comparison = paste0(disease, "_vs_Healthy"),
        direction = case_when(
          avg_log2FC > 0 ~ paste0("higher_in_", disease),
          avg_log2FC < 0 ~ "higher_in_Healthy",
          TRUE ~ "no_logFC_difference"
        ),
        .before = 1
      ) |>
      arrange(p_val_adj, desc(abs(avg_log2FC)))
  }
)

names(deg_by_disease) <- disease_levels

deg_disease_vs_healthy <- bind_rows(deg_by_disease)

deg_disease_vs_healthy |>
  count(comparison, name = "n_tested")






annot <- c(
  "0"  = "Homeo-core-FGFR3",
  "10" = "Homeo-core-FGFR3",
  "1"  = "Homeo-Glu-SLC1A2",
  "6"  = "Homeo-Glu-SLC1A2",
  "2"  = "Homeo-ADH-SOX5",
  "3"  = "Reactive-ECM-DPP10",
  "4"  = "Reactive-ECM-DCLK1",
  "9"  = "Reactive-ECM-DCLK1",
  "5"  = "Homeo-GF-HGF",
  "7"  = "Reactive-OSM-CHI3L1",
  "11" = "Reactive-OSM-CHI3L1",
  "8"  = "Oligo-myelin-MBP",
  "12" = "Homeo-MT-RANBP3L",
  "13" = "Reactive-inflam-CHI3L1",
  "14" = "Homeo-Wnt-WIF1",
  "15" = "Reactive-IFN-GBP2",
  "16" = "Microglia-DAM-SPP1",
  "17" = "Homeo-Urea-SLC14A1",
  "18" = "Oligo-preOL-ST18",
  "19" = "Reactive-IEG-NR4A3",
  "20" = "Reactive-OxS-HMOX1",
  "21" = "Reactive-HSP-SLC7A11",
  "22" = "Neuron-ExN-SLC17A7",
  "23" = "Oligo-myelin-CLDN11"
)

stopifnot(all(seu_f$RNA_snn_res.0.5 %in% names(annot)))

# Factor indexing uses level position, not labels — convert to character first
seu_f$L2_CELL_TYPE <- unname(annot[as.character(seu_f$RNA_snn_res.0.5)])

celltype_levels <- c(
  "Homeo-core-FGFR3",
  "Homeo-Glu-SLC1A2",
  "Homeo-ADH-SOX5",
  "Homeo-GF-HGF",
  "Homeo-MT-RANBP3L",
  "Homeo-Wnt-WIF1",
  "Homeo-Urea-SLC14A1",
  "Reactive-ECM-DPP10",
  "Reactive-ECM-DCLK1",
  "Reactive-OSM-CHI3L1",
  "Reactive-inflam-CHI3L1",
  "Reactive-IFN-GBP2",
  "Reactive-OxS-HMOX1",
  "Reactive-IEG-NR4A3",
  "Reactive-HSP-SLC7A11",
  "Oligo-myelin-MBP",
  "Oligo-preOL-ST18",
  "Oligo-myelin-CLDN11",
  "Microglia-DAM-SPP1",
  "Neuron-ExN-SLC17A7"
)

seu_f$L2_CELL_TYPE <- factor(seu_f$L2_CELL_TYPE, levels = celltype_levels)

# 检查
table(seu_f$RNA_snn_res.0.5, seu_f$L2_CELL_TYPE)
DimPlot(
  seu_f, 
  group.by = "L2_CELL_TYPE", 
  label = FALSE, 
  repel = TRUE,
  reduction = "umap", 
  cols = paletteer::paletteer_d("pals::glasbey", n = length(celltype_levels), type = "discrete")
)

# monocle3 ----
# 让 AI 参考 @code/src/r/monocle_demo.R 来跑即可
