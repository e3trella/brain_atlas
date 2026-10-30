import scanpy as sc
import sys
sys.path.append("./code")
from src.py import plot_utils

adata = sc.read_h5ad("data/processed/Atlas_20260721/output/rna_midas_epoch2000.h5ad")

# batch effect evaluation

plot_utils.plot_highlight_batches_umap(
    adata,
    batches=[
        'C4812Cd',
        'H5608A',
        'H5608Cin',
        'C5227Cin',
        'H5263Cd',
        'C5700Cin',
        'C5227A',
        'H4161A'
    ],
)


# Normalize the data
adata.layers['counts'] = adata.X.copy()  # 保留原始计数数据
sc.pp.normalize_total(adata, target_sum=1e4)  # 每个细胞总计数归一化到 10,000
sc.pp.log1p(adata)  # 对数据进行 log 转换
adata.layers['lognorm'] = adata.X.copy()

# annotation ----------
sc.pl.umap(
    adata,
    color="midas_leiden",
    legend_loc="on data",
    legend_fontsize=8,
    frameon=False,
)

L1_marker = {
    # ---------------- 泛神经元与兴奋/抑制性大类 ----------------
    'Neurons_Pan': ["SNAP25", "SYT1", "RBFOX3", "MAP2"],
    
    # 谷氨酸能 / 兴奋性神经元 (皮层为主)
    'Excitatory': ["SLC17A7", "SLC17A6", "CAMK2A", "CUX2", "RORB"],
    
    # GABA能 / 抑制性神经元 (通用)
    'Inhibitory': ["SLC32A1", "GAD1", "GAD2", "PVALB", "SST", "VIP", "LAMP5"],
    
    # 中脑特异：多巴胺能神经元 (Dopaminergic Neurons)
    'Dopaminergic': ["TH", "SLC6A3", "SLC18A2", "NR4A2", "LMX1A"],
    
    # ---------------- 神经胶质细胞类 ----------------
    'Astrocytes': ["SLC1A2", "AQP4", "ALDH1L1", "SOX9", "GLUL"],  # snRNA 表现优异的组合
    'OPCs': ["PDGFRA", "CSPG4", "VCAN"],
    'Oligodendrocytes': ["PLP1", "ST18", "MOBP", "MOG", "MBP"],
    'Bergmann': ["ADAMTSL1","GDF10", "TNC"],  # 若数据仅含皮层/中脑，通常无 Bergmann 胶质(主要在小脑)

    # ---------------- 髓系与免疫细胞 ----------------
    'Microglia': ["C1QA", "C1QB", "CSF1R", "HEXB", "LRRC25"],
    'Macrophage_CAMs': ["MRC1", "CD163", "MS4A7"],  # 血管周围/硬膜巨噬细胞
    'Immune_T_NK': ["PTPRC", "CD3D", "CD3E","CD8A", "SKAP1"],
    
    # ---------------- 血管与基质细胞类 ----------------
    'Endothelial': ["CLDN5", "FLT1", "PECAM1", "ABCB1"],
    'Pericytes': ["PDGFRB", "VTN", "HIGD1B"],
    'SMC': ["ACTA2", "TAGLN", "MYH11"],  # 血管平滑肌细胞
    'VLMC_Fibroblasts': ["DCN", "LUM", "COL1A1", "COL1A2"],  # 血管周边成纤维样细胞
    
    # ---------------- 上皮/脉络丛类 ----------------
    'Ependymal': ["FOXJ1", "PIFO"],
    'Choroid_plexus': ["TTR", "FOLR1", "PRLR"]
}

L1_marker, info = plot_utils.filter_existing_genes(adata, L1_marker)

sc.tl.dendrogram(adata, groupby='midas_leiden', use_rep='X_midas')
sc.pl.dotplot(
    adata,
    var_names=L1_marker,
    groupby='midas_leiden',
    standard_scale='var',
    dendrogram=True,
    show=True
)

plot_utils.plot_highlight_batches_umap(
    adata,
    batch_col = 'midas_leiden',
    batches=[
        "10"
    ],
)

cl_annotation = dict.fromkeys((str(i) for i in range(57)), "Unassigned")

## adata2000
n_inhi_c = ['26','27','51','28','60','31','56','29','53','47','48','52']
n_exci_c = ['6','35','61','43','9','30','11','7','13','45','14','50','5','10','12','46']
n_Dopa_c = ["36"]
n_panOnly_c = ['41','58']

astro_c = ['37', '3', '4', '39', '38','22']
oligo_c = ['20','16','15','19','21','23','40','17','25','18','24']
astr_oligo_micro = ['2']
opc_c   = ['32','33','34','57','55']
micro_c = ['0','8','49']
T_c     = ['1']
chor_c  = ['59']

endo_c  = ['42']
fibro_c = ['44']
astro_epend_c = ['54']
astro_n_c = ["22"]

cl_annotation.update({k: "n_inhi_c" for k in n_inhi_c})
cl_annotation.update({k: "n_exci_c" for k in n_exci_c})
cl_annotation.update({k: "n_Dopa_c" for k in n_Dopa_c})
cl_annotation.update({k: "n_panOnly_c" for k in n_panOnly_c})
cl_annotation.update({k: "astro_c" for k in astro_c})
cl_annotation.update({k: "micro_c" for k in micro_c})
cl_annotation.update({k: "oligo_c" for k in oligo_c})
cl_annotation.update({k: "astr_oligo_micro" for k in astr_oligo_micro})
cl_annotation.update({k: "opc_c" for k in opc_c})
cl_annotation.update({k: "T_c" for k in T_c})
cl_annotation.update({k: "chor_c" for k in chor_c})
cl_annotation.update({k: "endo_c" for k in endo_c})
cl_annotation.update({k: "fibro_c" for k in fibro_c})
cl_annotation.update({k: "astro_epend_c" for k in astro_epend_c})
cl_annotation.update({k: "astro_n_c" for k in astro_n_c})


adata.obs["L1_CELL_TYPE"] = adata.obs['midas_leiden'].map(cl_annotation)
sc.pl.umap(
    adata,
    color=['L1_CELL_TYPE'],
    frameon=False,
    legend_loc='on data',
    legend_fontsize=8
)

mask = adata.obs["L1_CELL_TYPE"].isin(["micro_c", "astr_oligo_micro","T_c"])
adata_sub = adata[mask].copy()

plot_utils.plot_highlight_batches_umap(
    adata,
    batch_col = 'midas_leiden',
    batches=[
        "10","1"
    ],
)
