import scanpy as sc
import sys
sys.path.append("./code")
from src.py import plot_utils
from src.py.midas_clustering import run_midas_clustering

adata = sc.read_h5ad("data/processed/Atlas_20260721/output/astro_reMIDAS.h5ad")
sc.pp.normalize_total(adata, target_sum=1e4)  # 每个细胞总计数归一化到 10,000
sc.pp.log1p(adata)  # 对数据进行 log 转换

# astro peeling iteratively ----
# workflow
run_midas_clustering(
    adata,
    cluster_key="midas_leiden_R1.0",
    n_neighbors=30,
    use_rep="X_midas",
    neighbors_key="midas_neighbors",
    resolution=1.0,
    random_state=42,
    plot=True,
    legend_loc="on data",
    legend_fontsize=8,
    legend_fontoutline=0,
    frameon=False,
)

L1_marker, info = plot_utils.filter_existing_genes(adata, L1_marker)

sc.tl.dendrogram(adata, groupby='midas_leiden_R1.0', use_rep='X_midas')
sc.pl.dotplot(
    adata,
    var_names=L1_marker,
    groupby='midas_leiden_R1.0',
    standard_scale='var',
    show=True
)

plot_utils.plot_highlight_batches_umap(
    adata,
    batch_col = 'batch',
    batches=[
        "H4161A",
        "C5227A",
        "H5109Cd",
        "Healthy_MM_437"
    ],
)

sc.pl.umap(
    adata,
    color=['midas_leiden_R1.0'],
    frameon=False,
    legend_loc='on data',
    legend_fontsize=8,
    legend_fontoutline=2,
    legend_fontweight='bold'
)

vc = adata.obs.loc[adata.obs['midas_leiden_R1.0'].isin(['12']), 'batch'].value_counts().reset_index()
print(f"\n非零 batch 数: {(vc > 0).sum()}  总细胞数: {vc.sum()}")

print(
    adata.obs.groupby('midas_leiden')['batch']
    .value_counts()
    .loc[lambda s: s > 10]
    .groupby(level=0)
    .size()
    .sort_values(ascending=False)
)
