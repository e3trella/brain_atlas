import pandas as pd
import scanpy as sc


def run_midas_clustering(
    adata,
    cluster_key,
    n_neighbors=30,
    use_rep="X_midas",
    neighbors_key="midas_neighbors",
    resolution=1.0,
    random_state=42,
    plot=True,
    legend_loc="on data",
    legend_fontsize=8,
    legend_fontoutline=2,
    frameon=False,
):
    """
    基于 Midas 表示完成 KNN -> Leiden -> UMAP 聚类流程，并绘制 UMAP。
    簇编号按簇内细胞数从大到小重新编号，0 对应细胞数最多的簇。

    Parameters
    ----------
    adata : AnnData
        输入数据，需包含 use_rep 指定的 obsm 表示
    cluster_key : str
        Leiden 聚类结果写入 adata.obs 的列名，同时用作 UMAP 着色列
    n_neighbors : int
        KNN 邻居数
    use_rep : str
        构建邻居图所用的 obsm key，默认 "X_midas"
    neighbors_key : str
        邻居图在 adata.uns 中的 key
    resolution : float
        Leiden 分辨率
    random_state : int
        随机种子
    plot : bool
        是否绘制 UMAP
    legend_loc : str
        传递给 sc.pl.umap 的图例位置
    legend_fontsize : int
        图例字号
    legend_fontoutline : int
        图例文字描边宽度
    frameon : bool
        是否绘制坐标轴边框

    Returns
    -------
    None
        adata 原地更新，写入邻居图、cluster_key 列及 UMAP 坐标
    """
    sc.pp.neighbors(
        adata,
        n_neighbors=n_neighbors,
        use_rep=use_rep,
        key_added=neighbors_key,
        random_state=random_state,
    )
    sc.tl.leiden(
        adata,
        neighbors_key=neighbors_key,
        key_added=cluster_key,
        resolution=resolution,
        random_state=random_state,
        flavor="igraph",
        n_iterations=2,
        directed=False,
    )
    # 按簇大小（细胞数）从大到小重新编号，0 对应最大的簇
    sizes = adata.obs[cluster_key].astype(str).value_counts()
    order = sizes.index.tolist()
    relabel = {old: str(new) for new, old in enumerate(order)}
    adata.obs[cluster_key] = pd.Categorical(
        adata.obs[cluster_key].astype(str).map(relabel),
        categories=[str(i) for i in range(len(order))],
    )
    sc.tl.umap(
        adata,
        neighbors_key=neighbors_key,
        random_state=random_state,
    )
    if plot:
        sc.pl.umap(
            adata,
            color=cluster_key,
            legend_loc=legend_loc,
            legend_fontsize=legend_fontsize,
            legend_fontoutline=legend_fontoutline,
            frameon=frameon,
        )
