import scanpy as sc
import pandas as pd
from functools import reduce
from typing import Optional, Literal

def find_markers_by_region(
    mdata, 
    filters=None,
    group_col="DISEASE_NEW", 
    target_group="FTLD", 
    reference_group="Healthy",
    layer="lognorm",
    method="wilcoxon",
    pts=True,
    min_logfc=0.25,
    pval_adj_thresh=0.05
):
    """
    在指定条件内，比较 target_group vs reference_group 的差异基因，返回 DataFrame。
    
    Parameters
    ----------
    mdata : MuData
        包含 RNA 数据的 MuData 对象
    filters : dict or None
        筛选条件字典，例如 {"RL_2": "Cerebrum", "cell_type": "Neuron"}。
        value 为 "all" 时跳过该列的筛选；为 None 时不做任何筛选。
        也支持 list，如 {"RL_2": ["Cerebrum", "Brainstem"]} 表示多值筛选。
    group_col : str
        分组列名，例如 DISEASE_NEW
    target_group : str
        目标组
    reference_group : str
        对照组
    layer : str
        RNA 数据层，例如 lognorm
    method : str
        差异分析方法，例如 "wilcoxon"
    pts : bool
        是否计算 pct.1 和 pct.2
    min_logfc : float
        logFC 绝对值阈值，低于此值的基因将被过滤
    pval_adj_thresh : float
        FDR 阈值，高于此值的基因将被过滤
    
    Returns
    -------
    pandas.DataFrame
        包含 marker gene 的 DataFrame，类似 Seurat FindMarkers
    """
    # 复制 RNA 数据
    adata = mdata['rna'].copy()
    adata.obs = mdata.obs.copy()

    # ── 修改1&2：支持多列、多值筛选，value=="all" 时跳过 ──────────────────────
    if filters is not None:
        for col, val in filters.items():
            if val == "all":          # 明确传 "all" → 不筛选该列
                continue
            elif isinstance(val, list):
                adata = adata[adata.obs[col].isin(val)].copy()
            else:
                adata = adata[adata.obs[col] == val].copy()
    # filters=None 时完全不筛选

    # 只保留目标组和对照组
    if reference_group != "rest":
        adata = adata[adata.obs[group_col].isin([target_group, reference_group])].copy()
    else:
        pass
    # 清理 category
    adata.obs[group_col] = adata.obs[group_col].astype("category").cat.remove_unused_categories()

    # 使用指定数据层
    adata.X = adata.layers[layer]

    # 差异分析
    sc.tl.rank_genes_groups(
        adata,
        groupby=group_col,
        groups=[target_group],
        reference=reference_group,
        method=method,
        pts=pts
    )

    # 获取 DataFrame
    df = sc.get.rank_genes_groups_df(adata, group=target_group)

    # ── 修改3：按 logFC 绝对值过滤 ────────────────────────────────────────────
    df = df[df["logfoldchanges"].abs() >= min_logfc]

    # ── 修改4：按校正 p 值过滤 ────────────────────────────────────────────────
    df = df[df["pvals_adj"] <= pval_adj_thresh]

    return df.reset_index(drop=True)



def get_common_degs(
    file_paths: list[str],
    direction: Literal['up', 'down', 'both'] = 'up',
    lfc_threshold: float = 0,
    pval_adj_threshold: Optional[float] = None,
) -> set[str] | dict[str, set[str]]:
    """
    从多个DEG文件中找出共同差异表达基因的交集。

    Parameters
    ----------
    file_paths : list[str]
        各疾病DEG文件的完整路径列表
    direction : 'up' | 'down' | 'both'
        筛选方向：上调、下调、或同时返回两者
    lfc_threshold : float
        log fold change 的绝对值阈值，默认 0
    pval_adj_threshold : float, optional
        adjusted p-value 上限，如 0.05; 为 None 时不过滤

    Returns
    -------
    set[str]
        direction 为 'up' 或 'down' 时，返回共有基因集合
    dict[str, set[str]]
        direction 为 'both' 时，返回 {'up': set, 'down': set}
    """
    if not file_paths:
        raise ValueError("file_paths 不能为空")

    def _filter(df: pd.DataFrame, direc: Literal['up', 'down']) -> set[str]:
        if direc == 'up':
            mask = df['logfoldchanges'] > lfc_threshold
        else:
            mask = df['logfoldchanges'] < -lfc_threshold
        if pval_adj_threshold is not None:
            mask &= df['pvals_adj'] < pval_adj_threshold
        return set(df.loc[mask, 'names'].tolist())

    directions = ['up', 'down'] if direction == 'both' else [direction]
    results = {}

    for direc in directions:
        gene_sets = []
        for path in file_paths:
            df = pd.read_csv(path)
            genes = _filter(df, direc)
            print(f"[{direc}] {path}: {len(genes)} genes")
            gene_sets.append(genes)

        common = reduce(lambda x, y: x & y, gene_sets)
        print(f"\n共有{direc}调基因数量: {len(common)}\n")
        results[direc] = common

    return results if direction == 'both' else results[direction]
