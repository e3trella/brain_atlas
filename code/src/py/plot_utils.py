import matplotlib.pyplot as plt
import muon as mu
import numpy as np
import matplotlib.colors as mcolors
import scanpy as sc

def plot_umap_highlight(
    mdata,
    umap_key="X_umap_midas",
    celltype_col="celltype_r2",
    target_celltypes=None,
    highlight_colors=None,
    background_color="lightgrey",
    bg_size=1,
    fg_size=2,
    bg_alpha=0.3,
    fg_alpha=0.9,
    figsize=(6,6),
    show_legend=True,
    title=None,
    save_path=None
):
    """
    绘制 UMAP, 高亮一个或多个细胞类型, 其余细胞置灰。

    Parameters
    ----------
    mdata : MuData 或 AnnData 对象
        包含 UMAP 和细胞类型信息
    umap_key : str
        存放 UMAP 坐标的 obsm key
    celltype_col : str
        存放细胞类型的 obs 列名
    target_celltypes : list of str
        要高亮的细胞类型列表
    highlight_colors : list of str
        高亮细胞类型对应颜色
        如果 None, 将自动生成颜色
    background_color : str
        背景细胞颜色
    bg_size : float
        背景点大小
    fg_size : float
        高亮点大小
    bg_alpha : float
        背景点透明度
    fg_alpha : float
        高亮点透明度
    figsize : tuple
        图像尺寸
    show_legend : bool
        是否显示图例
    title : str
        图标题
    save_path : str
        如果指定, 将图保存到文件
    """
    
    X = mdata.obsm[umap_key]
    celltype = mdata.obs[celltype_col]
    
    if target_celltypes is None:
        raise ValueError("请至少指定一个 target_celltypes")
    if isinstance(target_celltypes, str):
        target_celltypes = [target_celltypes]

    # 自动生成颜色
    if highlight_colors is None:
        cmap = plt.get_cmap("tab10")
        highlight_colors = [cmap(i % 10) for i in range(len(target_celltypes))]
    if len(highlight_colors) != len(target_celltypes):
        raise ValueError("highlight_colors 长度必须和 target_celltypes 一致")

    plt.figure(figsize=figsize)

    # 背景：所有非目标细胞
    mask_bg = ~celltype.isin(target_celltypes)
    plt.scatter(
        X[mask_bg, 0],
        X[mask_bg, 1],
        s=bg_size,
        color=background_color,
        alpha=bg_alpha,
        linewidths=0
    )

    # 高亮每个目标细胞类型
    for ct, color in zip(target_celltypes, highlight_colors):
        mask = celltype == ct
        plt.scatter(
            X[mask, 0],
            X[mask, 1],
            s=fg_size,
            color=color,
            alpha=fg_alpha,
            linewidths=0,
            label=ct
        )

    if show_legend:
        plt.legend(frameon=False)
    plt.xlabel("UMAP1")
    plt.ylabel("UMAP2")
    if title is None:
        if len(target_celltypes) == 1:
            title = f"Highlighted cell type: {target_celltypes[0]}"
        else:
            title = f"Highlighted cell types: {', '.join(target_celltypes)}"
    plt.title(title)
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=300)
        print(f"Figure saved to {save_path}")
    
    plt.show()

def filter_existing_genes(adata, gene_list):
    """
    过滤掉不存在于 adata.var_names 中的基因
    
    支持两种输入格式：
    1. 简单列表: ["INA", "SLC17A6", ...]
    2. 字典格式: {"Neurons": ["INA"], "Excitatory": ["SLC17A6", "SLC17A7"], ...}
    
    Parameters:
    -----------
    adata : AnnData
        单细胞数据对象
    gene_list : list or dict
        基因列表或细胞类型到基因列表的字典
    
    Returns:
    --------
    filtered : list or dict
        与输入相同格式, 但只包含存在的基因
    info : dict
        包含被过滤掉的基因信息
    """
    
    # 判断输入类型
    is_dict = isinstance(gene_list, dict)
    
    if is_dict:
        # 处理字典格式
        filtered = {}
        all_missing = []
        total_original = 0
        total_existing = 0
        
        for cell_type, markers in gene_list.items():
            existing_markers = [g for g in markers if g in adata.var_names]
            missing_markers = [g for g in markers if g not in adata.var_names]
            
            # 只保留有基因的类别
            if len(existing_markers) > 0:
                filtered[cell_type] = existing_markers
            
            all_missing.extend(missing_markers)
            total_original += len(markers)
            total_existing += len(existing_markers)
        
        print(f"原始细胞类型数: {len(gene_list)}")
        print(f"保留细胞类型数: {len(filtered)}")
        print(f"原始基因总数: {total_original}")
        print(f"存在基因总数: {total_existing}")
        print(f"缺失基因: {list(set(all_missing))}")  # 去重显示
        
        info = {
            "missing": list(set(all_missing)),
            "original_count": total_original,
            "existing_count": total_existing,
            "removed_cell_types": set(gene_list.keys()) - set(filtered.keys())
        }
        
    else:
        # 处理简单列表格式
        filtered = [g for g in gene_list if g in adata.var_names]
        missing = [g for g in gene_list if g not in adata.var_names]
        
        print(f"原始基因数: {len(gene_list)}")
        print(f"存在基因数: {len(filtered)}")
        print(f"缺失基因: {missing}")
        
        info = {
            "missing": missing,
            "original_count": len(gene_list),
            "existing_count": len(filtered)
        }
    
    return filtered, info

# 使用示例
# L1_marker_filtered, info = filter_existing_genes(mdata['rna'], L1_marker)

def embedding_split(mdata, basis, split_key, size=5):
    import matplotlib.pyplot as plt

    classes = mdata.obs[split_key].unique()
    fig, axes = plt.subplots(1, len(classes), figsize=(4*len(classes), 4))

    for ax, cls in zip(axes, classes):
        subset = mdata[mdata.obs[split_key] == cls]
        mu.pl.embedding(
            subset,
            basis=basis,
            color=split_key,
            size=size,
            frameon=False,
            ax=ax,
            legend_loc=None
        )
        ax.set_title(str(cls))
    return fig



def plot_highlight_batches_umap(
    adata_obj,
    batches,
    batch_col="batch",
    basis="X_umap",
    size=8,
    other_color="lightgray",
    title=None,
):
    """Plot a UMAP highlighting selected batches with distinct colors.

    Cells whose ``batch_col`` value is included in ``batches`` are colored by
    batch. All other cells are assigned to an ``"Other"`` group and shown in
    ``other_color``. The function writes a temporary ``"highlight_batch"``
    column to ``adata_obj.obs`` for plotting.

    Parameters
    ----------
    adata_obj : anndata.AnnData
        AnnData object containing cell metadata in ``.obs`` and embedding
        coordinates in ``.obsm``.
    batches : iterable of str
        Batch values to highlight. Values are compared against ``batch_col``
        after converting the column to string.
    batch_col : str, default "batch"
        Column in ``adata_obj.obs`` that stores batch labels.
    basis : str, default "X_umap"
        Key in ``adata_obj.obsm`` containing the embedding coordinates to plot.
    size : int or float, default 8
        Point size passed to ``scanpy.pl.embedding``.
    other_color : str, default "lightgray"
        Color used for cells not belonging to the selected batches.
    title : str or None, default None
        Plot title. If ``None``, a default title based on ``batch_col`` is used.

    Raises
    ------
    ValueError
        If ``batch_col`` is not in ``adata_obj.obs``, ``basis`` is not in
        ``adata_obj.obsm``, or any requested batch is absent from ``batch_col``.
    """
    # Validate that the required metadata column and embedding are available.
    if batch_col not in adata_obj.obs.columns:
        raise ValueError(f"{batch_col!r} is not in adata_obj.obs.columns")

    if basis not in adata_obj.obsm:
        raise ValueError(f"{basis!r} is not in adata_obj.obsm")

    # Convert selected batches to a list so their order is stable for coloring.
    batches = list(batches)
    missing_batches = sorted(set(batches) - set(adata_obj.obs[batch_col].astype(str)))
    if missing_batches:
        raise ValueError(f"These batches are not present in {batch_col!r}: {missing_batches}")

    # Create a plotting column: selected batches keep their label; all other
    # cells are grouped as "Other".
    highlight_col = "highlight_batch"
    batch_values = adata_obj.obs[batch_col].astype(str)
    adata_obj.obs[highlight_col] = np.where(
        batch_values.isin(batches),
        batch_values,
        "Other",
    )
    adata_obj.obs[highlight_col] = adata_obj.obs[highlight_col].astype("category")

    # Build a palette with a neutral color for background cells and one distinct
    # color per selected batch.
    cmap = list(mcolors.TABLEAU_COLORS.values()) + list(mcolors.CSS4_COLORS.values())
    palette = {"Other": other_color}
    palette.update({batch: cmap[i % len(cmap)] for i, batch in enumerate(batches)})

    sc.pl.embedding(
        adata_obj,
        basis=basis,
        color=highlight_col,
        palette=palette,
        groups=batches,
        size=size,
        frameon=False,
        title=title or f"Highlighted batches in {batch_col}",
    )
