import os
import scanpy as sc
import anndata as ad
from pathlib import Path

def ensure_dir_exists(file_path):
    """
    确保文件所在的目录存在，不存在则创建
    """
    if isinstance(file_path, str):
        file_path = Path(file_path)
    
    dir_path = file_path.parent
    
    if not dir_path.exists():
        dir_path.mkdir(parents=True, exist_ok=True)
        print(f"✅ 已创建目录: {dir_path}")
    else:
        print(f"📁 目录已存在: {dir_path}")

# 生成subset_ list
def generate_subset_path_list(
    base_dir, 
    subset_count, 
    subset_prefix='subset_',
    file_prefix='cell_names',
    file_suffix='.csv'
):
    """
    生成 subset_ list 的路径列表
    """
    subset_paths = []
    for i in range(subset_count):
        subset_name = f"{subset_prefix}{i}"
        file_name = f"{file_prefix}{file_suffix}"
        subset_path = Path(base_dir) / subset_name / file_name
        subset_paths.append(subset_path)
    
    return subset_paths

def export_mudata_rna_to_r(
    mdata,
    output_path: str,
    counts_layer: str = "counts",
    include_layers: list = None,
    include_obsm: list = None  # 新增参数
):
    """
    Export RNA modality from MuData to h5ad for R analysis.
    """
    # Extract RNA modality
    rna = mdata.mod["rna"].copy()
    
    # Use MuData .obs as metadata
    rna.obs = mdata.obs.copy()
    
    # 复制 obsm 降维数据
    if include_obsm is not None:
        for key in include_obsm:
            if key in mdata.obsm:
                rna.obsm[key] = mdata.obsm[key]
    else:
        # 默认导出所有 obsm
        rna.obsm = mdata.obsm.copy()
    
    # 复制 uns 中的相关数据（如 neighbors, umap 参数）
    for key in ['neighbors', 'umap']:
        if key in mdata.uns:
            rna.uns[key] = mdata.uns[key]
    
    # Handle layers
    if counts_layer in rna.layers:
        if include_layers is not None:
            layers_to_keep = {counts_layer: rna.layers[counts_layer]}
            for layer in include_layers:
                if layer in rna.layers and layer != counts_layer:
                    layers_to_keep[layer] = rna.layers[layer]
            rna.layers = layers_to_keep
        else:
            rna.layers = {counts_layer: rna.layers[counts_layer]}
    else:
        raise ValueError(f"Counts layer '{counts_layer}' not found in RNA modality")
    
    rna.write_h5ad(output_path)
    
    print(f"Saved RNA modality to: {output_path}")
    print(f"  - Observations (cells): {rna.n_obs}")
    print(f"  - Variables (genes): {rna.n_vars}")
    print(f"  - Layers: {list(rna.layers.keys())}")
    print(f"  - obsm: {list(rna.obsm.keys())}")
    
    return output_path