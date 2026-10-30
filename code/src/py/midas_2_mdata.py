import pandas as pd
import scanpy as sc
import numpy as np
import anndata as ad
import muon as mu
import scipy.sparse as sp
import mudata as md
import gc
import warnings
from typing import Union, Dict, Any
warnings.filterwarnings('ignore')


def create_mudata_from_midas(h5ad_path, feat_dir, cell_names_dir):
    """
    从 MIDAS 输出创建 MuData 对象
    
    Parameters:
        h5ad_path: midas_out 中的 h5ad 文件路径
        feat_dir: feature 文件所在目录
        cell_names_dir: cell_names.csv 所在目录
    
    Returns:
        MuData 对象,包含可用的模态(RNA/ATAC/ADT)和 joint latent
    """
    
    # 1. 读取 h5ad 文件
    midas_adata = sc.read_h5ad(h5ad_path)
    
    # 定义模态配置
    modalities = {
        'rna': {'obsm_key': 'x_bc_rna', 'feat_file': 'feat_names_rna.csv'},
        'atac': {'obsm_key': 'x_bc_atac', 'feat_file': 'feat_names_atac.csv'},
        'adt': {'obsm_key': 'x_bc_adt', 'feat_file': 'feat_names_adt.csv'}
    }
    
    # 读取 cell names
    cell_names = pd.read_csv(f"{cell_names_dir}/cell_names.csv")['x'].values.tolist()
    
    # 检测可用的模态并创建 AnnData
    adatas = {}
    for mod_name, mod_config in modalities.items():
        if mod_config['obsm_key'] in midas_adata.obsm:
            feat = pd.read_csv(f"{feat_dir}/{mod_config['feat_file']}")['x'].values.tolist()
            adata = ad.AnnData(X=midas_adata.obsm[mod_config['obsm_key']])
            adata.obs_names = cell_names
            adata.var_names = feat
            adata.obs = midas_adata.obs.copy()
            adatas[mod_name] = adata
    
    # 合并成 MuData
    mdata = mu.MuData(adatas)
    
    # 传入 obs
    mdata.obs = midas_adata.obs.copy()
    
    # 将 z_c_joint 和 z_u_joint 放入 mudata.obsm
    if 'z_c_joint' in midas_adata.obsm:
        mdata.obsm['z_c_joint'] = midas_adata.obsm['z_c_joint']
    if 'z_u_joint' in midas_adata.obsm:
        mdata.obsm['z_u_joint'] = midas_adata.obsm['z_u_joint']
    
    return mdata


# 使用示例
# mdata = create_mudata_from_midas(
#     h5ad_path='midas_out/w1.h5ad',
#     feat_dir='dataset/teadog_mosaic_mtx/data/feat',
#     cell_names_dir='dataset/teadog_mosaic_mtx/data/w1'
# )

# mdata


def build_adata_mdata(outputs: Dict[str, Any], consume_input: bool = True) -> Union[ad.AnnData, md.MuData]:
    """
    将基于 Batch 分块的多模态预测字典，安全、极速地组装为 AnnData 或 MuData。
    
    动态特性:
        - 自动检测存在的模态 (基于第一个 batch 的 'x_bc' keys)。
        - 自动检测潜空间特征 (动态提取包含 'joint' 的键，不再写死 z_c, z_u, s)。
        - 单模态 -> 返回 AnnData，潜空间存放于 adata.obsm。
        - 多模态 -> 返回 MuData，潜空间存放于 mdata.obsm。
        
    参数:
        outputs: 模型的输出字典。
        consume_input: 是否在提取数据的同时销毁原输入 (赋值为 None)，以达到极致的内存节省。
    """
    # 兼容处理：检查输入是否直接是 pred 字典，或是包含 pred 的大字典
    pred = outputs.get('pred', outputs)
    if not pred:
        raise ValueError("输入字典缺少有效数据或数据为空。")

    batch_keys = list(pred.keys())
    first_batch = batch_keys[0]
    first_v = pred[first_batch]
    
    # ---------------------------------------------------------
    # 0. 动态嗅探 (Dynamic Sniffing)
    # ---------------------------------------------------------
    if 'x_bc' not in first_v:
        raise KeyError(f"Batch '{first_batch}' 中未找到 'x_bc' 数据。")
        
    # 动态探测组学模态
    mods = list(first_v['x_bc'].keys())
    
    # 动态探测潜空间变量 (提取所有包含 'joint' 键的顶层字典)
    latent_keys =[k for k in first_v.keys() if isinstance(first_v[k], dict) and 'joint' in first_v[k]]
    
    print(f"[*] 探测到 {len(mods)} 个组学模态: {mods}")
    print(f"[*] 探测到 {len(latent_keys)} 个潜空间特征: {latent_keys}")
    
    # 初始化数据容器
    all_latents = {k:[] for k in latent_keys}
    all_x_bc = {m:[] for m in mods}

    # ---------------------------------------------------------
    # 1. 动态提取、前置稀疏化与显式内存切断
    # ---------------------------------------------------------
    print(f"[*] 正在执行数据提取与稀疏化 (Memory-safe: {consume_input})...")
    for k in batch_keys:
        print(f"    - 处理 batch: {k} ...")
        v = pred[k]
        
        # 提取低维潜空间特征 (保持密集 NumPy 格式)
        for lat_key in latent_keys:
            all_latents[lat_key].append(v[lat_key]['joint'])
            if consume_input:
                v[lat_key]['joint'] = None  # 释放潜空间内存
                
        # 提取高维组学矩阵：动态遍历所有模态
        for m in mods:
            X = v['x_bc'][m]
            # 优化：判断是否已是稀疏矩阵，避免重复转换
            if not sp.issparse(X):
                X = sp.csr_matrix(X)
            all_x_bc[m].append(X)
            
            # 内存核心：立刻切断该 batch 中对应模态 Dense 矩阵的引用
            if consume_input:
                v['x_bc'][m] = None 
                
        if consume_input:
            gc.collect()

    print("[*] 数据提取与稀疏化完成。")

    # ---------------------------------------------------------
    # 2. 生成全局细胞 Metadata
    # ---------------------------------------------------------
    # 优先使用第一个模态的矩阵行数来计算每个 batch 的细胞数量
    cell_counts =[matrix.shape[0] for matrix in all_x_bc[mods[0]]]
    obs_batch = np.repeat(batch_keys, cell_counts)
    obs_df = pd.DataFrame({'batch': pd.Categorical(obs_batch)})
    # ---------------------------------------------------------
    # 3. 特征的极速合并 (Dense Arrays & Sparse Matrices)
    # ---------------------------------------------------------
    print("[*] 正在合并低维潜空间...")
    merged_latents = {}
    for lat_key, matrix_list in all_latents.items():
        if matrix_list:
            merged_latents[lat_key] = np.concatenate(matrix_list, axis=0)
    del all_latents

    print("[*] 正在合并高维稀疏矩阵...")
    merged_x_bc = {}
    for m in mods:
        merged_x_bc[m] = sp.vstack(all_x_bc[m], format='csr')
    del all_x_bc
    gc.collect()

    # ---------------------------------------------------------
    # 4. 动态构建 AnnData 或 MuData
    # ---------------------------------------------------------
    # 分支 A：单模态情况 -> 返回 AnnData
    if len(mods) == 1:
        print("[*] 仅检测到单模态，正在构建 AnnData 对象...")
        m = mods[0]
        adata = ad.AnnData(X=merged_x_bc[m], obs=obs_df.copy())
        
        # 挂载潜空间
        for lat_key, lat_matrix in merged_latents.items():
            adata.obsm[lat_key] = lat_matrix
            
        del merged_x_bc, merged_latents
        gc.collect()
        
        print(f"[*] 组装完成！单模态 AnnData 包含细胞数: {adata.n_obs}")
        return adata

    # 分支 B：多模态情况 -> 返回 MuData
    else:
        print("[*] 检测到多模态，正在构建 MuData 对象...")
        adatas = {}
        for m in mods:
            adatas[m] = ad.AnnData(X=merged_x_bc[m], obs=obs_df.copy())
            
        mdata = md.MuData(adatas)
        
        # 将潜空间挂载到全局 (Top-level) obsm
        for lat_key, lat_matrix in merged_latents.items():
            mdata.obsm[lat_key] = lat_matrix
            
        # 将批次信息提升到全局，确保跨模态可视化可用
        mdata.obs['batch'] = obs_df['batch'].values
        
        del merged_x_bc, merged_latents, adatas
        gc.collect()
        
        print(f"[*] 组装完成！多模态 MuData 包含细胞数: {mdata.n_obs}")
        return mdata

