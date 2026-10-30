"""
从 MIDAS 输出的 RNA 与 ATAC 两个 h5ad 中按同一组细胞取子集，并输出诊断信息。
输出的子集文件供 R 中 schard::h5ad2seurat 读取。

运行前请修改下方四个路径。cells_keep.txt 为每行一个 barcode，不含表头。
"""
import sys
import numpy as np
import pandas as pd
import anndata as ad
import scipy.sparse as sp

RNA_IN = "data/processed/Atlas_20260721/output/rna/rna_midas_epoch2000.h5ad"      # 请改为实际 RNA 文件路径
ATAC_IN = "data/processed/Atlas_20260721/output/atac/atac_midas_epoch2000.h5ad"

CELLS = "data/processed/Atlas_20260721/output/wnn/pure_astro/pure_astro_cells_keep.txt"
RNA_OUT = "data/processed/Atlas_20260721/output/wnn/pure_astro/pure_astro_rna.h5ad"
ATAC_OUT = "data/processed/Atlas_20260721/output/wnn/pure_astro/pure_astro_atac.h5ad"

SPARSIFY_FRAC = 0.5          # 稠密矩阵非零比例低于该值才转为 CSR
DESCRIBE_SAMPLE = 1_000_000  # 诊断统计的采样上限

param_names = ["RNA_IN", "ATAC_IN", "CELLS", "RNA_OUT", "ATAC_OUT",
               "SPARSIFY_FRAC", "DESCRIBE_SAMPLE"]
print("运行参数：")
for _name in param_names:
    print(f"  {_name} = {globals()[_name]}")

cells = pd.Index(pd.read_csv(CELLS, header=None)[0].astype(str)).drop_duplicates()
print(f"目标细胞数：{len(cells)}")


def subset(path, cells, label):
    a = ad.read_h5ad(path, backed="r")
    idx = a.obs_names.get_indexer(cells)
    if (idx < 0).any():
        miss = cells[idx < 0]
        print(f"[{label}] 有 {len(miss)} 个细胞不在文件中，示例：{list(miss[:5])}")
        sys.exit(1)
    # backed 模式下按行读取要求索引递增，先排序读取，再恢复为给定顺序
    sub = a[np.sort(idx)].to_memory()
    a.file.close()
    sub = sub[cells].copy()
    # 全量邻接图与 uns 中的聚类信息对子集已无意义，删除以减小文件并避免读取兼容问题
    for k in list(sub.obsp.keys()):
        del sub.obsp[k]
    sub.uns = {}
    return sub


def describe(sub, label):
    X = sub.X
    print(f"\n[{label}] shape={sub.shape}, 稀疏={sp.issparse(X)}, dtype={X.dtype}")
    vals = X.data if sp.issparse(X) else np.asarray(X).ravel()
    s = vals[:DESCRIBE_SAMPLE]
    print(f"  非零值 min={s.min():.4g}, max={s.max():.4g}, 是否全为整数={np.allclose(s, np.round(s))}")
    nnz = X.nnz if sp.issparse(X) else np.count_nonzero(X)
    print(f"  非零比例={nnz / (X.shape[0] * X.shape[1]):.4f}")
    print(f"  前 5 个特征名：{list(sub.var_names[:5])}")
    print(f"  obs 列：{list(sub.obs.columns)}")
    print(f"  obsm：{list(sub.obsm.keys())}")
    print(f"  已有 layers：{list(sub.layers.keys())}")


def prepare(sub, label):
    X = sub.X
    if not sp.issparse(X):
        frac = np.count_nonzero(X) / X.size
        if frac < SPARSIFY_FRAC:
            X = sp.csr_matrix(X)
            print(f"[{label}] 稠密矩阵非零比例 {frac:.3f}，已转为稀疏 CSR")
        else:
            print(f"[{label}] 警告：矩阵稠密且非零比例 {frac:.3f}，R 中内存占用约为 "
                  f"{X.shape[0] * X.shape[1] * 8 / 1e9:.1f} GB 的双精度矩阵")
    sub.X = X
    # 与 schard 读取方式保持一致：把 X 复制到 counts layer
    sub.layers["counts"] = X
    return sub


rna = subset(RNA_IN, cells, "RNA")
atac = subset(ATAC_IN, cells, "ATAC")
assert (rna.obs_names == atac.obs_names).all(), "两个子集的细胞顺序不一致"

describe(rna, "RNA")
describe(atac, "ATAC")

rna = prepare(rna, "RNA")
atac = prepare(atac, "ATAC")

rna.write_h5ad(RNA_OUT)
atac.write_h5ad(ATAC_OUT)
print(f"\n已写出：\n  {RNA_OUT}\n  {ATAC_OUT}")