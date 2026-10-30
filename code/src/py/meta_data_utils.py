import pandas as pd
from mudata import MuData

def merge_metadata(obs: pd.DataFrame, df: pd.DataFrame, on: str, suffixes=('_orig', '_new')):
    """
    将外部 DataFrame 合并到 mdata.obs 中，并保证细胞顺序和索引不乱。
    
    参数:
    ----------
    obs : pd.DataFrame
        mdata.obs 数据框。
    df : pd.DataFrame
        包含新信息的 DataFrame。
    on : str
        合并的键(Key), 这一列必须同时存在于 obs 和 df 中。
        例如：'sample_id' 或 'patient_id'。
    suffixes : tuple
        如果列名冲突，加的后缀。
        
    返回:
    ----------
    pd.DataFrame : 更新后的 obs 数据框
    """
    
    # 1. 检查 Key 是否存在
    if on not in obs.columns:
        raise ValueError(f"错误：列 '{on}' 不在 obs 中。")
    if on not in df.columns:
        raise ValueError(f"错误：列 '{on}' 不在输入的 DataFrame 中。")
    # 2. 检查 df 中是否有重复的 on 值
    duplicated = df[on].duplicated()
    if duplicated.any():
        dup_values = df.loc[duplicated, on].unique()
        raise ValueError(
            f"警告：df 中的列 '{on}' 包含重复值: {dup_values[:5]}... "
            f"这可能导致数据合并不确定。请先去重。"
        )
    # 2. 准备工作：保留原始的 Cell Barcode (Index)
    # 这是一个关键步骤，防止 merge 后 index 丢失
    original_index = obs.index.copy()
    
    # 临时将 index 变成一列，方便 merge 后找回
    obs_temp = obs.copy()
    obs_temp['__original_barcode_index__'] = obs_temp.index
    
    # 3. 执行合并 (Left Join)
    # how='left' 保证 mdata 里的细胞一个都不会少
    merged_obs = pd.merge(
        obs_temp, 
        df, 
        on=on, 
        how='left', 
        suffixes=suffixes
    )
    # 5. 【新增】检查是否有未匹配的行
    on_null = merged_obs[on].isna()
    if on_null.any():
        print(f"注意：{on_null.sum()} 个细胞在 df 中没有匹配到 {on} 值")
    # 4. 恢复 Index 并核对顺序
    # 重新设置 Index
    merged_obs = merged_obs.set_index('__original_barcode_index__')
    # 删除 index 的名字（美观）
    merged_obs.index.name = None
    
    # 【最关键一步】强制按照原始 obs 的顺序重排
    # 如果 df 里有重复匹配导致行数增加，这里会自动去重或报错，保证数据完整性
    # 这里的 reindex 会自动对齐数据，如果有细胞没匹配上，新列会自动填 NaN
    merged_obs = merged_obs.reindex(original_index)
    
    print(f"成功合并!obs 新增列: {[col for col in df.columns if col != on]}")
    
    return merged_obs


# ...existing code...


def optimize_dataframe_to_categorical(df, columns=None, cardinality_threshold=0.5):
    """
    Convert pandas DataFrame columns to categorical type with memory optimization report.
    
    Parameters:
    -----------
    df : pd.DataFrame
        The input dataframe (e.g., mdata.obs).
    columns : list, optional
        Specific columns to convert. If None, auto-selects 'object' and 'string' columns.
    cardinality_threshold : float
        Only convert if (unique values / total rows) < threshold. 
        High cardinality columns (like unique Cell IDs) should stay as strings.
    """
    df_optimized = df.copy()
    
    # 计算初始内存占用 (Deep memory usage for accurate string size)
    start_mem = df_optimized.memory_usage(deep=True).sum() / 1024**2
    
    # 自动识别需要转换的列
    if columns is None:
        columns = df_optimized.select_dtypes(include=['object', 'string']).columns.tolist()
    
    print(f"--- Starting Optimization for {len(columns)} columns ---")
    
    for col in columns:
        n_unique = df_optimized[col].nunique()
        n_total = len(df_optimized)
        ratio = n_unique / n_total
        
        # 批判性评估：如果唯一值太多（如接近 1.2M），转换反而会增加负担
        if ratio < cardinality_threshold:
            # 显式转换为 category，处理空值并保留原始索引
            df_optimized[col] = df_optimized[col].astype('category')
            # print(f"Column '{col}': Converted. (Unique: {n_unique})")
        else:
            print(f"Column '{col}': Skipped. High cardinality (Ratio: {ratio:.2f})")

    # 计算优化后内存占用
    end_mem = df_optimized.memory_usage(deep=True).sum() / 1024**2
    savings = (1 - end_mem / start_mem) * 100
    
    print("--- Optimization Complete ---")
    print(f"Memory Usage: {start_mem:.2f} MB -> {end_mem:.2f} MB ({savings:.1%}- savings)")
    
    return df_optimized

# 使用示例
# mdata.obs = optimize_dataframe_to_categorical(mdata.obs)
import pandas as pd
import numpy as np

def update_metadata_custom(
    df, 
    subset_col,
    subset_list, 
    target_col, 
    source_col=None, 
    mapping_dict=None,
    overwrite_na=False, 
    inplace=False
):
    """
    针对大规模单细胞元数据进行批量、条件化填充的高效函数。    
    参数:
    ----
    df : pd.DataFrame
        输入的 obs 数据框。
    subset_col : str
        用于过滤的列名，例如 'SUBSET_ID'。
    subset_list : list
        需要修改的子集列表 (例如 ['subset_181', 'subset_182'])。
    target_col : str
        目标修改列名 (例如 'REGION_L4')。
    source_col : str, optional
        来源列名 (例如 'tissue')。如果不提供 mapping_dict，则从此列提取值。
    mapping_dict : dict, optional
        多对一映射需求。例如 {'subset_51': 'Lung', 'subset_52': 'Liver'}。
        如果提供此参数，将优先于 source_col 执行。
    overwrite_na : bool, default False
        若来源值为 NA，是否覆盖目标列中的已有值。
    inplace : bool, default False
        是否在原位修改。
        
    返回:
    ----
    pd.DataFrame (如果 inplace=False) 或 None
    """
    
    # 1. 内存模式选择
    working_df = df if inplace else df.copy()
    
    # 防御性编程：确保目标列存在
    if target_col not in working_df.columns:
        working_df[target_col] = np.nan

    # 2. 构建基础布尔掩码 (Mask)
    subset_mask = working_df[subset_col].isin(subset_list)
    if not subset_mask.any():
        print(f"⚠️ 未找到符合条件的行 (Subset: {subset_list})。操作已跳过。")
        return None if inplace else working_df

    # 3. 提取来源数据
    if mapping_dict:
        # map 操作在底层由 C 语言实现，速度极快
        fill_values = working_df.loc[subset_mask, subset_col].map(mapping_dict)
        source_info = "mapping_dict"
    elif source_col:
        fill_values = working_df.loc[subset_mask, source_col]
        source_info = f"column '{source_col}'"
    else:
        raise ValueError("必须提供 'source_col' 或 'mapping_dict' 其中之一。")

    # 4. 高效处理 NA 与掩码对齐 (核心优化区)
    if not overwrite_na:
        valid_mask = fill_values.notna()
        # 更新掩码：原在子集内，且来源值不为 NA。使用底层布尔赋值避免耗时的 Index 匹配
        final_mask = subset_mask.copy()
        final_mask[subset_mask] = valid_mask 
        
        final_fill_values = fill_values[valid_mask]
        update_info = f"保留目标列原值 (跳过 {(~valid_mask).sum()} 个 NA 源值)"
    else:
        final_mask = subset_mask
        final_fill_values = fill_values
        update_info = f"执行覆盖 (包含 {fill_values.isna().sum()} 个 NA 源值)"

    # 5. Categorical 安全降级与扩展 (彻底杜绝类型碰撞)
    # [降维打击]：剥离来源数据的 Categorical 外壳，保留纯 object 与 np.nan
    if isinstance(final_fill_values.dtype, pd.CategoricalDtype):
        safe_fill_values = final_fill_values.astype(object)
    else:
        safe_fill_values = final_fill_values

    # [目标扩展]：如果目标列是 Categorical，必须先注册新标签
    if isinstance(working_df[target_col].dtype, pd.CategoricalDtype):
        # 提取去重后的新标签，排除 NaN
        new_labels = pd.Series(safe_fill_values).dropna().unique()
        existing_cats = working_df[target_col].cat.categories
        to_add = [val for val in new_labels if val not in existing_cats]
        
        if to_add:
            working_df[target_col] = working_df[target_col].cat.add_categories(to_add)

    # 6. 底层内存连续赋值 (消除 OOM 风险)
    # 强制使用 .values 提取 NumPy 数组，进行精确的等长位置替换
    working_df.loc[final_mask, target_col] = safe_fill_values.values

    # 7. 统计与质控输出
    print(f"✅ 任务完成: [{target_col}]")
    print(f"   - 来源: {source_info} | 策略: {update_info}")
    print(f"   - 更新规模: 精确修改了 {final_mask.sum()} 个细胞/行。")

    if not inplace:
        return working_df



def read_cell_names(subset_paths):
    """
    读取多个 subset_ 的 cell_names.csv 文件，并合并成一个 DataFrame。
    
    参数:
    ----
    subset_paths : list of str or Path
        包含所有 subset_ 的 cell_names.csv 文件路径的列表。
        
    返回:
    ----
    pd.DataFrame
        合并后的 DataFrame，包含所有细胞的名称和对应的 subset_ 信息。
    """
    df_list = []
    
    for path in subset_paths:
        df = pd.read_csv(path, index_col=1)  # 假设细胞名称在第二列
        df_list.append(df)
    
    combined_df = pd.concat(df_list, ignore_index=False)  # 保持原索引不变
    
    return combined_df