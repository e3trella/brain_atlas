import torch

# 1. 加载原始文件

def rename_and_save_model(input_path, output_path):
    ckpt = torch.load(input_path, map_location='cpu')

    old_sd = ckpt["net"]   # 重点在这里
    new_sd = {}

    # 重命名 ATAC
    rename_map = {
        "encoder.after_concat.atac.net.0.weight": "encoder.transform_concat.atac.net.0.weight",
        "encoder.after_concat.atac.net.0.bias":   "encoder.transform_concat.atac.net.0.bias",
        "decoder.before_concat.atac.net.0.weight": "decoder.transform_concat.atac.net.0.weight",
        "decoder.before_concat.atac.net.0.bias":   "decoder.transform_concat.atac.net.0.bias",
    }

    # 丢弃 RNA
    drop_keys = {
        "decoder.before_concat.rna.net.0.weight",
        "decoder.before_concat.rna.net.0.bias",
    }

    for k, v in old_sd.items():
        # 删除
        if k in drop_keys:
            continue
        
        # 重命名
        if k in rename_map:
            new_sd[rename_map[k]] = v
            continue
        
        # 其他不动
        new_sd[k] = v

    # 更新回 checkpoint
    ckpt["net"] = new_sd

    # 保存
    torch.save(ckpt, output_path)

