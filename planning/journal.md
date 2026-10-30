# parameters
configs_fxl.toml
香玲师姐模型的参数

fxl_gene.csv
香玲师姐做MIDAS整合时scRNA部分用的Gene。计算的高变基因

gene_datasets.csv
182个数据集的Gene的并集

gene_paper.csv
从数据集来源的paper中扒下来的Gene

gene_20260221.csv
fxl_gene.csv + gene_paper.csv，用作2026.02.21开始训练的scRNA部分的feats

Gemini_KIMI_GPT.csv
利用Gemini、KIMI、ChatGPT得到的GENE

gene_20260302.csv
包含以下几部分
1. Gemini_KIMI_GPT
2. gene_datasets中count大于158的基因
3. gene_20260221
4. @Ageing promotes microglial accumulation of slow-degrading synaptic proteins 中提到的1726个蛋白
用作2026.03.02开始训练的scRNA部分的feats(实际训练的开始时间是2026.03.03)

