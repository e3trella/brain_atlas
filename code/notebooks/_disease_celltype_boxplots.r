# 疾病 × 细胞类型细胞比例箱线图
#
# 功能：
#   1. 读取纯免疫细胞 WNN Seurat 对象的 meta.data，按 batch 统计 mg_level3
#      各细胞类型的细胞数，并计算每个细胞类型在该 batch 全部细胞中的比例
#      （分母为该 batch 全部目标细胞，先算比例再按疾病筛选）。
#   2. 剔除细胞总数少于 50 的 batch；把 Healthy 重编码为 HC，并把 DISEASE
#      固定为 factor 顺序（HC 在最前，其后为各疾病）。
#   3. 通过命令行参数 --disease 与 --celltype（逗号分隔，可传多个）选择要画
#      的疾病与细胞类型；默认画全部疾病和全部细胞类型。
#   4. 对每个疾病 × 细胞类型组合，画出所有包含该疾病的 TITLE 分组的箱线图
#      （x 轴显示该 TITLE 对应的疾病名，分面使用 TITLE 的简称），并在每个
#      分面内用 Wilcoxon 检验做两两比较、标注未校正的显著性。
#   5. 每张图输出到 results/plots/{disease}_{celltype}.pdf。
#
# 用法示例：
#   Rscript code/notebooks/_disease_celltype_boxplots.r --disease AD,PD --celltype Microglia

library(ggplot2)

plot_dir <- "results/plots"
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

obj <- readRDS("data/processed/Atlas_20260721/output/wnn/pure_imm/candidate_pure_imm_wnn_seurat.rds")

title_short <- c(
  "A comparative atlas of single-cell chromatin accessibility in the human brain" = "Chromatin atlas",
  "A lymphocyte–microglia–astrocyte axis in chronic active multiple sclerosis" = "Chronic active MS",
  "A single nuclear transcriptomic characterisation of mechanisms responsible for impaired angiogenesis and blood-brain barrier function in alzheimer’s disease" = "AD angiogen.",
  "A single-nuclei paired multiomic analysis of the human midbrain reveals age- and parkinson’s disease–associated glial changes" = "Midbrain aging & PD",
  "Altered human oligodendrocyte heterogeneity in multiple sclerosis" = "MS oligodendrocyte",
  "Cell type specific transcriptomic differences in depression show similar patterns between males and females but implicate distinct cell types and genes" = "MDD sex diff",
  "Dysregulation of brain and choroid plexus cell types in severe COVID-19" = "Severe COVID-19",
  "Human striatal glia differentially contribute to AD- and PD-specific neurodegeneration" = "Striatal glia",
  "Huntington disease oligodendrocyte maturation deficits revealed by single-nucleus RNAseq are rescued by thiamine-biotin supplementation" = "Huntington OL",
  "Neuronal vulnerability and multilineage diversity in multiple sclerosis" = "MS neuronal",
  "Single cell transcriptome profiling of the human alcohol-dependent brain" = "Alcohol dependence",
  "Single nucleus multiomics identifies ZEB1 and MAFB as candidate regulators of alzheimer’s disease-specific cis-regulatory elements" = "AD ZEB1/MAFB",
  "Single nucleus transcriptomics of ventral midbrain identifies glial activation associated with chronic opioid use disorder" = "Opioid use disorder",
  "Single-cell dissection of the human motor and prefrontal cortices in ALS and FTLD" = "ALS/FTLD cortex",
  "Single-cell genomics identifies cell type-specific molecular changes in autism" = "ASD",
  "Single-cell sequencing of human midbrain reveals glial activation and a Parkinson-specific neuronal state" = "Midbrain PD",
  "Single-cell transcriptomic and chromatin dynamics of the human brain in PTSD" = "PTSD",
  "Single-cell transcriptomic and proteomic analysis of Parkinson’s disease brains" = "PD omics",
  "Single-nucleus RNA-seq of normal-appearing brain regions in relapsing-remitting vs Secondary progressive multiple sclerosis: Implications for the efficacy of fingolimod" = "RRMS vs SPMS",
  "Single-nucleus multiomics reveals the disrupted regulatory programs in three brain regions of sporadic early-onset Alzheimer’s disease" = "Early-onset AD",
  "Single-nucleus transcriptomics of the prefrontal cortex in major depressive disorder implicates oligodendrocyte precursor cells and excitatory neurons" = "MDD prefrontal cortex",
  "Unravelling cell type-specific responses to Parkinson's Disease at single cell resolution" = "PD cell types",
  "pTDP-43 levels correlate with cell type–specific molecular alterations in the prefrontal cortex of C9orf72 ALS/FTD patients" = "C9orf72 ALS/FTD"
)

diseases <- c("AD", "TDP", "FTLD", "PD", "Huntington", "ALS", "MS",
              "COVID-19", "Alcohol", "Drug", "MDD", "PTSD", "ASD")
disease_levels <- c("HC", diseases)

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default) {
  i <- match(paste0("--", flag), args)
  if (!is.na(i) && length(args) >= i + 1) return(strsplit(args[i + 1], ",")[[1]])
  eq <- grep(paste0("^--", flag, "="), args, value = TRUE)
  if (length(eq)) return(strsplit(sub(paste0("^--", flag, "="), "", eq[1]), ",")[[1]])
  default
}

meta <- obj@meta.data |>
  dplyr::filter(!is.na(TITLE), !is.na(mg_level3), DISEASE %in% c("Healthy", diseases)) |>
  dplyr::count(batch, TITLE, DISEASE, mg_level3, name = "n") |>
  dplyr::group_by(batch) |>
  dplyr::filter(sum(n) >= 50) |>
  dplyr::mutate(prop = n / sum(n)) |>
  dplyr::ungroup() |>
  dplyr::mutate(
    TITLE_short = unname(title_short[TITLE]),
    DISEASE = factor(ifelse(DISEASE == "Healthy", "HC", as.character(DISEASE)),
                     levels = disease_levels)
  )

all_cell_types <- sort(unique(meta$mg_level3))
sel_diseases <- get_arg("disease", diseases)
sel_cell_types <- get_arg("celltype", all_cell_types)

if (length(bad <- setdiff(sel_diseases, diseases))) {
  stop("Unknown disease: ", paste(bad, collapse = ", "),
       "\nValid: ", paste(diseases, collapse = ", "), call. = FALSE)
}
if (length(bad <- setdiff(sel_cell_types, all_cell_types))) {
  stop("Unknown cell type: ", paste(bad, collapse = ", "),
       "\nValid: ", paste(all_cell_types, collapse = ", "), call. = FALSE)
}

combos <- expand.grid(
  disease = sel_diseases,
  cell_type = sel_cell_types,
  stringsAsFactors = FALSE
)

for (i in seq_len(nrow(combos))) {
  dis <- combos$disease[i]
  ct <- combos$cell_type[i]

  titles <- unique(meta$TITLE[meta$DISEASE == dis])
  if (!length(titles)) next

  d <- meta |>
    dplyr::filter(TITLE %in% titles, mg_level3 == ct)

  pw <- d |>
    dplyr::group_by(TITLE) |>
    dplyr::group_modify(function(x, key) {
      lv <- levels(droplevels(x$DISEASE))
      if (length(lv) < 2) return(tibble::tibble())
      tryCatch(
        rstatix::wilcox_test(x, prop ~ DISEASE) |>
          rstatix::add_significance("p") |>
          dplyr::mutate(xmin = match(group1, lv),
                        xmax = match(group2, lv),
                        y = max(x$prop) + 0.10 * max(x$prop) * dplyr::row_number()),
        error = function(e) tibble::tibble()
      )
    }) |>
    dplyr::ungroup() |>
    dplyr::left_join(dplyr::distinct(d, TITLE, TITLE_short), by = "TITLE")

  p <- ggplot2::ggplot(d, ggplot2::aes(x = DISEASE, y = prop)) +
    ggplot2::geom_boxplot(outlier.shape = NA, width = 0.6) +
    ggplot2::geom_jitter(width = 0.12, size = 0.7, alpha = 0.4)

  if (nrow(pw) > 0) {
    p <- p + ggpubr::stat_pvalue_manual(
      pw, xmin = "xmin", xmax = "xmax", y.position = "y",
      label = "p.signif", tip.length = 0.01, size = 4, vjust = 0.4
    )
  }

  ncol <- min(6, length(titles))
  nrow <- ceiling(length(titles) / ncol)

  p <- p +
    ggplot2::facet_wrap(~ TITLE_short, ncol = ncol, scales = "free",
                        labeller = ggplot2::label_wrap_gen(10)) +
    ggplot2::scale_y_continuous(
      labels = function(x) formatC(x, format = "g", digits = 2),
      expand = ggplot2::expansion(mult = c(0.02, 0.20))
    ) +
    ggplot2::labs(x = NULL, y = sprintf("Proportion of %s per batch", ct),
                  title = sprintf("%s - %s", dis, ct)) +
    ggplot2::theme_bw(base_size = 15) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5, size = 11),
      axis.text.y = ggplot2::element_text(size = 9),
      axis.title.y = ggplot2::element_text(size = 13),
      strip.text = ggplot2::element_text(size = 10, lineheight = 0.85),
      panel.spacing.x = grid::unit(0.6, "lines")
    )

  out_file <- file.path(
    plot_dir,
    sprintf("%s_%s.pdf", gsub("[^A-Za-z0-9_.-]", "_", dis),
            gsub("[^A-Za-z0-9_.-]", "_", ct))
  )
  ggplot2::ggsave(out_file, p, device = "pdf",
                  width = 2.8 * ncol + 0.8, height = 3.2 * nrow + 1.0)
  message(sprintf("[%d/%d] %s", i, nrow(combos), out_file))
}
