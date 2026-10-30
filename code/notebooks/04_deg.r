# 04_deg.r
# 目的：对每个疾病，使用来源和年龄匹配的 Healthy donor 计算差异基因。
# 方法：1) donor 级 1:1 匹配（每个疾病独立匹配，年龄差不超过 15 岁）；
#          优先在相同 TITLE 内匹配，再为剩余病例跨 TITLE 补配；
#       2) 单细胞 Wilcoxon（仅纳入当前疾病的配对 donor，不在检验中使用 pair_id）；
#       3) donor 级 pseudobulk DESeq2（配对设计 ~ pair_id + group）。
# 输出：results/deg_age_matched/*.csv

# ---------------- 参数 ----------------
input_rds <- "data/processed/Atlas_20260721/output/astro_e2000_subFrontal.rds"
output_dir <- "results/deg_age_matched"
diseases <- c("ALS", "FTLD", "AD", "ASD", "COVID-19", "PD", "MDD", "MS", "PTSD", "TDP")
max_age_diff <- 10
min_pairs_deseq2 <- 2L

# ---------------- 输入与校验 ----------------
seu_f <- readRDS(input_rds)

sample_meta <- seu_f@meta.data |>
  dplyr::distinct(batch, DISEASE, AGE, Gender, TITLE) |>
  dplyr::rename(
    sample_id = batch,
    age = AGE,
    sex = Gender,
    title = TITLE
  ) |>
  dplyr::mutate(age = suppressWarnings(as.numeric(age)))

if (anyDuplicated(sample_meta$sample_id)) {
  stop("Some sample_id values have inconsistent donor metadata.")
}
if (anyNA(sample_meta$age)) {
  stop("AGE contains missing or non-numeric values.")
}

healthy <- sample_meta |>
  dplyr::filter(DISEASE == "Healthy")

# ---------------- 来源优先的年龄配对 ----------------
# 在给定候选集中最大化配对数，并在此基础上最小化总年龄差。
match_by_age <- function(cases, controls, max_age_diff) {
  nr <- nrow(cases)
  nc <- nrow(controls)
  if (nr == 0L || nc == 0L) {
    return(tibble::tibble(
      case_sample = character(),
      control_sample = character(),
      case_age = double(),
      control_age = double(),
      case_title = character(),
      control_title = character(),
      same_title = logical(),
      age_diff = double()
    ))
  }

  age_dist <- abs(outer(cases$age, controls$age, `-`))
  unmatched_cost <- max_age_diff + 1
  forbidden_cost <- 1e6
  cost <- matrix(0, nrow = nr + nc, ncol = nc + nr)
  cost[seq_len(nr), seq_len(nc)] <- ifelse(
    age_dist <= max_age_diff,
    age_dist,
    forbidden_cost
  )
  cost[seq_len(nr), nc + seq_len(nr)] <- unmatched_cost
  cost[nr + seq_len(nc), seq_len(nc)] <- unmatched_cost

  assignment <- clue::solve_LSAP(cost)
  control_index <- as.integer(assignment[seq_len(nr)])
  candidate_index <- pmin(control_index, nc)
  keep <- control_index <= nc &
    age_dist[cbind(seq_len(nr), candidate_index)] <= max_age_diff

  tibble::tibble(
    case_sample = cases$sample_id[keep],
    control_sample = controls$sample_id[control_index[keep]],
    case_age = cases$age[keep],
    control_age = controls$age[control_index[keep]],
    case_title = cases$title[keep],
    control_title = controls$title[control_index[keep]],
    same_title = case_title == control_title,
    age_diff = abs(case_age - control_age)
  )
}

match_disease <- function(disease) {
  cases <- sample_meta |>
    dplyr::filter(.data$DISEASE == .env$disease)

  if (nrow(cases) == 0L) {
    return(tibble::tibble())
  }

  # 第一阶段：各 TITLE 内独立匹配，保证同来源优先。
  same_title_pairs <- purrr::map_dfr(unique(cases$title), function(current_title) {
    match_by_age(
      cases |> dplyr::filter(.data$title == .env$current_title),
      healthy |> dplyr::filter(.data$title == .env$current_title),
      max_age_diff = max_age_diff
    )
  })

  # 第二阶段：仅对第一阶段未匹配的病例和未使用的 Healthy 跨 TITLE 补配。
  remaining_cases <- cases |>
    dplyr::filter(!.data$sample_id %in% same_title_pairs$case_sample)
  remaining_controls <- healthy |>
    dplyr::filter(!.data$sample_id %in% same_title_pairs$control_sample)
  cross_title_pairs <- match_by_age(
    remaining_cases,
    remaining_controls,
    max_age_diff = max_age_diff
  )

  dplyr::bind_rows(same_title_pairs, cross_title_pairs) |>
    dplyr::arrange(dplyr::desc(.data$same_title), .data$case_title, .data$case_sample) |>
    dplyr::mutate(
      disease = .env$disease,
      pair_id = paste0(.env$disease, "_pair_", dplyr::row_number()),
      max_age_diff_allowed = max_age_diff,
      .before = 1
    )
}

all_pairs <- purrr::map_dfr(diseases, match_disease)

# ---------------- 每种疾病独立分析 ----------------
run_disease_de <- function(disease) {
  pairs <- all_pairs |>
    dplyr::filter(.data$disease == .env$disease)

  if (nrow(pairs) == 0L) {
    return(list(
      single = tibble::tibble(),
      pseudo = tibble::tibble(),
      status = tibble::tibble(
        disease = disease,
        findmarkers_status = "not_run: no matched pairs",
        pseudobulk_status = "not_run: no matched pairs"
      )
    ))
  }

  case_samples <- pairs$case_sample
  control_samples <- pairs$control_sample
  matched_samples <- c(case_samples, control_samples)

  # 只纳入当前疾病配对表中的 donor；不复用其他疾病的 Healthy donor。
  seu_sub <- subset(seu_f, subset = batch %in% matched_samples)
  seu_sub$group <- factor(
    ifelse(seu_sub$batch %in% case_samples, disease, "Healthy"),
    levels = c("Healthy", disease)
  )
  Seurat::Idents(seu_sub) <- "group"

  single <- Seurat::FindMarkers(
    seu_sub,
    ident.1 = disease,
    ident.2 = "Healthy",
    test.use = "wilcox",
    logfc.threshold = 0,
    min.pct = 0.1,
    max.cells.per.ident = 2000,
    only.pos = FALSE
  ) |>
    tibble::rownames_to_column("gene") |>
    dplyr::mutate(
      disease = disease,
      comparison = paste0(disease, "_vs_Healthy"),
      n_pairs = nrow(pairs),
      .before = 1
    ) |>
    dplyr::arrange(.data$p_val_adj, dplyr::desc(abs(.data$avg_log2FC)))

  if (nrow(pairs) < min_pairs_deseq2) {
    return(list(
      single = single,
      pseudo = tibble::tibble(),
      status = tibble::tibble(
        disease = disease,
        findmarkers_status = "completed",
        pseudobulk_status = paste0("not_run: fewer than ", min_pairs_deseq2, " pairs")
      )
    ))
  }

  counts <- Seurat::GetAssayData(seu_sub, assay = "RNA", layer = "counts")
  sample_ids <- as.character(seu_sub$batch)
  bulk <- vapply(
    matched_samples,
    function(sample_id) {
      Matrix::rowSums(counts[, sample_ids == sample_id, drop = FALSE])
    },
    numeric(nrow(counts))
  )
  rownames(bulk) <- rownames(counts)

  pair_map <- pairs |>
    tidyr::pivot_longer(
      cols = c("case_sample", "control_sample"),
      names_to = "role",
      values_to = "sample_id"
    ) |>
    dplyr::transmute(
      sample_id,
      pair_id = factor(pair_id),
      group = factor(
        ifelse(role == "case_sample", .env$disease, "Healthy"),
        levels = c("Healthy", .env$disease)
      )
    ) |>
    dplyr::arrange(match(sample_id, matched_samples)) |>
    tibble::column_to_rownames("sample_id")

  bulk <- bulk[, rownames(pair_map), drop = FALSE]
  keep_gene <- Matrix::rowSums(bulk) > 0
  bulk <- bulk[keep_gene, , drop = FALSE]

  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = round(bulk),
    colData = pair_map,
    design = ~ pair_id + group
  )
  dds <- DESeq2::DESeq(dds, quiet = TRUE)
  coefficient <- grep("^group_.*_vs_Healthy$", DESeq2::resultsNames(dds), value = TRUE)
  if (length(coefficient) != 1L) {
    stop("Could not identify the disease-vs-Healthy DESeq2 coefficient for ", disease, ".")
  }

  result_raw <- DESeq2::results(dds, name = coefficient) |>
    as.data.frame() |>
    tibble::rownames_to_column("gene")

  result_shrunk <- DESeq2::lfcShrink(
    dds,
    coef = coefficient,
    type = "apeglm",
    svalue = TRUE,
    quiet = TRUE
  )

  pseudo <- result_raw |>
    dplyr::mutate(
      log2FC_shrunk = result_shrunk$log2FoldChange[match(.data$gene, rownames(result_shrunk))],
      svalue = result_shrunk$svalue[match(.data$gene, rownames(result_shrunk))],
      disease = disease,
      comparison = paste0(disease, "_vs_Healthy"),
      n_pairs = nrow(pairs),
      .before = 1
    ) |>
    tibble::as_tibble() |>
    dplyr::arrange(.data$padj, dplyr::desc(abs(.data$log2FoldChange)))

  list(
    single = single,
    pseudo = pseudo,
    status = tibble::tibble(
      disease = disease,
      findmarkers_status = "completed",
      pseudobulk_status = "completed"
    )
  )
}

results_by_disease <- stats::setNames(lapply(diseases, run_disease_de), diseases)
deg_single <- lapply(results_by_disease, `[[`, "single")
deg_pseudo <- lapply(results_by_disease, `[[`, "pseudo")
deg_single_combined <- dplyr::bind_rows(deg_single)
deg_pseudo_combined <- dplyr::bind_rows(deg_pseudo)
status_tab <- dplyr::bind_rows(lapply(results_by_disease, `[[`, "status"))

# ---------------- 保存结果 ----------------
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

cell_counts <- table(seu_f$batch)
all_pairs <- all_pairs |>
  dplyr::mutate(
    n_cells_case = as.integer(cell_counts[match(.data$case_sample, names(cell_counts))]),
    n_cells_control = as.integer(cell_counts[match(.data$control_sample, names(cell_counts))])
  )

summary_tab <- tibble::tibble(disease = diseases) |>
  dplyr::left_join(
    all_pairs |>
      dplyr::group_by(.data$disease) |>
      dplyr::summarise(
        n_pairs = dplyr::n(),
        n_same_title_pairs = sum(.data$same_title),
        pct_same_title_pairs = mean(.data$same_title),
        median_age_diff = stats::median(.data$age_diff),
        max_age_diff = max(.data$age_diff),
        n_cells_case = sum(.data$n_cells_case),
        n_cells_healthy_matched = sum(.data$n_cells_control),
        .groups = "drop"
      ),
    by = "disease"
  ) |>
  dplyr::mutate(
    n_sig_findmarkers = purrr::map_int(
      .data$disease,
      ~ sum(deg_single[[.x]]$p_val_adj < 0.05, na.rm = TRUE)
    ),
    n_sig_pseudobulk = purrr::map_int(
      .data$disease,
      ~ sum(deg_pseudo[[.x]]$padj < 0.05, na.rm = TRUE)
    )
  ) |>
  dplyr::left_join(status_tab, by = "disease")

readr::write_csv(all_pairs, file.path(output_dir, "matched_pairs_all_diseases.csv"))
readr::write_csv(deg_single_combined, file.path(output_dir, "deg_findmarkers_wilcox_all_diseases.csv"))
readr::write_csv(deg_pseudo_combined, file.path(output_dir, "deg_pseudobulk_deseq2_all_diseases.csv"))
readr::write_csv(summary_tab, file.path(output_dir, "summary.csv"))
readr::write_csv(status_tab, file.path(output_dir, "status.csv"))

summary_tab
