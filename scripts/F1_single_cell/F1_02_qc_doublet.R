# F1.2 固定QC与doublet ------------------------------------------------------
#
# 生物学目的：排除明显低质量细胞和主算法预测双细胞，同时保留每个公开细胞的
# 决策轨迹。环境RNA评估要等F1.4获得可靠粗谱系标签后再按样本运行DecontX。
#
# 主要输入：01_all_cells_raw_or_initial.rds。
# 主要输出：02_all_cells_qc_filtered.rds、逐细胞/逐样本QC和doublet摘要。

file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), winslash = "/"))
source(file.path(script_dir, "F1_config.R"), encoding = "UTF-8")
source(file.path(script_dir, "F1_utils.R"), encoding = "UTF-8")

args <- f1_parse_args()
config <- f1_build_config(args$project_root)
stage <- "F1.2"
outputs <- c(
  config$paths$object_02,
  file.path(config$paths$qc_dir, "qc_cell_decision_audit.tsv"),
  file.path(config$paths$qc_dir, "sample_qc_summary.tsv"),
  file.path(config$paths$qc_dir, "doublet_summary_by_sample.tsv")
)

if (!args$execute) {
  f1_stage_dry_run(
    stage,
    "逐样本执行冻结QC、scDblFinder主判定、DoubletFinder敏感性，并写出仅去nCount下界的预注册QC敏感性mask",
    c(config$paths$object_01),
    outputs,
    config$packages[[stage]]
  )
  quit(save = "no", status = 0)
}

f1_require_packages(config$packages[[stage]], stage)
f1_check_f0_ready(config)
if (!file.exists(config$paths$object_01)) stop("请先完成F1.1：", config$paths$object_01)
f1_prepare_directories(config)
set.seed(config$seed)
if (!isTRUE(config$qc$no_ncount_sensitivity_enabled)) {
  stop("预注册的no-nCount QC敏感性mask被关闭；请先恢复F1_config.R中的冻结设置。")
}
f1_append_log(config, stage, "开始逐样本固定QC和doublet识别")

get_df_function <- function(candidates) {
  ns <- asNamespace("DoubletFinder")
  hit <- candidates[vapply(candidates, exists, logical(1), envir = ns, inherits = FALSE)]
  if (!length(hit)) stop("当前DoubletFinder版本缺少函数：", paste(candidates, collapse = "/"))
  get(hit[[1]], envir = ns, inherits = FALSE)
}

run_doubletfinder_sensitivity <- function(counts, sample_id, config) {
  n_cells <- ncol(counts)
  if (n_cells < config$doublet$minimum_cells_for_doubletfinder) {
    return(list(
      status = "not_evaluable_small_sample",
      pANN = setNames(rep(NA_real_, n_cells), colnames(counts)),
      class = setNames(rep(NA_character_, n_cells), colnames(counts)),
      pK = NA_real_, nExp = NA_integer_, expected_rate = NA_real_, error = NA_character_
    ))
  }

  tryCatch({
    set.seed(config$seed)
    x <- Seurat::CreateSeuratObject(counts = counts, min.cells = 0, min.features = 0)
    # 下面的LogNormalize/ScaleData/PCA只为DoubletFinder构造其算法所需空间，
    # 不会保存为F1主聚类对象，也不改变F1.3的SCTransform主线。
    x <- Seurat::NormalizeData(x, verbose = FALSE)
    x <- Seurat::FindVariableFeatures(
      x,
      selection.method = "vst",
      nfeatures = min(2000L, nrow(x)),
      verbose = FALSE
    )
    x <- Seurat::ScaleData(x, features = SeuratObject::VariableFeatures(x), verbose = FALSE)
    max_pc <- min(max(config$doublet$doubletfinder_pcs), ncol(x) - 1L, length(SeuratObject::VariableFeatures(x)) - 1L)
    pcs <- seq_len(max_pc)
    x <- Seurat::RunPCA(x, features = SeuratObject::VariableFeatures(x), npcs = max_pc, verbose = FALSE)

    param_sweep <- get_df_function(c("paramSweep", "paramSweep_v3"))
    summarize_sweep <- get_df_function(c("summarizeSweep"))
    find_pk <- get_df_function(c("find.pK"))
    doublet_finder <- get_df_function(c("doubletFinder", "doubletFinder_v3"))

    sweep <- param_sweep(x, PCs = pcs, sct = FALSE)
    sweep_stats <- summarize_sweep(sweep, GT = FALSE)
    pk_table <- find_pk(sweep_stats)
    bc_column <- intersect(c("BCmetric", "BCmvn"), colnames(pk_table))
    if (!length(bc_column)) {
      stop(
        "DoubletFinder pK表缺少BCmetric/BCmvn；实际字段为：",
        paste(colnames(pk_table), collapse = ", ")
      )
    }
    pk_numeric <- suppressWarnings(as.numeric(as.character(pk_table$pK)))
    bc_metric <- suppressWarnings(as.numeric(as.character(pk_table[[bc_column[[1]]]])))
    valid <- is.finite(pk_numeric) & is.finite(bc_metric)
    if (!any(valid)) stop("DoubletFinder BCmetric表没有可用pK。")
    best <- which(valid & bc_metric == max(bc_metric[valid]))
    selected_pk <- min(pk_numeric[best])

    # 上样/回收信息未公开，使用与scDblFinder自动估计相同量级的10x经验率。
    expected_rate <- min(0.20, config$doublet$scdblfinder_dbr_per_1k * n_cells / 1000)
    n_exp <- max(1L, as.integer(round(expected_rate * n_cells)))
    x <- doublet_finder(
      x,
      PCs = pcs,
      pN = config$doublet$doubletfinder_pN,
      pK = selected_pk,
      nExp = n_exp,
      # DoubletFinder 2.0.6以NULL表示首次计算；FALSE会被误作元数据列索引。
      reuse.pANN = NULL,
      sct = FALSE
    )
    meta <- x[[]]
    pann_col <- tail(grep("^pANN", colnames(meta), value = TRUE), 1)
    class_col <- tail(grep("^DF.classifications", colnames(meta), value = TRUE), 1)
    if (!length(pann_col) || !length(class_col)) stop("DoubletFinder未生成pANN/class字段。")
    list(
      status = "completed_sensitivity_only",
      pANN = setNames(as.numeric(meta[[pann_col]]), rownames(meta)),
      class = setNames(as.character(meta[[class_col]]), rownames(meta)),
      pK = selected_pk,
      nExp = n_exp,
      expected_rate = expected_rate,
      error = NA_character_
    )
  }, error = function(e) {
    list(
      status = "failed_sensitivity_only",
      pANN = setNames(rep(NA_real_, n_cells), colnames(counts)),
      class = setNames(rep(NA_character_, n_cells), colnames(counts)),
      pK = NA_real_, nExp = NA_integer_, expected_rate = NA_real_, error = conditionMessage(e)
    )
  })
}

plot_sample_qc <- function(decision, sample_id, out_dir) {
  long <- rbind(
    data.frame(cell_id = decision$cell_id_final, metric = "nCount_RNA", value = decision$nCount_RNA, fixed_qc_pass = decision$fixed_qc_pass),
    data.frame(cell_id = decision$cell_id_final, metric = "nFeature_RNA", value = decision$nFeature_RNA, fixed_qc_pass = decision$fixed_qc_pass),
    data.frame(cell_id = decision$cell_id_final, metric = "mt_percent", value = decision$mt_percent, fixed_qc_pass = decision$fixed_qc_pass),
    data.frame(cell_id = decision$cell_id_final, metric = "HB_percent", value = decision$HB_percent, fixed_qc_pass = decision$fixed_qc_pass)
  )
  p1 <- ggplot2::ggplot(long, ggplot2::aes(x = fixed_qc_pass, y = value, fill = fixed_qc_pass)) +
    ggplot2::geom_violin(scale = "width", trim = TRUE) +
    ggplot2::facet_wrap(~metric, scales = "free_y", ncol = 2) +
    ggplot2::labs(
      title = paste0(sample_id, ": before and after fixed QC"),
      x = "Pass fixed QC",
      y = NULL
    ) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(legend.position = "none")
  p2 <- ggplot2::ggplot(
    decision,
    ggplot2::aes(x = nCount_RNA, y = nFeature_RNA, color = mt_percent)
  ) +
    ggplot2::geom_point(size = 0.35, alpha = 0.55) +
    ggplot2::scale_color_viridis_c() +
    ggplot2::labs(
      title = paste0(sample_id, ": library complexity and mitochondrial fraction"),
      color = "mt_percent"
    ) +
    ggplot2::theme_bw(base_size = 10)
  ggplot2::ggsave(
    file.path(out_dir, paste0(sample_id, "_QC_review.pdf")),
    patchwork::wrap_plots(p1, p2, ncol = 2),
    width = 10,
    height = 8,
    limitsize = FALSE
  )
}

initial <- readRDS(config$paths$object_01)
initial <- f1_join_assay(initial, "RNA")
all_counts <- SeuratObject::LayerData(initial, assay = "RNA", layer = "counts")
raw_signature_before <- f1_counts_signature(all_counts)
sample_ids <- sort(unique(as.character(initial$sample_id)))

sample_objects <- vector("list", length(sample_ids))
names(sample_objects) <- sample_ids
decision_rows <- vector("list", length(sample_ids))
threshold_rows <- vector("list", length(sample_ids))
summary_rows <- vector("list", length(sample_ids))
doublet_rows <- vector("list", length(sample_ids))

for (i in seq_along(sample_ids)) {
  sample_id <- sample_ids[[i]]
  message(sprintf("[%s] 处理 %s (%d/%d)", stage, sample_id, i, length(sample_ids)))
  cells <- colnames(initial)[initial$sample_id == sample_id]
  raw <- all_counts[, cells, drop = FALSE]

  detected_cells <- Matrix::rowSums(raw > 0)
  keep_feature <- detected_cells >= config$qc$min_cells_per_feature
  working <- raw[keep_feature, , drop = FALSE]
  if (!nrow(working)) stop(sample_id, "在min.cells=3后没有可用feature。")

  ncount <- Matrix::colSums(working)
  nfeature <- Matrix::colSums(working > 0)
  mt_genes <- grepl("^MT-", toupper(rownames(working)))
  hb_genes <- toupper(rownames(working)) %in% config$qc$globin_panel
  mt_counts <- if (any(mt_genes)) Matrix::colSums(working[mt_genes, , drop = FALSE]) else rep(0, ncol(working))
  hb_counts <- if (any(hb_genes)) Matrix::colSums(working[hb_genes, , drop = FALSE]) else rep(0, ncol(working))
  mt_percent <- ifelse(ncount > 0, 100 * mt_counts / ncount, NA_real_)
  hb_percent <- ifelse(ncount > 0, 100 * hb_counts / ncount, NA_real_)

  pass_nfeature_low <- nfeature >= config$qc$nfeature_min_inclusive
  pass_nfeature_high <- nfeature < config$qc$nfeature_max_exclusive
  pass_ncount <- ncount > config$qc$ncount_min_exclusive
  pass_mt <- mt_percent <= config$qc$mt_max_inclusive
  pass_hb <- hb_percent < config$qc$hb_max_exclusive
  fixed_pass <- pass_nfeature_low & pass_nfeature_high & pass_ncount & pass_mt & pass_hb
  # 预注册敏感性只去掉nCount>1000；其他四条QC规则保持完全相同。
  # 这里只产出mask，后续平行分析必须在该mask上重新运行doublet及其余步骤。
  no_ncount_sensitivity_pass <-
    pass_nfeature_low & pass_nfeature_high & pass_mt & pass_hb
  if (any(fixed_pass & !no_ncount_sensitivity_pass, na.rm = TRUE)) {
    stop(sample_id, "的主QC通过细胞不是no-nCount敏感性mask的子集。")
  }
  if (!any(fixed_pass, na.rm = TRUE)) stop(sample_id, "应用冻结QC后没有剩余细胞，需人工复核。")

  fixed_cells <- cells[fixed_pass]
  fixed_counts <- working[, fixed_cells, drop = FALSE]
  # 工作对象继续保留min.cells=3定义的feature；算法只跳过过滤后已全零的无信息行。
  algorithm_counts <- fixed_counts[Matrix::rowSums(fixed_counts) > 0, , drop = FALSE]
  sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = algorithm_counts))
  set.seed(config$seed)
  bp_param <- if (
    .Platform$OS.type != "windows" &&
      config$doublet$scdblfinder_workers > 1L
  ) {
    # 每个样本使用独立的可复现随机流；Windows仍回退到SerialParam。
    BiocParallel::MulticoreParam(
      workers = config$doublet$scdblfinder_workers,
      RNGseed = config$seed + i,
      progressbar = TRUE
    )
  } else {
    BiocParallel::SerialParam(progressbar = TRUE)
  }
  sce <- scDblFinder::scDblFinder(
    sce,
    dbr = NULL,
    dbr.per1k = config$doublet$scdblfinder_dbr_per_1k,
    BPPARAM = bp_param,
    verbose = TRUE
  )
  sc_meta <- as.data.frame(SummarizedExperiment::colData(sce))
  if (!all(c("scDblFinder.score", "scDblFinder.class") %in% colnames(sc_meta))) {
    stop(sample_id, "的scDblFinder未返回score/class。")
  }
  sc_score <- setNames(as.numeric(sc_meta$scDblFinder.score), rownames(sc_meta))
  sc_class <- setNames(as.character(sc_meta$scDblFinder.class), rownames(sc_meta))
  if (any(!tolower(sc_class) %in% c("singlet", "doublet"))) {
    stop(sample_id, "的scDblFinder出现未知分类。")
  }

  df_result <- run_doubletfinder_sensitivity(algorithm_counts, sample_id, config)
  retained_cells <- names(sc_class)[tolower(sc_class) == "singlet"]
  retained_counts <- fixed_counts[, retained_cells, drop = FALSE]

  meta <- initial[[]][retained_cells, , drop = FALSE]
  meta$nCount_RNA <- as.numeric(ncount[retained_cells])
  meta$nFeature_RNA <- as.numeric(nfeature[retained_cells])
  meta$mt_percent <- as.numeric(mt_percent[retained_cells])
  meta$HB_percent <- as.numeric(hb_percent[retained_cells])
  meta$scDblFinder_score <- sc_score[retained_cells]
  meta$scDblFinder_class <- sc_class[retained_cells]
  meta$DoubletFinder_pANN <- df_result$pANN[retained_cells]
  meta$DoubletFinder_class <- df_result$class[retained_cells]
  sample_objects[[sample_id]] <- Seurat::CreateSeuratObject(
    counts = retained_counts,
    assay = "RNA",
    project = "GSE183904_F1_QC",
    meta.data = meta,
    min.cells = 0,
    min.features = 0
  )

  sc_score_all <- setNames(rep(NA_real_, length(cells)), cells)
  sc_class_all <- setNames(rep(NA_character_, length(cells)), cells)
  df_pann_all <- setNames(rep(NA_real_, length(cells)), cells)
  df_class_all <- setNames(rep(NA_character_, length(cells)), cells)
  sc_score_all[names(sc_score)] <- sc_score
  sc_class_all[names(sc_class)] <- sc_class
  df_pann_all[names(df_result$pANN)] <- df_result$pANN
  df_class_all[names(df_result$class)] <- df_result$class
  failure_reasons <- f1_collapse_reasons(
    nFeature_RNA_below_500 = as.list(pass_nfeature_low),
    nFeature_RNA_at_least_6000 = as.list(pass_nfeature_high),
    nCount_RNA_at_most_1000 = as.list(pass_ncount),
    mt_percent_above_20 = as.list(pass_mt),
    HB_percent_at_least_5 = as.list(pass_hb)
  )
  final_include <- fixed_pass & !is.na(sc_class_all[cells]) & tolower(sc_class_all[cells]) == "singlet"
  final_reason <- ifelse(
    !fixed_pass,
    paste0("fixed_QC:", failure_reasons),
    ifelse(tolower(sc_class_all[cells]) == "doublet", "scDblFinder_doublet", "included")
  )
  decision <- data.frame(
    cell_id_final = cells,
    original_barcode = initial$original_barcode[match(cells, colnames(initial))],
    sample_id = sample_id,
    nCount_RNA = as.numeric(ncount[cells]),
    nFeature_RNA = as.numeric(nfeature[cells]),
    mt_percent = as.numeric(mt_percent[cells]),
    HB_percent = as.numeric(hb_percent[cells]),
    pass_nFeature_RNA_min = pass_nfeature_low,
    pass_nFeature_RNA_max = pass_nfeature_high,
    pass_nCount_RNA_min = pass_ncount,
    pass_mt_percent_max = pass_mt,
    pass_HB_percent_max = pass_hb,
    fixed_qc_pass = fixed_pass,
    fixed_qc_pass_no_ncount_sensitivity = no_ncount_sensitivity_pass,
    fixed_qc_failure_reasons = failure_reasons,
    scDblFinder_score = as.numeric(sc_score_all[cells]),
    scDblFinder_class = as.character(sc_class_all[cells]),
    DoubletFinder_pANN = as.numeric(df_pann_all[cells]),
    DoubletFinder_class = as.character(df_class_all[cells]),
    final_main_include = final_include,
    final_exclusion_reason = final_reason,
    stringsAsFactors = FALSE
  )
  decision_rows[[sample_id]] <- decision
  plot_sample_qc(decision, sample_id, config$paths$qc_dir)

  globin_present <- config$qc$globin_panel[config$qc$globin_panel %in% toupper(rownames(working))]
  threshold_rows[[sample_id]] <- cbind(
    data.frame(
      sample_id = sample_id,
      working_feature_count = nrow(working),
      globin_panel_present = paste(globin_present, collapse = "|"),
      globin_panel_missing = paste(setdiff(config$qc$globin_panel, globin_present), collapse = "|"),
      rule = "nFeature>=500 & nFeature<6000 & nCount>1000 & mt<=20 & HB<5",
      no_ncount_sensitivity_rule = "nFeature>=500 & nFeature<6000 & mt<=20 & HB<5",
      stringsAsFactors = FALSE
    ),
    as.data.frame(c(
      f1_summary_stats(ncount, "nCount_RNA"),
      f1_summary_stats(nfeature, "nFeature_RNA"),
      f1_summary_stats(mt_percent, "mt_percent"),
      f1_summary_stats(hb_percent, "HB_percent")
    ), check.names = FALSE)
  )
  summary_rows[[sample_id]] <- data.frame(
    sample_id = sample_id,
    public_cells = length(cells),
    fail_nFeature_min = sum(!pass_nfeature_low),
    fail_nFeature_max = sum(!pass_nfeature_high),
    fail_nCount_min = sum(!pass_ncount),
    fail_mt_max = sum(!pass_mt),
    fail_HB_max = sum(!pass_hb),
    fixed_qc_union_excluded = sum(!fixed_pass),
    fixed_qc_passed = sum(fixed_pass),
    no_ncount_sensitivity_qc_passed = sum(no_ncount_sensitivity_pass, na.rm = TRUE),
    no_ncount_sensitivity_extra_vs_main =
      sum(no_ncount_sensitivity_pass & !fixed_pass, na.rm = TRUE),
    scDblFinder_doublets = sum(tolower(sc_class) == "doublet"),
    final_retained = length(retained_cells),
    final_retained_fraction = length(retained_cells) / length(cells),
    stringsAsFactors = FALSE
  )
  doublet_rows[[sample_id]] <- data.frame(
    sample_id = sample_id,
    input_cells = length(fixed_cells),
    scDblFinder_version = as.character(utils::packageVersion("scDblFinder")),
    scDblFinder_workers = config$doublet$scdblfinder_workers,
    scDblFinder_expected_rate_source = "automatic_10x_rate_dbr_per1k_0.008",
    scDblFinder_doublets = sum(tolower(sc_class) == "doublet"),
    DoubletFinder_version = as.character(utils::packageVersion("DoubletFinder")),
    DoubletFinder_status = df_result$status,
    DoubletFinder_pN = config$doublet$doubletfinder_pN,
    DoubletFinder_pK = df_result$pK,
    DoubletFinder_nExp = df_result$nExp,
    DoubletFinder_expected_rate = df_result$expected_rate,
    DoubletFinder_only_doublets = sum(
      tolower(df_result$class[names(sc_class)]) == "doublet" & tolower(sc_class) == "singlet",
      na.rm = TRUE
    ),
    DoubletFinder_error = df_result$error,
    final_deletion_rule = "scDblFinder_doublet_only",
    seed = config$seed,
    stringsAsFactors = FALSE
  )
  rm(raw, working, fixed_counts, algorithm_counts, retained_counts, sce, sc_meta)
  gc(verbose = FALSE)
}

filtered <- sample_objects[[1]]
if (length(sample_objects) > 1L) filtered <- merge(filtered, y = sample_objects[-1], merge.data = FALSE)
filtered <- f1_join_assay(filtered, "RNA")
f1_assert_integer_counts(SeuratObject::LayerData(filtered, assay = "RNA", layer = "counts"), "F1.2 filtered RNA counts")

decisions <- do.call(rbind, decision_rows)
rownames(decisions) <- decisions$cell_id_final
meta_to_add <- decisions[colnames(initial), setdiff(colnames(decisions), c("cell_id_final", "original_barcode", "sample_id")), drop = FALSE]
initial <- Seurat::AddMetaData(initial, metadata = meta_to_add)
raw_signature_after <- f1_counts_signature(SeuratObject::LayerData(initial, assay = "RNA", layer = "counts"))
if (!identical(unname(raw_signature_before), unname(raw_signature_after))) {
  stop("向01对象追加QC metadata后raw counts签名改变，已停止保存。")
}

f1_write_tsv(decisions, file.path(config$paths$qc_dir, "qc_cell_decision_audit.tsv"))
f1_write_tsv(do.call(rbind, threshold_rows), file.path(config$paths$qc_dir, "qc_thresholds_by_sample.tsv"))
f1_write_tsv(do.call(rbind, summary_rows), file.path(config$paths$qc_dir, "sample_qc_summary.tsv"))
f1_write_tsv(do.call(rbind, doublet_rows), file.path(config$paths$qc_dir, "doublet_summary_by_sample.tsv"))

limitations <- data.frame(
  assessment_item = c("barcode_rank_knee", "Cell_Ranger_cell_calling", "emptyDrops", "SoupX", "CellBender"),
  required_input = c("raw_droplet_matrix", "raw_or_filtered_feature_barcode_and_pipeline_outputs", "raw_droplet_matrix", "raw_droplet_matrix", "raw_droplet_matrix"),
  input_available = FALSE,
  status = "not_evaluable_input_limited",
  reason = "GSE183904_public_input_contains_called_or_retained_cell_count_matrices_only",
  permitted_substitute = c("library_size_rank_only", "none", "none", "DecontX_after_coarse_annotation", "none"),
  substitute_interpretation_limit = c(
    "not_a_true_barcode_knee", "cannot_reconstruct_cell_calling", "cannot_test_empty_droplets",
    "does_not_model_empty_droplet_background", "cannot_run_cell_calling_background_branch"
  ),
  stringsAsFactors = FALSE
)
f1_write_tsv(limitations, file.path(config$paths$qc_dir, "qc_input_limitation_audit.tsv"))
f1_save_rds_atomic(initial, config$paths$object_01, compress = FALSE)
f1_save_rds_atomic(filtered, config$paths$object_02, compress = FALSE)
f1_save_session_info(config, "F1_02_qc_doublet")
f1_append_log(
  config,
  stage,
  sprintf("完成：公开细胞%d，固定QC+scDblFinder后保留%d；raw counts保持主矩阵", ncol(initial), ncol(filtered))
)
message("F1.2完成：", config$paths$object_02)
