# F1.6 inferCNV、CopyKAT与恶性上皮判定 ------------------------------------
#
# 生物学目的：在可靠上皮细胞内，用大片段CNV模式、CopyKAT非整倍体结果、
# 上皮marker和组织背景联合区分恶性、非恶性与不确定细胞。
#
# 重要边界：inferCNV和CopyKAT都来自RNA表达，只能称互补方法稳健性，不能称
# 独立DNA验证。正式CNV计算需要额外参数 --approve-cnv-execution，防止误触。

file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), winslash = "/"))
source(file.path(script_dir, "F1_config.R"), encoding = "UTF-8")
source(file.path(script_dir, "F1_utils.R"), encoding = "UTF-8")

args <- f1_parse_args()
config <- f1_build_config(args$project_root)
stage <- "F1.6"
outputs <- c(
  config$paths$malignancy_review_template,
  config$paths$object_05,
  config$paths$object_06a,
  config$paths$object_06b,
  file.path(config$paths$malignancy_dir, "F1_CNV_reference_ambient_preflight.tsv"),
  file.path(config$paths$malignancy_dir, "infercnv_subcluster_membership.tsv"),
  file.path(config$paths$malignancy_dir, "infercnv_subcluster_summary.tsv"),
  file.path(config$paths$malignancy_dir, "copykat_cell_calls.tsv"),
  file.path(config$paths$malignancy_dir, "copykat_three_arm_group_comparison.tsv"),
  file.path(config$paths$malignancy_dir, "copykat_self_estimated_diploid_composition.tsv"),
  file.path(config$paths$malignancy_dir, "copykat_secretory_spike_test.tsv"),
  file.path(config$paths$malignancy_dir, "infercnv_ambient_sensitivity_comparison.tsv"),
  file.path(config$paths$malignancy_dir, "malignant_cell_calling_summary.tsv"),
  file.path(config$paths$malignancy_dir, "F1_final_report.md")
)

if (!args$execute) {
  f1_stage_dry_run(
    stage,
    paste0(
      "按样本运行raw-count inferCNV Leiden子聚类和当前样本内CopyKAT三臂，",
      "必要时运行corrected inferCNV敏感性并形成06a/06b对象"
    ),
    c(
      config$paths$object_04, config$paths$object_03,
      config$paths$epithelial_review_approved, config$paths$gene_order,
      config$paths$marker_panel, config$paths$ambient_summary
    ),
    outputs,
    config$packages[[stage]]
  )
  cat("CNV安全开关：正式计算还需 --approve-cnv-execution。\n")
  quit(save = "no", status = 0)
}

f1_require_packages(config$packages[[stage]], stage)
suppressPackageStartupMessages(
  library("copykat", character.only = TRUE)
)
copykat_required_data <- c("full.anno", "DNA.hg20", "cyclegenes")
copykat_package_environment <- as.environment("package:copykat")
copykat_missing_data <- copykat_required_data[
  !vapply(
    copykat_required_data,
    exists,
    logical(1),
    envir = copykat_package_environment,
    inherits = FALSE
  )
]
if (length(copykat_missing_data)) {
  stop(
    "CopyKAT包未正确挂载官方LazyData：",
    paste(copykat_missing_data, collapse = ", ")
  )
}
copykat_data_source <- paste0(
  "copykat_", as.character(utils::packageVersion("copykat")),
  "_LazyData_full.anno_DNA.hg20_cyclegenes"
)
f1_check_f0_ready(config)
if (!file.exists(config$paths$object_04)) stop("请先完成F1.5：", config$paths$object_04)
if (!file.exists(config$paths$object_03)) stop("缺少全细胞注释对象：", config$paths$object_03)
if (!file.exists(config$paths$epithelial_review_approved)) {
  stop("F1.6前必须审核上皮cluster并建立：", config$paths$epithelial_review_approved)
}
if (!file.exists(config$paths$gene_order)) stop("缺少inferCNV gene order：", config$paths$gene_order)
if (!file.exists(config$paths$ambient_summary)) {
  stop("缺少注释后DecontX逐样本摘要，请先完整完成F1.4：", config$paths$ambient_summary)
}
f1_prepare_directories(config)
set.seed(config$seed)
f1_append_log(config, stage, "开始准备上皮CNV证据与联合恶性判定")

epithelial <- readRDS(config$paths$object_04)
all_cells <- readRDS(config$paths$object_03)
f1_require_columns(
  epithelial[[]],
  c("sample_id", "patient_id", "group_analysis", "epithelial_cluster_id"),
  "F1.5 epithelial object metadata"
)
f1_require_columns(
  all_cells[[]],
  c(
    "sample_id", "patient_id", "group_analysis", "cell_type_major",
    "annotation_confidence", "retained_cell_ambient_contamination_estimate",
    "decontX_evaluation_status"
  ),
  "F1.4 all-cell object metadata"
)
epithelial_review <- f1_read_tsv(config$paths$epithelial_review_approved)
f1_require_columns(
  epithelial_review,
  c("epithelial_cluster", "include_in_malignancy", "contamination_status", "epithelial_subtype", "review_reason"),
  "F1_epithelial_cluster_review_approved.tsv"
)
clusters <- sort(unique(as.character(epithelial$epithelial_cluster_id)))
if (anyDuplicated(as.character(epithelial_review$epithelial_cluster)) ||
    !setequal(as.character(epithelial_review$epithelial_cluster), clusters)) {
  stop("批准的上皮cluster复核表必须与F1.5当前cluster一对一对应。")
}
epithelial_review$include_in_malignancy <- f1_as_logical(
  epithelial_review$include_in_malignancy,
  "include_in_malignancy"
)
approved_clusters <- as.character(epithelial_review$epithelial_cluster[epithelial_review$include_in_malignancy])
if (!length(approved_clusters)) stop("没有上皮cluster获准进入恶性判定。")
epithelial_review_index <- match(epithelial$epithelial_cluster_id, epithelial_review$epithelial_cluster)
epithelial$epithelial_subtype <- epithelial_review$epithelial_subtype[epithelial_review_index]
epithelial$epithelial_contamination_status <- epithelial_review$contamination_status[epithelial_review_index]
epithelial$epithelial_review_reason <- epithelial_review$review_reason[epithelial_review_index]
candidate_cells <- colnames(epithelial)[epithelial$epithelial_cluster_id %in% approved_clusters]

resource_summary <- data.frame(
  item = c(
    "candidate_epithelial_cells", "candidate_samples", "candidate_clusters",
    "requested_sample_workers", "requested_infercnv_threads_per_sample",
    "requested_copykat_cores_per_sample", "maximum_requested_nested_cores"
  ),
  value = c(
    length(candidate_cells), length(unique(epithelial$sample_id[match(candidate_cells, colnames(epithelial))])),
    length(approved_clusters), config$cnv$sample_workers,
    config$cnv$infercnv_threads, config$cnv$copykat_cores,
    config$cnv$sample_workers *
      max(config$cnv$infercnv_threads, config$cnv$copykat_cores)
  ),
  interpretation = c(
    "实际CNV观察细胞数", "逐样本运行单元", "人工批准进入CNV的上皮cluster数",
    "同时运行的样本数", "每个inferCNV任务的线程数",
    "每个CopyKAT任务的核数", "按样本并行和样本内并行估算的峰值核数"
  ),
  stringsAsFactors = FALSE
)
f1_write_tsv(resource_summary, file.path(config$paths$malignancy_dir, "F1_CNV_actual_input_resource_summary.tsv"))

# 新方法改变了observation分组和HMM设置，使用新目录；旧cnv_by_sample原样保留为历史结果。
cnv_dir <- file.path(
  config$paths$malignancy_dir,
  "cnv_by_sample_v2_single_observation_i6"
)
dir.create(cnv_dir, recursive = TRUE, showWarnings = FALSE)
reference_summary_path <- file.path(cnv_dir, "infercnv_reference_summary.tsv")
subcluster_membership_path <- file.path(cnv_dir, "infercnv_subcluster_membership_all.tsv")
reference_manifest_path <- file.path(cnv_dir, "cnv_reference_cell_manifest.tsv")
copykat_baseline_audit_path <- file.path(
  cnv_dir,
  "copykat_self_estimated_diploid_composition_all.tsv"
)
infercnv_reference_policy <- paste(
  "same_sample",
  "same_patient_normal_gastric",
  "balanced_other_patient_normal_gastric",
  sep = "_then_"
)
copykat_input_scope <- "all_QC_singlet_cells_from_current_sample_only"
existing_scores <- all(file.exists(c(
  file.path(cnv_dir, "infercnv_cell_scores_all.tsv"),
  reference_summary_path,
  subcluster_membership_path,
  reference_manifest_path
)))
existing_copykat <- all(file.exists(c(
  file.path(cnv_dir, "copykat_cell_calls_all.tsv"),
  copykat_baseline_audit_path
)))
if (existing_scores) {
  existing_reference <- f1_read_tsv(reference_summary_path)
  provenance_fields <- c(
    "analysis_contract_version", "infercnv_observation_group",
    "infercnv_analysis_mode", "infercnv_internal_subclustering",
    "infercnv_partition_method", "infercnv_k_nn", "infercnv_leiden_resolution",
    "infercnv_leiden_method", "infercnv_leiden_function",
    "infercnv_inspect_subclusters", "infercnv_HMM", "infercnv_HMM_type",
    "infercnv_HMM_report_by", "infercnv_BayesMaxPNormal",
    "infercnv_reassignCNVs", "infercnv_denoise",
    "infercnv_cutoff", "infercnv_scipen", "infercnv_bitmap_type",
    "reference_minimum_cells", "reference_maximum_cells",
    "infercnv_reference_policy", "copykat_input_scope",
    "copykat_external_reference", "copykat_internal_min_genes_per_cell",
    "copykat_data_source"
  )
  provenance_ok <- isTRUE(
    all(provenance_fields %in% colnames(existing_reference)) &&
      all(existing_reference$analysis_contract_version == "F1_06_v2_single_observation_i6") &&
      all(existing_reference$infercnv_observation_group == config$cnv$infercnv_observation_group) &&
      all(existing_reference$infercnv_analysis_mode == config$cnv$infercnv_analysis_mode) &&
      all(tolower(as.character(existing_reference$infercnv_internal_subclustering)) == "true") &&
      all(
        existing_reference$infercnv_partition_method ==
          config$cnv$infercnv_tumor_subcluster_partition_method
      ) &&
      all(as.integer(existing_reference$infercnv_k_nn) == config$cnv$infercnv_k_nn) &&
      all(
        as.character(existing_reference$infercnv_leiden_resolution) ==
          as.character(config$cnv$infercnv_leiden_resolution)
      ) &&
      all(existing_reference$infercnv_leiden_method == config$cnv$infercnv_leiden_method) &&
      all(existing_reference$infercnv_leiden_function == config$cnv$infercnv_leiden_function) &&
      all(tolower(as.character(existing_reference$infercnv_inspect_subclusters)) == "true") &&
      all(
        tolower(as.character(existing_reference$infercnv_HMM)) ==
          tolower(as.character(config$cnv$infercnv_hmm))
      ) &&
      all(existing_reference$infercnv_HMM_type == config$cnv$infercnv_hmm_type) &&
      all(existing_reference$infercnv_HMM_report_by == config$cnv$infercnv_hmm_report_by) &&
      all(
        as.numeric(existing_reference$infercnv_BayesMaxPNormal) ==
          config$cnv$infercnv_bayes_max_p_normal
      ) &&
      all(
        tolower(as.character(existing_reference$infercnv_reassignCNVs)) ==
          tolower(as.character(config$cnv$infercnv_reassign_cnvs))
      ) &&
      all(
        tolower(as.character(existing_reference$infercnv_denoise)) ==
          tolower(as.character(config$cnv$infercnv_denoise))
      ) &&
      all(as.numeric(existing_reference$infercnv_cutoff) == config$cnv$infercnv_cutoff) &&
      all(as.integer(existing_reference$infercnv_scipen) == config$cnv$infercnv_scipen) &&
      all(existing_reference$infercnv_bitmap_type == config$cnv$infercnv_bitmap_type) &&
      all(
        as.integer(existing_reference$reference_minimum_cells) ==
          config$cnv$minimum_reference_cells
      ) &&
      all(
        as.integer(existing_reference$reference_maximum_cells) ==
          config$cnv$maximum_reference_cells
      ) &&
      all(existing_reference$infercnv_reference_policy == infercnv_reference_policy) &&
      all(
        existing_reference$copykat_input_scope == copykat_input_scope
      ) &&
      all(tolower(as.character(existing_reference$copykat_external_reference)) == "prohibited") &&
      all(
        as.integer(existing_reference$copykat_internal_min_genes_per_cell) ==
          config$cnv$copykat_internal_min_genes_per_cell
      ) &&
      all(existing_reference$copykat_data_source == copykat_data_source)
  )
  if (!provenance_ok) {
    existing_scores <- FALSE
    existing_copykat <- FALSE
    f1_append_log(
      config,
      stage,
      paste0(
        "既有CNV汇总与当前单一observation组、正式i6 HMM或",
        "CopyKAT三臂规则不一致，禁止复用。"
      )
    )
  }
}
all_cells <- f1_join_assay(all_cells, "RNA")
full_counts <- SeuratObject::LayerData(all_cells, assay = "RNA", layer = "counts")
f1_assert_integer_counts(full_counts, "F1.6 full-cell RNA counts")
gene_order <- f1_read_tsv(config$paths$gene_order)
f1_require_columns(gene_order, c("gene_symbol", "chromosome", "start", "end"), "inferCNV gene order")
gene_order <- gene_order[grepl("^chr([1-9]|1[0-9]|2[0-2])$", gene_order$chromosome), , drop = FALSE]
gene_order_duplicate_count <- sum(duplicated(gene_order$gene_symbol))
gene_order <- gene_order[!duplicated(gene_order$gene_symbol), , drop = FALSE]
matched_gene_order <- gene_order[gene_order$gene_symbol %in% rownames(full_counts), , drop = FALSE]
if (!nrow(matched_gene_order)) stop("inferCNV gene order与RNA gene symbol没有交集。")
mapping_report <- data.frame(
  gene_order_source = config$paths$gene_order,
  source_chr1_22_genes = nrow(gene_order),
  RNA_features = nrow(full_counts),
  matched_genes = nrow(matched_gene_order),
  matched_fraction_of_gene_order = nrow(matched_gene_order) / nrow(gene_order),
  duplicate_gene_symbols_removed = gene_order_duplicate_count,
  chromosomes = paste(unique(matched_gene_order$chromosome), collapse = "|"),
  input_assay_layer = "RNA_counts",
  status = if (nrow(matched_gene_order) >= 1000) "PASS_BASIC_ID_COMPATIBILITY" else "FAIL_LIKELY_ID_MISMATCH",
  stringsAsFactors = FALSE
)
f1_write_tsv(mapping_report, file.path(config$paths$malignancy_dir, "gene_order_mapping_report.tsv"))
if (mapping_report$status == "FAIL_LIKELY_ID_MISMATCH") {
  stop("inferCNV gene-order匹配少于1000个基因，提示gene ID或输入明显不匹配。")
}

f1_valid_patient_id <- function(x) {
  x <- trimws(as.character(x))
  !is.na(x) & nzchar(x) &
    !tolower(x) %in% c("unknown", "na", "n/a", "none", "not_available")
}

f1_round_robin_cap <- function(cells, meta, maximum) {
  cells <- unique(intersect(as.character(cells), rownames(meta)))
  if (!length(cells) || length(cells) <= maximum) return(cells)
  source_sample <- as.character(meta$sample_id[match(cells, rownames(meta))])
  by_sample <- split(sort(cells), source_sample)
  by_sample <- by_sample[sort(names(by_sample))]
  position <- stats::setNames(rep(1L, length(by_sample)), names(by_sample))
  selected <- character()
  while (length(selected) < maximum) {
    added <- FALSE
    for (source_id in names(by_sample)) {
      index <- position[[source_id]]
      if (index <= length(by_sample[[source_id]])) {
        selected <- c(selected, by_sample[[source_id]][[index]])
        position[[source_id]] <- index + 1L
        added <- TRUE
        if (length(selected) >= maximum) break
      }
    }
    if (!added) break
  }
  selected
}

f1_reference_pool <- function(meta, mask, minimum, maximum) {
  mask <- !is.na(mask) & mask
  confidence_ok <- tolower(as.character(meta$annotation_confidence)) %in% c("high", "medium")
  tnk <- rownames(meta)[mask & confidence_ok & meta$cell_type_major %in% "T/NK"]
  bplasma <- rownames(meta)[mask & confidence_ok & meta$cell_type_major %in% "B/Plasma"]
  use_bplasma <- length(tnk) < minimum
  candidates <- if (use_bplasma) unique(c(tnk, bplasma)) else tnk
  selected <- f1_round_robin_cap(candidates, meta, maximum)
  list(
    cells = selected,
    sufficient = length(selected) >= minimum,
    lineage_rule = if (use_bplasma) "T_NK_plus_B_Plasma" else "T_NK",
    available_before_cap = length(candidates)
  )
}

choose_reference_cells <- function(sample_id, all_cells, minimum, maximum) {
  meta <- all_cells[[]]
  target_rows <- as.character(meta$sample_id) %in% sample_id
  if (!any(target_rows)) stop(sample_id, "在全细胞对象中不存在。")
  target_patient_values <- unique(as.character(meta$patient_id[target_rows]))
  if (length(target_patient_values) != 1L) {
    stop(sample_id, "对应多个patient_id，不能可靠选择inferCNV reference。")
  }
  target_patient <- target_patient_values[[1]]
  target_patient_valid <- f1_valid_patient_id(target_patient)
  normal_gastric <- as.character(meta$group_analysis) %in% "Normal_Gastric"
  same_sample <- as.character(meta$sample_id) %in% sample_id
  same_patient_normal <- rep(FALSE, nrow(meta))
  if (target_patient_valid) {
    same_patient_normal <- normal_gastric &
      as.character(meta$patient_id) %in% target_patient &
      !same_sample
  }
  source_patient_valid <- f1_valid_patient_id(meta$patient_id)
  other_patient_normal <- normal_gastric & !same_sample & source_patient_valid
  if (target_patient_valid) {
    other_patient_normal <- other_patient_normal &
      !as.character(meta$patient_id) %in% target_patient
  }

  tier_masks <- list(
    same_sample = same_sample,
    same_patient_normal_gastric = same_patient_normal,
    balanced_other_patient_normal_gastric = other_patient_normal
  )
  attempts <- lapply(
    tier_masks,
    function(mask) f1_reference_pool(meta, mask, minimum, maximum)
  )
  selected_tier <- names(attempts)[vapply(attempts, function(x) x$sufficient, logical(1))]
  if (length(selected_tier)) {
    selected_tier <- selected_tier[[1]]
    selected <- attempts[[selected_tier]]
    status <- "evaluable"
  } else {
    available <- vapply(attempts, function(x) length(x$cells), integer(1))
    selected_tier <- names(which.max(available))[[1]]
    selected <- attempts[[selected_tier]]
    status <- "not_evaluable_all_reference_tiers_below_minimum"
  }
  selected_cells <- selected$cells
  source_samples <- sort(unique(as.character(meta$sample_id[match(selected_cells, rownames(meta))])))
  source_patients <- sort(unique(as.character(meta$patient_id[match(selected_cells, rownames(meta))])))
  list(
    cells = selected_cells,
    source = paste(selected_tier, selected$lineage_rule, sep = "__"),
    status = status,
    target_patient_id = target_patient,
    source_samples = source_samples,
    source_patients = source_patients,
    available_before_cap = selected$available_before_cap,
    tier_counts = vapply(attempts, function(x) length(x$cells), integer(1))
  )
}

choose_copykat_known_normals <- function(sample_id, all_cells, minimum, maximum) {
  meta <- all_cells[[]]
  pool <- f1_reference_pool(
    meta,
    as.character(meta$sample_id) %in% sample_id,
    minimum,
    maximum
  )
  if (!pool$sufficient) {
    return(list(
      cells = character(),
      candidate_cells = pool$cells,
      source = "automatic_within_current_sample",
      available_same_sample_candidates = pool$available_before_cap
    ))
  }
  list(
    cells = pool$cells,
    candidate_cells = pool$cells,
    source = paste("known_same_sample", pool$lineage_rule, sep = "__"),
    available_same_sample_candidates = pool$available_before_cap
  )
}

f1_reference_manifest <- function(
    target_sample_id,
    method,
    cells,
    source,
    all_cells
) {
  meta <- all_cells[[]]
  cells <- intersect(cells, rownames(meta))
  if (!length(cells)) {
    return(data.frame(
      target_sample_id = character(),
      method = character(),
      reference_cell_id = character(),
      reference_sample_id = character(),
      reference_patient_id = character(),
      reference_group_analysis = character(),
      reference_cell_type_major = character(),
      annotation_confidence = character(),
      reference_source = character(),
      stringsAsFactors = FALSE
    ))
  }
  index <- match(cells, rownames(meta))
  data.frame(
    target_sample_id = target_sample_id,
    method = method,
    reference_cell_id = cells,
    reference_sample_id = as.character(meta$sample_id[index]),
    reference_patient_id = as.character(meta$patient_id[index]),
    reference_group_analysis = as.character(meta$group_analysis[index]),
    reference_cell_type_major = as.character(meta$cell_type_major[index]),
    annotation_confidence = as.character(meta$annotation_confidence[index]),
    reference_source = source,
    stringsAsFactors = FALSE
  )
}

extract_infercnv_subcluster_membership <- function(
    infer_object,
    obs_cells,
    obs_cluster,
    sample_id
) {
  nested <- infer_object@tumor_subclusters$subclusters
  observation_groups <- names(infer_object@observation_grouped_cell_indices)
  if (is.null(nested) || !length(nested) || !length(observation_groups)) {
    stop(sample_id, " inferCNV结果缺少observation subcluster结构。")
  }
  rows <- list()
  row_index <- 0L
  for (group_name in observation_groups) {
    group_subclusters <- nested[[group_name]]
    if (is.null(group_subclusters) || !length(group_subclusters)) {
      stop(sample_id, " 的inferCNV observation group缺少subcluster：", group_name)
    }
    local_names <- names(group_subclusters)
    if (is.null(local_names) || any(!nzchar(local_names))) {
      local_names <- paste0(group_name, "_s", seq_along(group_subclusters))
    }
    for (i in seq_along(group_subclusters)) {
      indices <- group_subclusters[[i]]
      cells <- names(indices)
      if (is.null(cells) || any(!nzchar(cells))) {
        cells <- colnames(infer_object@expr.data)[as.integer(indices)]
      }
      cells <- intersect(as.character(cells), obs_cells)
      if (!length(cells)) next
      row_index <- row_index + 1L
      rows[[row_index]] <- data.frame(
        cell_id_final = cells,
        sample_id = sample_id,
        # 原表达cluster只是细胞注释；inferCNV可以按CNV模式把多个表达cluster的细胞归到同一组。
        epithelial_cluster = unname(obs_cluster[cells]),
        infercnv_parent_group = group_name,
        infercnv_subcluster_local = local_names[[i]],
        infercnv_subcluster_id = paste(sample_id, local_names[[i]], sep = "__"),
        infercnv_evaluation_status = "evaluable",
        stringsAsFactors = FALSE
      )
    }
  }
  if (!length(rows)) stop(sample_id, "未提取到任何observation subcluster成员。")
  membership <- do.call(rbind, rows)
  if (anyDuplicated(membership$cell_id_final)) {
    stop(sample_id, " inferCNV subcluster成员出现重复cell ID。")
  }
  if (!setequal(membership$cell_id_final, obs_cells)) {
    missing <- setdiff(obs_cells, membership$cell_id_final)
    extra <- setdiff(membership$cell_id_final, obs_cells)
    stop(
      sample_id, " inferCNV subcluster未完整覆盖观察细胞；missing=",
      length(missing), "，extra=", length(extra)
    )
  }
  membership[match(obs_cells, membership$cell_id_final), , drop = FALSE]
}

infercnv_cell_scores <- function(
    infer_object,
    reference_cells,
    sample_id,
    epithelial_cluster,
    subcluster_membership
) {
  expr <- infer_object@expr.data
  reference_cells <- intersect(reference_cells, colnames(expr))
  if (!length(reference_cells)) stop(sample_id, " inferCNV结果中没有reference cell。")
  ref_center <- Matrix::rowMeans(expr[, reference_cells, drop = FALSE])
  scores <- numeric(ncol(expr))
  chunk_starts <- seq.int(1L, ncol(expr), by = 250L)
  for (start in chunk_starts) {
    index <- start:min(start + 249L, ncol(expr))
    block <- as.matrix(expr[, index, drop = FALSE])
    scores[index] <- colMeans(abs(sweep(block, 1, ref_center, FUN = "-")))
  }
  names(scores) <- colnames(expr)
  threshold <- as.numeric(stats::quantile(scores[reference_cells], 0.95, na.rm = TRUE, type = 7))
  cells <- intersect(names(epithelial_cluster), names(scores))
  data.frame(
    cell_id_final = cells,
    sample_id = sample_id,
    epithelial_cluster = unname(epithelial_cluster[cells]),
    infercnv_parent_group = subcluster_membership$infercnv_parent_group[
      match(cells, subcluster_membership$cell_id_final)
    ],
    infercnv_subcluster_local = subcluster_membership$infercnv_subcluster_local[
      match(cells, subcluster_membership$cell_id_final)
    ],
    infercnv_subcluster_id = subcluster_membership$infercnv_subcluster_id[
      match(cells, subcluster_membership$cell_id_final)
    ],
    infercnv_cell_burden = unname(scores[cells]),
    reference_background_P95 = threshold,
    above_reference_P95 = unname(scores[cells]) > threshold,
    stringsAsFactors = FALSE
  )
}

run_infercnv_matrix <- function(
    counts,
    input_cells,
    obs_cells,
    obs_cluster,
    reference_cells,
    sample_id,
    gene_order,
    out_dir,
    file_prefix,
    config
) {
  if (
    anyDuplicated(input_cells) ||
      !setequal(input_cells, c(obs_cells, reference_cells)) ||
      length(intersect(obs_cells, reference_cells))
  ) {
    stop(sample_id, " inferCNV观察细胞与reference cell集合不合法。")
  }
  if (!setequal(colnames(counts), input_cells)) {
    stop(sample_id, " inferCNV count矩阵未完整覆盖指定输入细胞。")
  }
  if (!setequal(rownames(counts), gene_order$gene_symbol)) {
    stop(sample_id, " inferCNV count矩阵与gene order基因集合不一致。")
  }
  counts <- counts[gene_order$gene_symbol, input_cells, drop = FALSE]
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  gene_order_file <- file.path(out_dir, paste0(file_prefix, "_gene_order.tsv"))
  data.table::fwrite(
    gene_order[, c("gene_symbol", "chromosome", "start", "end")],
    gene_order_file,
    sep = "\t",
    quote = FALSE,
    col.names = FALSE
  )
  annotation_file <- file.path(out_dir, paste0(file_prefix, "_annotations.tsv"))
  annotation <- data.frame(
    cell = input_cells,
    # 每个样本只设一个观察组，避免原表达cluster人为限制inferCNV内部子聚类。
    group = c(
      rep(config$cnv$infercnv_observation_group, length(obs_cells)),
      rep("reference", length(reference_cells))
    ),
    stringsAsFactors = FALSE
  )
  data.table::fwrite(annotation, annotation_file, sep = "\t", quote = FALSE, col.names = FALSE)

  old_scipen <- getOption("scipen")
  old_bitmap_type <- getOption("bitmapType")
  if (
    identical(config$cnv$infercnv_bitmap_type, "cairo") &&
      !isTRUE(capabilities("cairo"))
  ) {
    stop("当前R环境不支持cairo，不能在无X11服务器上绘制inferCNV子聚类PNG。")
  }
  options(scipen = max(config$cnv$infercnv_scipen, old_scipen))
  options(bitmapType = config$cnv$infercnv_bitmap_type)
  on.exit(options(scipen = old_scipen), add = TRUE)
  on.exit(options(bitmapType = old_bitmap_type), add = TRUE)
  infer_object <- infercnv::CreateInfercnvObject(
    raw_counts_matrix = counts,
    annotations_file = annotation_file,
    delim = "\t",
    gene_order_file = gene_order_file,
    ref_group_names = "reference"
  )
  set.seed(config$seed)
  infer_result <- infercnv::run(
    infer_object,
    cutoff = config$cnv$infercnv_cutoff,
    min_cells_per_gene = config$cnv$infercnv_min_cells_per_gene,
    window_length = config$cnv$infercnv_window_length,
    out_dir = out_dir,
    cluster_by_groups = TRUE,
    cluster_references = TRUE,
    # Leiden只帮助把热图中的相似CNV模式归组，不把自动分组解释为真实肿瘤亚克隆。
    analysis_mode = config$cnv$infercnv_analysis_mode,
    tumor_subcluster_partition_method =
      config$cnv$infercnv_tumor_subcluster_partition_method,
    k_nn = config$cnv$infercnv_k_nn,
    leiden_resolution = config$cnv$infercnv_leiden_resolution,
    leiden_method = config$cnv$infercnv_leiden_method,
    leiden_function = config$cnv$infercnv_leiden_function,
    inspect_subclusters = config$cnv$infercnv_inspect_subclusters,
    denoise = config$cnv$infercnv_denoise,
    HMM = config$cnv$infercnv_hmm,
    HMM_type = config$cnv$infercnv_hmm_type,
    HMM_report_by = config$cnv$infercnv_hmm_report_by,
    BayesMaxPNormal = config$cnv$infercnv_bayes_max_p_normal,
    reassignCNVs = config$cnv$infercnv_reassign_cnvs,
    num_threads = config$cnv$infercnv_threads,
    no_plot = FALSE,
    output_format = config$cnv$infercnv_output_format,
    useRaster = TRUE,
    write_expr_matrix = config$cnv$infercnv_write_expr_matrix,
    save_rds = TRUE,
    save_final_rds = TRUE,
    resume_mode = config$cnv$infercnv_resume_mode
  )
  final_object_path <- file.path(out_dir, paste0(file_prefix, "_final_object.rds"))
  saveRDS(infer_result, final_object_path, compress = FALSE)
  subcluster_membership <- extract_infercnv_subcluster_membership(
    infer_result,
    obs_cells,
    obs_cluster,
    sample_id
  )
  subcluster_membership_path <- file.path(
    out_dir,
    paste0(file_prefix, "_subcluster_membership.tsv")
  )
  data.table::fwrite(
    subcluster_membership,
    subcluster_membership_path,
    sep = "\t",
    quote = FALSE,
    na = "NA"
  )

  # 最终RDS包含重绘热图所需的expr.data；另存最终细胞顺序和输入分组，避免以后
  # 只能依赖默认图片。默认诊断图仍保留在完整out_dir中。
  final_cell_order <- colnames(infer_result@expr.data)
  membership_index <- match(final_cell_order, subcluster_membership$cell_id_final)
  plot_order <- data.frame(
    plot_order = seq_along(final_cell_order),
    cell_id = final_cell_order,
    annotation_group = annotation$group[match(final_cell_order, annotation$cell)],
    cell_role = ifelse(final_cell_order %in% reference_cells, "reference", "observation"),
    sample_id = sample_id,
    epithelial_cluster = subcluster_membership$epithelial_cluster[membership_index],
    infercnv_subcluster_id = subcluster_membership$infercnv_subcluster_id[membership_index],
    stringsAsFactors = FALSE
  )
  plot_order_path <- file.path(out_dir, paste0(file_prefix, "_final_plot_cell_order.tsv"))
  data.table::fwrite(plot_order, plot_order_path, sep = "\t", quote = FALSE, na = "NA")
  plot_source_manifest <- data.frame(
    artifact = c(
      "final_infercnv_object", "input_annotations", "input_gene_order",
      "final_plot_cell_order", "subcluster_membership", "infercnv_output_directory"
    ),
    path = c(
      final_object_path, annotation_file, gene_order_file,
      plot_order_path, subcluster_membership_path, out_dir
    ),
    purpose = c(
      "contains_final_expr_data_for_publication_redraw",
      "cell_groups_used_for_inference",
      "chromosome_order_used_for_inference",
      "exact_cell_order_for_heatmap_redraw",
      "observation_cell_to_infercnv_subcluster_mapping",
      "default_diagnostic_plots_and_intermediate_outputs"
    ),
    analysis_mode = config$cnv$infercnv_analysis_mode,
    internal_subclustering = identical(config$cnv$infercnv_analysis_mode, "subclusters"),
    observation_group = config$cnv$infercnv_observation_group,
    HMM = config$cnv$infercnv_hmm,
    HMM_type = config$cnv$infercnv_hmm_type,
    HMM_report_by = config$cnv$infercnv_hmm_report_by,
    stringsAsFactors = FALSE
  )
  data.table::fwrite(
    plot_source_manifest,
    file.path(out_dir, paste0(file_prefix, "_plot_source_manifest.tsv")),
    sep = "\t",
    quote = FALSE,
    na = "NA"
  )
  scores <- infercnv_cell_scores(
    infer_result,
    reference_cells,
    sample_id,
    obs_cluster,
    subcluster_membership
  )
  rm(infer_object, infer_result)
  scores
}

load_corrected_infercnv_counts <- function(input_cells, all_cells, genes, ambient_summary) {
  source_sample <- as.character(all_cells$sample_id[match(input_cells, colnames(all_cells))])
  if (anyNA(source_sample)) stop("corrected inferCNV输入细胞无法完整映射sample_id。")
  blocks <- lapply(unique(source_sample), function(sample_id) {
    row <- ambient_summary[ambient_summary$sample_id == sample_id, , drop = FALSE]
    if (nrow(row) != 1L) stop(sample_id, "在ambient摘要中不是唯一一行。")
    path <- as.character(row$corrected_counts_path[[1]])
    if (is.na(path) || !nzchar(path) || !file.exists(path)) {
      stop(sample_id, "缺少DecontX corrected counts，不能运行corrected inferCNV。")
    }
    cells <- input_cells[source_sample == sample_id]
    corrected <- readRDS(path)
    if (!all(cells %in% colnames(corrected))) {
      stop(sample_id, "的corrected矩阵未同时覆盖所需观察/reference细胞。")
    }
    present <- intersect(genes, rownames(corrected))
    corrected <- methods::as(corrected[present, cells, drop = FALSE], "dgCMatrix")
    entries <- Matrix::summary(corrected)
    Matrix::sparseMatrix(
      i = match(present[entries$i], genes),
      j = entries$j,
      x = entries$x,
      dims = c(length(genes), length(cells)),
      dimnames = list(genes, cells)
    )
  })
  combined <- do.call(cbind, blocks)
  combined[, input_cells, drop = FALSE]
}

empty_copykat_calls <- function(cells, sample_id, arm_id, note) {
  data.frame(
    cell_id_final = as.character(cells),
    sample_id = rep(sample_id, length(cells)),
    copykat_arm = rep(arm_id, length(cells)),
    copykat_call = rep("uncalled", length(cells)),
    copykat_returned = rep(FALSE, length(cells)),
    copykat_note = rep(note, length(cells)),
    stringsAsFactors = FALSE
  )
}

extract_copykat_calls <- function(result, candidate_cells, sample_id, arm_id) {
  if (!length(candidate_cells)) {
    return(empty_copykat_calls(candidate_cells, sample_id, arm_id, "no_evaluation_cells"))
  }
  prediction <- result$prediction
  if (is.null(prediction)) stop(sample_id, " CopyKAT结果缺少prediction。")
  prediction <- as.data.frame(prediction, stringsAsFactors = FALSE)
  f1_require_columns(prediction, c("cell.names", "copykat.pred"), paste0(sample_id, " CopyKAT prediction"))
  raw_call <- setNames(tolower(as.character(prediction$copykat.pred)), as.character(prediction$cell.names))
  call <- setNames(rep("uncalled", length(candidate_cells)), candidate_cells)
  matched <- intersect(candidate_cells, names(raw_call))
  normalized <- ifelse(
    grepl("aneuploid", raw_call[matched]), "aneuploid",
    ifelse(grepl("diploid", raw_call[matched]), "diploid", "uncalled")
  )
  call[matched] <- normalized
  data.frame(
    cell_id_final = candidate_cells,
    sample_id = sample_id,
    copykat_arm = arm_id,
    copykat_call = unname(call[candidate_cells]),
    copykat_returned = candidate_cells %in% names(raw_call),
    copykat_note = ifelse(candidate_cells %in% names(raw_call), "returned_by_copykat", "filtered_or_uncalled_by_copykat"),
    stringsAsFactors = FALSE
  )
}

run_copykat_arm <- function(
    copykat_input,
    sample_cells,
    evaluation_cells,
    sample_id,
    arm_id,
    arm_dir,
    known_normal_cells = character(),
    seed,
    config
) {
  sample_cells <- as.character(sample_cells)
  evaluation_cells <- intersect(as.character(evaluation_cells), sample_cells)
  known_normal_cells <- intersect(as.character(known_normal_cells), sample_cells)
  if (length(intersect(evaluation_cells, known_normal_cells))) {
    stop(sample_id, " ", arm_id, "的评价细胞不能同时作为known normal。")
  }
  dir.create(arm_dir, recursive = TRUE, showWarnings = FALSE)
  sam_name <- paste(sample_id, arm_id, sep = "__")
  prediction_path <- file.path(arm_dir, paste0(sam_name, "_copykat_prediction.txt"))
  clustering_path <- file.path(arm_dir, paste0(sam_name, "_copykat_clustering_results.rds"))
  raw_cna_path <- file.path(
    arm_dir,
    paste0(sam_name, "_copykat_CNA_raw_results_gene_by_cell.txt")
  )
  recovery_path <- file.path(arm_dir, paste0(sam_name, "_postprediction_recovery.tsv"))

  # CopyKAT 1.1.0偶尔会在prediction已写完后，于CNA输出阶段发生行数错误。
  # 只有prediction完整覆盖本次输入时才恢复分类；缺失的CNA矩阵和热图会明确降级记录。
  recover_prediction <- function(error_message) {
    if (!grepl("^arguments imply differing number of rows:", as.character(error_message))) {
      return(NULL)
    }
    recovery_inputs <- c(prediction_path, clustering_path, raw_cna_path)
    if (!all(file.exists(recovery_inputs)) || any(file.info(recovery_inputs)$size <= 0)) {
      return(NULL)
    }
    prediction <- f1_read_tsv(prediction_path)
    f1_require_columns(
      prediction,
      c("cell.names", "copykat.pred"),
      paste0(sample_id, " ", arm_id, " CopyKAT saved prediction")
    )
    prediction_ids <- as.character(prediction$cell.names)
    exact_cell_match <- length(prediction_ids) == length(sample_cells) &&
      !anyDuplicated(prediction_ids) && setequal(prediction_ids, sample_cells)
    if (!exact_cell_match || anyNA(prediction$copykat.pred)) return(NULL)
    recovery <- data.frame(
      sample_id = sample_id,
      copykat_arm = arm_id,
      recovery_status = "prediction_recovered_CNA_and_default_heatmap_unavailable",
      original_error = as.character(error_message),
      copykat_input_cells = length(sample_cells),
      prediction_rows = nrow(prediction),
      exact_cell_id_match = exact_cell_match,
      prediction_path = prediction_path,
      clustering_path = clustering_path,
      raw_cna_path = raw_cna_path,
      limitation = paste0(
        "官方prediction已保存；最终CNA_results和默认热图因包内后处理维度错误不可评估。"
      ),
      stringsAsFactors = FALSE
    )
    f1_write_tsv(recovery, recovery_path)
    list(
      result = list(prediction = prediction, recovery = recovery),
      status = "prediction_recovered_CNA_and_default_heatmap_unavailable",
      recovery_artifacts = c(recovery_inputs, recovery_path)
    )
  }

  old_wd <- getwd()
  old_bitmap_type <- getOption("bitmapType")
  if (
    identical(config$cnv$infercnv_bitmap_type, "cairo") &&
      !isTRUE(capabilities("cairo"))
  ) {
    stop("当前R环境不支持cairo，不能在无X11服务器上绘制CopyKAT热图。")
  }
  options(bitmapType = config$cnv$infercnv_bitmap_type)
  setwd(arm_dir)
  copykat_run <- tryCatch({
    set.seed(seed)
    copykat_args <- list(
      rawmat = copykat_input,
      id.type = "S",
      cell.line = "no",
      ngene.chr = config$cnv$copykat_ngene_chr,
      LOW.DR = 0.05,
      UP.DR = 0.10,
      win.size = config$cnv$copykat_win_size,
      KS.cut = config$cnv$copykat_ks_cut,
      sam.name = sam_name,
      distance = "euclidean",
      output.seg = "FALSE",
      plot.genes = "TRUE",
      genome = "hg20",
      n.cores = config$cnv$copykat_cores
    )
    # A臂必须完全省略norm.cell.names；B/C臂才显式传入同样本known normal。
    if (length(known_normal_cells)) {
      copykat_args$norm.cell.names <- known_normal_cells
    }
    unsupported <- setdiff(names(copykat_args), names(formals(copykat::copykat)))
    if (length(unsupported)) {
      stop("当前CopyKAT版本不支持参数：", paste(unsupported, collapse = ", "))
    }
    list(
      result = do.call(copykat::copykat, copykat_args),
      status = "complete",
      recovery_artifacts = character()
    )
  }, error = function(e) {
    recovered <- recover_prediction(conditionMessage(e))
    if (is.null(recovered)) stop(e)
    recovered
  }, finally = {
    setwd(old_wd)
    options(bitmapType = old_bitmap_type)
  })

  result_path <- file.path(arm_dir, paste0(sam_name, "_result.rds"))
  saveRDS(copykat_run$result, result_path, compress = FALSE)
  all_calls <- extract_copykat_calls(
    copykat_run$result, sample_cells, sample_id, arm_id
  )
  evaluation_calls <- all_calls[
    match(evaluation_cells, all_calls$cell_id_final),
    ,
    drop = FALSE
  ]
  list(
    all_calls = all_calls,
    evaluation_calls = evaluation_calls,
    status = copykat_run$status,
    known_normal_cells = known_normal_cells,
    required_artifacts = c(result_path, copykat_run$recovery_artifacts)
  )
}

audit_copykat_self_estimated_baseline <- function(
    all_calls,
    candidate_cells,
    sample_id,
    run_status,
    suspect_fraction
) {
  diploid_cells <- all_calls$cell_id_final[all_calls$copykat_call == "diploid"]
  epithelial_fraction <- if (length(diploid_cells)) {
    mean(diploid_cells %in% candidate_cells)
  } else {
    NA_real_
  }
  baseline_status <- if (!identical(run_status, "complete")) {
    "not_evaluable_CNA_and_default_heatmap_unavailable"
  } else if (!length(diploid_cells)) {
    "not_evaluable_no_predicted_diploid_cells"
  } else if (epithelial_fraction > suspect_fraction) {
    "not_evaluable_baseline_suspect"
  } else {
    "evaluable"
  }
  data.frame(
    sample_id = sample_id,
    copykat_arm = "A_self_estimated",
    copykat_run_status = run_status,
    predicted_diploid_cells = length(diploid_cells),
    predicted_diploid_candidate_epithelial_cells =
      sum(diploid_cells %in% candidate_cells),
    predicted_diploid_candidate_epithelial_fraction = epithelial_fraction,
    suspect_fraction_rule = suspect_fraction,
    copykat_primary_baseline_status = baseline_status,
    copykat_primary_evidence_usable = identical(baseline_status, "evaluable"),
    stringsAsFactors = FALSE
  )
}

f1_split_marker_field <- function(x) {
  x <- as.character(x %||% "")
  if (is.na(x) || !nzchar(x) || identical(x, "NA")) return(character())
  unique(trimws(strsplit(x, ",", fixed = TRUE)[[1]]))
}

f1_lineage_markers <- function(panel_row) {
  unique(c(
    f1_split_marker_field(panel_row$positive_markers[[1]]),
    f1_split_marker_field(panel_row$supporting_markers[[1]])
  ))
}

f1_marker_detection_summary <- function(counts, cells, genes) {
  cells <- intersect(cells, colnames(counts))
  present <- intersect(genes, rownames(counts))
  if (!length(cells) || !length(present)) return("not_evaluable")
  detection <- Matrix::rowMeans(counts[present, cells, drop = FALSE] > 0) * 100
  paste0(present, "=", sprintf("%.1f%%", detection), collapse = "|")
}

summarize_selected_ambient <- function(cells) {
  values <- as.numeric(
    all_cells$retained_cell_ambient_contamination_estimate[
      match(cells, colnames(all_cells))
    ]
  )
  finite <- values[is.finite(values)]
  c(
    evaluable_cells = length(finite),
    contamination_median = if (length(finite)) stats::median(finite) else NA_real_,
    contamination_P90 = if (length(finite)) {
      as.numeric(stats::quantile(finite, 0.90, names = FALSE, type = 7))
    } else {
      NA_real_
    }
  )
}

sample_ids <- sort(unique(as.character(
  epithelial$sample_id[match(candidate_cells, colnames(epithelial))]
)))
preflight_rows <- lapply(sample_ids, function(sample_id) {
  sample_cells <- colnames(all_cells)[as.character(all_cells$sample_id) %in% sample_id]
  obs_cells <- candidate_cells[
    epithelial$sample_id[match(candidate_cells, colnames(epithelial))] %in% sample_id
  ]
  reference <- choose_reference_cells(
    sample_id,
    all_cells,
    config$cnv$minimum_reference_cells,
    config$cnv$maximum_reference_cells
  )
  reference_cells <- setdiff(intersect(reference$cells, colnames(full_counts)), obs_cells)
  copykat_reference <- choose_copykat_known_normals(
    sample_id,
    all_cells,
    config$cnv$minimum_reference_cells,
    config$cnv$maximum_reference_cells
  )
  observation_ambient <- summarize_selected_ambient(obs_cells)
  infercnv_reference_ambient <- summarize_selected_ambient(reference_cells)
  copykat_reference_ambient <- summarize_selected_ambient(copykat_reference$cells)
  data.frame(
    sample_id = sample_id,
    observation_cells = length(obs_cells),
    copykat_input_cells = length(sample_cells),
    infercnv_reference_status = reference$status,
    infercnv_reference_source = reference$source,
    infercnv_reference_cells = length(reference_cells),
    copykat_A_baseline_mode = "self_estimated_without_norm.cell.names",
    copykat_B_reference_status = if (length(copykat_reference$cells)) {
      "evaluable"
    } else {
      "not_evaluable_insufficient_same_sample_immune_reference"
    },
    copykat_B_reference_source = copykat_reference$source,
    copykat_known_normal_cells = length(copykat_reference$cells),
    copykat_C_applicability = if (
      unique(as.character(all_cells$group_analysis[match(sample_cells, colnames(all_cells))])) %in%
        "Normal_Gastric"
    ) {
      "applicable_two_fold_normal_epithelial_holdout"
    } else {
      "not_applicable_non_normal_gastric"
    },
    observation_ambient_evaluable_cells = observation_ambient[["evaluable_cells"]],
    observation_contamination_median = observation_ambient[["contamination_median"]],
    observation_contamination_P90 = observation_ambient[["contamination_P90"]],
    infercnv_reference_ambient_evaluable_cells =
      infercnv_reference_ambient[["evaluable_cells"]],
    infercnv_reference_contamination_median =
      infercnv_reference_ambient[["contamination_median"]],
    infercnv_reference_contamination_P90 =
      infercnv_reference_ambient[["contamination_P90"]],
    copykat_reference_ambient_evaluable_cells =
      copykat_reference_ambient[["evaluable_cells"]],
    copykat_reference_contamination_median =
      copykat_reference_ambient[["contamination_median"]],
    copykat_reference_contamination_P90 =
      copykat_reference_ambient[["contamination_P90"]],
    stringsAsFactors = FALSE
  )
})
preflight <- do.call(rbind, preflight_rows)
preflight_path <- file.path(
  config$paths$malignancy_dir,
  "F1_CNV_reference_ambient_preflight.tsv"
)
f1_write_tsv(preflight, preflight_path)
if (
  config$cnv$require_explicit_execution_approval &&
    !(existing_scores && existing_copykat) &&
    !args$approve_cnv_execution
) {
  stop(
    "已完成F1.6输入、gene-order、reference与ambient预检，但尚未运行高资源CNV步骤。",
    "批准设备后重新运行并增加： --approve-cnv-execution；预检表：", preflight_path
  )
}

if (existing_scores && existing_copykat) {
  infer_scores_all <- f1_read_tsv(file.path(cnv_dir, "infercnv_cell_scores_all.tsv"))
  copykat_calls_all <- f1_read_tsv(file.path(cnv_dir, "copykat_cell_calls_all.tsv"))
  reference_summary_all <- f1_read_tsv(file.path(cnv_dir, "infercnv_reference_summary.tsv"))
  subcluster_membership_all <- f1_read_tsv(subcluster_membership_path)
  reference_manifest_all <- f1_read_tsv(reference_manifest_path)
  copykat_baseline_audit_all <- f1_read_tsv(copykat_baseline_audit_path)
  f1_append_log(config, stage, "复用已保存的逐样本inferCNV与CopyKAT结果")
} else {
  f1_character_md5 <- function(x) {
    path <- tempfile("f1_cnv_signature_")
    on.exit(unlink(path), add = TRUE)
    writeLines(enc2utf8(as.character(x)), path, useBytes = TRUE)
    toupper(unname(tools::md5sum(path)))
  }

  run_cnv_sample <- function(sample_id) {
    message("[", stage, "] CNV分析：", sample_id)
    sample_cells <- colnames(all_cells)[as.character(all_cells$sample_id) %in% sample_id]
    obs_cells <- candidate_cells[
      epithelial$sample_id[match(candidate_cells, colnames(epithelial))] %in% sample_id
    ]
    obs_cluster <- setNames(
      as.character(epithelial$epithelial_cluster_id[match(obs_cells, colnames(epithelial))]),
      obs_cells
    )
    reference <- choose_reference_cells(
      sample_id,
      all_cells,
      config$cnv$minimum_reference_cells,
      config$cnv$maximum_reference_cells
    )
    reference_cells <- setdiff(intersect(reference$cells, colnames(full_counts)), obs_cells)
    infercnv_evaluable <- identical(reference$status, "evaluable") &&
      length(reference_cells) >= config$cnv$minimum_reference_cells
    infercnv_status <- if (infercnv_evaluable) "evaluable" else "not_evaluable_insufficient_reference"
    sample_dir <- file.path(cnv_dir, sample_id)
    infer_dir <- file.path(sample_dir, "infercnv")
    copykat_dir <- file.path(sample_dir, "copykat")
    dir.create(infer_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(copykat_dir, recursive = TRUE, showWarnings = FALSE)

    input_cells <- c(obs_cells, reference_cells)
    checkpoint_path <- file.path(sample_dir, paste0(sample_id, "_F1_06_sample_result.rds"))
    checkpoint_signature <- list(
      format_version = "F1_06_sample_parallel_v2_single_observation_i6_copykat_ABC",
      sample_id = sample_id,
      all_sample_cells_md5 = f1_character_md5(sort(sample_cells)),
      observation_cells_md5 = f1_character_md5(sort(obs_cells)),
      reference_cells_md5 = f1_character_md5(sort(reference_cells)),
      infercnv_genes_md5 = f1_character_md5(matched_gene_order$gene_symbol),
      approved_clusters = paste(sort(approved_clusters), collapse = "|"),
      reference_source = reference$source,
      reference_status = reference$status,
      infercnv_observation_group = config$cnv$infercnv_observation_group,
      infercnv_analysis_mode = config$cnv$infercnv_analysis_mode,
      infercnv_partition_method = config$cnv$infercnv_tumor_subcluster_partition_method,
      infercnv_k_nn = config$cnv$infercnv_k_nn,
      infercnv_leiden_resolution = config$cnv$infercnv_leiden_resolution,
      infercnv_cutoff = config$cnv$infercnv_cutoff,
      infercnv_scipen = config$cnv$infercnv_scipen,
      infercnv_bitmap_type = config$cnv$infercnv_bitmap_type,
      infercnv_hmm = config$cnv$infercnv_hmm,
      infercnv_hmm_type = config$cnv$infercnv_hmm_type,
      infercnv_hmm_report_by = config$cnv$infercnv_hmm_report_by,
      infercnv_bayes_max_p_normal = config$cnv$infercnv_bayes_max_p_normal,
      infercnv_reassign_cnvs = config$cnv$infercnv_reassign_cnvs,
      infercnv_denoise = config$cnv$infercnv_denoise,
      copykat_ngene_chr = config$cnv$copykat_ngene_chr,
      copykat_internal_min_genes_per_cell =
        config$cnv$copykat_internal_min_genes_per_cell,
      copykat_data_source = copykat_data_source,
      copykat_win_size = config$cnv$copykat_win_size,
      copykat_ks_cut = config$cnv$copykat_ks_cut,
      copykat_arms = "A_self_estimated|B_same_sample_immune|C_normal_gastric_two_fold_holdout",
      copykat_holdout_seed = config$cnv$copykat_holdout_seed,
      copykat_baseline_suspect_epithelial_fraction =
        config$cnv$copykat_baseline_suspect_epithelial_fraction
    )
    if (file.exists(checkpoint_path)) {
      checkpoint <- tryCatch(readRDS(checkpoint_path), error = function(e) NULL)
      artifacts_exist <- !is.null(checkpoint) &&
        all(file.exists(as.character(checkpoint$required_artifacts)))
      if (
        !is.null(checkpoint) &&
          identical(checkpoint$signature, checkpoint_signature) &&
          artifacts_exist &&
          is.list(checkpoint$result)
      ) {
        message("[", stage, "] 复用样本检查点：", sample_id)
        checkpoint$result$checkpoint_reused <- TRUE
        return(checkpoint$result)
      }
    }

    if (infercnv_evaluable) {
      infer_counts <- full_counts[matched_gene_order$gene_symbol, input_cells, drop = FALSE]
      infer_score_row <- run_infercnv_matrix(
        counts = infer_counts,
        input_cells = input_cells,
        obs_cells = obs_cells,
        obs_cluster = obs_cluster,
        reference_cells = reference_cells,
        sample_id = sample_id,
        gene_order = matched_gene_order,
        out_dir = infer_dir,
        file_prefix = paste0(sample_id, "_infercnv_raw"),
        config = config
      )
      infer_score_row$infercnv_evaluation_status <- infercnv_status
      rm(infer_counts)
    } else {
      # reference不足只影响该样本的inferCNV；保留占位状态并继续其他样本及CopyKAT。
      infer_score_row <- data.frame(
        cell_id_final = obs_cells,
        sample_id = sample_id,
        epithelial_cluster = unname(obs_cluster[obs_cells]),
        infercnv_parent_group = config$cnv$infercnv_observation_group,
        infercnv_subcluster_local = "not_evaluable_sample",
        infercnv_subcluster_id = paste0(sample_id, "__not_evaluable_sample"),
        infercnv_cell_burden = NA_real_,
        reference_background_P95 = NA_real_,
        above_reference_P95 = NA,
        infercnv_evaluation_status = infercnv_status,
        stringsAsFactors = FALSE
      )
    }

    # 三个CopyKAT臂都查看当前样本全部QC后singlet，并在同一个样本任务内顺序运行。
    # inferCNV即使使用配对或其他正常胃样本reference，也绝不传给CopyKAT。
    copykat_input <- full_counts[, sample_cells, drop = FALSE]
    copykat_reference <- choose_copykat_known_normals(
      sample_id,
      all_cells,
      config$cnv$minimum_reference_cells,
      config$cnv$maximum_reference_cells
    )
    stale_failure_path <- file.path(sample_dir, paste0(sample_id, "_F1_06_failure.tsv"))

    # A臂：官方默认的样本内自估基线，必须完全省略norm.cell.names。
    copykat_a <- run_copykat_arm(
      copykat_input = copykat_input,
      sample_cells = sample_cells,
      evaluation_cells = obs_cells,
      sample_id = sample_id,
      arm_id = "A_self_estimated",
      arm_dir = file.path(copykat_dir, "A_self_estimated"),
      known_normal_cells = character(),
      seed = config$seed,
      config = config
    )
    copykat_a_audit <- audit_copykat_self_estimated_baseline(
      copykat_a$all_calls,
      obs_cells,
      sample_id,
      copykat_a$status,
      config$cnv$copykat_baseline_suspect_epithelial_fraction
    )
    copykat_a_eval <- copykat_a$evaluation_calls[
      match(obs_cells, copykat_a$evaluation_calls$cell_id_final),
      ,
      drop = FALSE
    ]

    # B臂：同样本免疫细胞足够时，显式提供known normal；不足时不把它伪装成A臂。
    if (length(copykat_reference$cells) >= config$cnv$minimum_reference_cells) {
      copykat_b <- run_copykat_arm(
        copykat_input = copykat_input,
        sample_cells = sample_cells,
        evaluation_cells = obs_cells,
        sample_id = sample_id,
        arm_id = "B_same_sample_immune",
        arm_dir = file.path(copykat_dir, "B_same_sample_immune"),
        known_normal_cells = copykat_reference$cells,
        seed = config$seed,
        config = config
      )
      copykat_b_eval <- copykat_b$evaluation_calls[
        match(obs_cells, copykat_b$evaluation_calls$cell_id_final),
        ,
        drop = FALSE
      ]
      copykat_b_status <- copykat_b$status
    } else {
      copykat_b <- list(required_artifacts = character(), known_normal_cells = character())
      copykat_b_eval <- empty_copykat_calls(
        obs_cells,
        sample_id,
        "B_same_sample_immune",
        "not_evaluable_insufficient_same_sample_immune_reference"
      )
      copykat_b_status <- "not_evaluable_insufficient_same_sample_immune_reference"
    }

    # C臂：只在Normal_Gastric中做两折留出。每个细胞只在未作为known normal的那一折评价。
    sample_group <- unique(as.character(
      all_cells$group_analysis[match(sample_cells, colnames(all_cells))]
    ))
    if (length(sample_group) != 1L || is.na(sample_group)) {
      stop(sample_id, "无法唯一确定group_analysis。")
    }
    copykat_c_reference_sets <- list()
    copykat_c_required_artifacts <- character()
    if (identical(sample_group, "Normal_Gastric") && length(obs_cells) >= 2L) {
      set.seed(config$cnv$copykat_holdout_seed)
      shuffled <- sample(sort(obs_cells), length(obs_cells), replace = FALSE)
      fold_assignment <- setNames(
        rep(1:2, length.out = length(shuffled)),
        shuffled
      )
      copykat_c_rows <- vector("list", 2L)
      copykat_c_status_rows <- character(2L)
      for (fold in 1:2) {
        held_out <- names(fold_assignment)[fold_assignment == fold]
        held_in <- setdiff(obs_cells, held_out)
        known_normal <- unique(c(copykat_reference$candidate_cells, held_in))
        arm_id <- paste0("C_normal_gastric_holdout_fold", fold)
        copykat_c_fold <- run_copykat_arm(
          copykat_input = copykat_input,
          sample_cells = sample_cells,
          evaluation_cells = held_out,
          sample_id = sample_id,
          arm_id = arm_id,
          arm_dir = file.path(copykat_dir, arm_id),
          known_normal_cells = known_normal,
          seed = config$cnv$copykat_holdout_seed,
          config = config
        )
        fold_calls <- copykat_c_fold$evaluation_calls
        fold_calls$copykat_C_holdout_fold <- fold
        fold_calls$copykat_C_holdout_status <- copykat_c_fold$status
        copykat_c_rows[[fold]] <- fold_calls
        copykat_c_status_rows[[fold]] <- copykat_c_fold$status
        copykat_c_reference_sets[[as.character(fold)]] <- known_normal
        copykat_c_required_artifacts <- c(
          copykat_c_required_artifacts,
          copykat_c_fold$required_artifacts
        )
        rm(copykat_c_fold)
        gc(verbose = FALSE)
      }
      copykat_c_eval <- do.call(rbind, copykat_c_rows)
      copykat_c_eval <- copykat_c_eval[
        match(obs_cells, copykat_c_eval$cell_id_final),
        ,
        drop = FALSE
      ]
      if (anyNA(copykat_c_eval$cell_id_final) || anyDuplicated(copykat_c_eval$cell_id_final)) {
        stop(sample_id, "的CopyKAT C臂未一对一覆盖正常胃候选上皮。")
      }
      copykat_c_status <- paste(unique(copykat_c_status_rows), collapse = "|")
    } else {
      copykat_c_status <- if (identical(sample_group, "Normal_Gastric")) {
        "not_evaluable_fewer_than_two_normal_epithelial_cells"
      } else {
        "not_applicable_non_normal_gastric"
      }
      copykat_c_eval <- empty_copykat_calls(
        obs_cells, sample_id, "C_normal_gastric_holdout", copykat_c_status
      )
      copykat_c_eval$copykat_C_holdout_fold <- NA_integer_
      copykat_c_eval$copykat_C_holdout_status <- copykat_c_status
    }

    copykat_row <- data.frame(
      cell_id_final = obs_cells,
      sample_id = sample_id,
      copykat_primary_raw_call = copykat_a_eval$copykat_call,
      copykat_call = if (isTRUE(copykat_a_audit$copykat_primary_evidence_usable)) {
        copykat_a_eval$copykat_call
      } else {
        rep("uncalled", length(obs_cells))
      },
      copykat_primary_returned = copykat_a_eval$copykat_returned,
      copykat_primary_note = copykat_a_eval$copykat_note,
      copykat_primary_run_status = copykat_a$status,
      copykat_primary_baseline_status =
        copykat_a_audit$copykat_primary_baseline_status,
      copykat_primary_evidence_usable =
        copykat_a_audit$copykat_primary_evidence_usable,
      copykat_B_immune_call = copykat_b_eval$copykat_call,
      copykat_B_immune_returned = copykat_b_eval$copykat_returned,
      copykat_B_immune_status = copykat_b_status,
      copykat_C_holdout_call = copykat_c_eval$copykat_call,
      copykat_C_holdout_returned = copykat_c_eval$copykat_returned,
      copykat_C_holdout_fold = copykat_c_eval$copykat_C_holdout_fold,
      copykat_C_holdout_status = copykat_c_eval$copykat_C_holdout_status,
      stringsAsFactors = FALSE
    )
    infercnv_reference_manifest <- f1_reference_manifest(
      sample_id,
      if (infercnv_evaluable) "infercnv_reference" else "infercnv_reference_candidate_not_used",
      reference_cells,
      reference$source,
      all_cells
    )
    copykat_b_reference_manifest <- f1_reference_manifest(
      sample_id,
      "copykat_B_same_sample_immune_known_normal",
      copykat_reference$cells,
      copykat_reference$source,
      all_cells
    )
    copykat_c_reference_manifest <- if (length(copykat_c_reference_sets)) {
      do.call(rbind, lapply(names(copykat_c_reference_sets), function(fold) {
        f1_reference_manifest(
          sample_id,
          paste0("copykat_C_fold", fold, "_known_normal"),
          copykat_c_reference_sets[[fold]],
          paste0("same_sample_immune_plus_held_in_normal_epithelium_fold", fold),
          all_cells
        )
      }))
    } else {
      f1_reference_manifest(
        sample_id,
        "copykat_C_not_run",
        character(),
        copykat_c_status,
        all_cells
      )
    }
    observation_ambient <- summarize_selected_ambient(obs_cells)
    infercnv_reference_ambient <- summarize_selected_ambient(reference_cells)
    copykat_reference_ambient <- summarize_selected_ambient(copykat_reference$cells)
    reference_row <- data.frame(
      analysis_contract_version = "F1_06_v2_single_observation_i6",
      sample_id = sample_id,
      target_patient_id = reference$target_patient_id,
      observation_cells = length(obs_cells),
      reference_cells = length(reference_cells),
      reference_source = reference$source,
      reference_source_samples = if (length(reference$source_samples)) {
        paste(reference$source_samples, collapse = "|")
      } else {
        "none"
      },
      reference_source_patients = if (length(reference$source_patients)) {
        paste(reference$source_patients, collapse = "|")
      } else {
        "none"
      },
      reference_lineages = paste(sort(unique(all_cells$cell_type_major[match(reference_cells, colnames(all_cells))])), collapse = "|"),
      observation_ambient_evaluable_cells = observation_ambient[["evaluable_cells"]],
      observation_contamination_median = observation_ambient[["contamination_median"]],
      observation_contamination_P90 = observation_ambient[["contamination_P90"]],
      infercnv_reference_ambient_evaluable_cells =
        infercnv_reference_ambient[["evaluable_cells"]],
      infercnv_reference_contamination_median =
        infercnv_reference_ambient[["contamination_median"]],
      infercnv_reference_contamination_P90 =
        infercnv_reference_ambient[["contamination_P90"]],
      same_sample_reference_candidates = reference$tier_counts[["same_sample"]],
      same_patient_normal_gastric_reference_candidates =
        reference$tier_counts[["same_patient_normal_gastric"]],
      other_patient_normal_gastric_reference_candidates =
        reference$tier_counts[["balanced_other_patient_normal_gastric"]],
      reference_minimum_cells = config$cnv$minimum_reference_cells,
      reference_maximum_cells = config$cnv$maximum_reference_cells,
      infercnv_reference_policy = infercnv_reference_policy,
      infercnv_reference_pairing_rule = "exact_patient_id_and_Normal_Gastric_only",
      infercnv_observation_group = config$cnv$infercnv_observation_group,
      copykat_input_scope = copykat_input_scope,
      copykat_input_cells = length(sample_cells),
      copykat_known_normal_cells = length(copykat_reference$cells),
      copykat_reference_ambient_evaluable_cells =
        copykat_reference_ambient[["evaluable_cells"]],
      copykat_reference_contamination_median =
        copykat_reference_ambient[["contamination_median"]],
      copykat_reference_contamination_P90 =
        copykat_reference_ambient[["contamination_P90"]],
      copykat_same_sample_normal_candidates =
        copykat_reference$available_same_sample_candidates,
      copykat_baseline_mode = "A_self_estimated_without_norm.cell.names",
      copykat_primary_run_status = copykat_a$status,
      copykat_primary_baseline_status =
        copykat_a_audit$copykat_primary_baseline_status,
      copykat_primary_evidence_usable =
        copykat_a_audit$copykat_primary_evidence_usable,
      copykat_B_immune_status = copykat_b_status,
      copykat_B_known_normal_cells = length(copykat_reference$cells),
      copykat_C_holdout_status = copykat_c_status,
      copykat_C_holdout_seed = config$cnv$copykat_holdout_seed,
      copykat_external_reference = "prohibited",
      copykat_internal_min_genes_per_cell =
        config$cnv$copykat_internal_min_genes_per_cell,
      copykat_data_source = copykat_data_source,
      infercnv_evaluation_status = infercnv_status,
      infercnv_cutoff = config$cnv$infercnv_cutoff,
      infercnv_scipen = config$cnv$infercnv_scipen,
      infercnv_bitmap_type = config$cnv$infercnv_bitmap_type,
      infercnv_HMM = config$cnv$infercnv_hmm,
      infercnv_HMM_type = config$cnv$infercnv_hmm_type,
      infercnv_HMM_report_by = config$cnv$infercnv_hmm_report_by,
      infercnv_BayesMaxPNormal = config$cnv$infercnv_bayes_max_p_normal,
      infercnv_reassignCNVs = config$cnv$infercnv_reassign_cnvs,
      infercnv_denoise = config$cnv$infercnv_denoise,
      infercnv_analysis_mode = config$cnv$infercnv_analysis_mode,
      infercnv_partition_method =
        config$cnv$infercnv_tumor_subcluster_partition_method,
      infercnv_k_nn = config$cnv$infercnv_k_nn,
      infercnv_leiden_resolution =
        as.character(config$cnv$infercnv_leiden_resolution),
      infercnv_leiden_method = config$cnv$infercnv_leiden_method,
      infercnv_leiden_function = config$cnv$infercnv_leiden_function,
      infercnv_inspect_subclusters = config$cnv$infercnv_inspect_subclusters,
      infercnv_internal_subclustering = identical(
        config$cnv$infercnv_analysis_mode,
        "subclusters"
      ),
      infercnv_min_cells_per_gene = config$cnv$infercnv_min_cells_per_gene,
      infercnv_window_length = config$cnv$infercnv_window_length,
      infercnv_write_expr_matrix = config$cnv$infercnv_write_expr_matrix,
      infercnv_output_format = config$cnv$infercnv_output_format,
      stringsAsFactors = FALSE
    )
    required_artifacts <- c(
      copykat_a$required_artifacts,
      copykat_b$required_artifacts,
      copykat_c_required_artifacts
    )
    if (infercnv_evaluable) {
      required_artifacts <- c(
        required_artifacts,
        file.path(infer_dir, paste0(sample_id, "_infercnv_raw_final_object.rds")),
        file.path(infer_dir, paste0(sample_id, "_infercnv_raw_subcluster_membership.tsv")),
        file.path(infer_dir, paste0(sample_id, "_infercnv_raw_plot_source_manifest.tsv"))
      )
    }
    copykat_statuses <- c(copykat_a$status, copykat_b_status, copykat_c_status)
    sample_execution_status <- if (any(grepl("CNA_and_default_heatmap_unavailable", copykat_statuses))) {
      "completed_with_copykat_CNA_unavailable"
    } else {
      "completed"
    }
    sample_result <- list(
      infer_scores = infer_score_row,
      copykat_calls = copykat_row,
      copykat_baseline_audit = copykat_a_audit,
      reference_summary = reference_row,
      reference_manifest = rbind(
        infercnv_reference_manifest,
        copykat_b_reference_manifest,
        copykat_c_reference_manifest
      ),
      sample_execution_status = sample_execution_status,
      log_note = if (identical(sample_execution_status, "completed_with_copykat_CNA_unavailable")) {
        paste0(
          sample_id,
          "的inferCNV/CopyKAT已完成；至少一个CopyKAT臂只恢复了官方prediction，",
          "对应最终CNA_results/默认热图不可评估，且不作为主恶性支持。"
        )
      } else if (infercnv_evaluable) {
        paste0(sample_id, "的inferCNV和CopyKAT A/B/C适用分支计算完成。")
      } else {
        paste0(sample_id, "的合格免疫reference不足，inferCNV记为not_evaluable；CopyKAT适用分支已完成。")
      },
      checkpoint_reused = FALSE
    )
    saveRDS(
      list(
        signature = checkpoint_signature,
        required_artifacts = required_artifacts,
        result = sample_result
      ),
      checkpoint_path,
      compress = TRUE
    )
    if (file.exists(stale_failure_path)) unlink(stale_failure_path)
    rm(
      copykat_input, copykat_reference, copykat_a, copykat_b,
      copykat_a_eval, copykat_b_eval, copykat_c_eval
    )
    gc(verbose = FALSE)
    sample_result
  }

  run_cnv_sample_safely <- function(sample_id) {
    tryCatch(
      list(ok = TRUE, result = run_cnv_sample(sample_id), error = ""),
      error = function(e) {
        error_message <- conditionMessage(e)
        failure_dir <- file.path(cnv_dir, sample_id)
        dir.create(failure_dir, recursive = TRUE, showWarnings = FALSE)
        failure_path <- file.path(
          failure_dir,
          paste0(sample_id, "_F1_06_failure.tsv")
        )
        f1_write_tsv(
          data.frame(
            sample_id = sample_id,
            failed_at = f1_now(),
            error_message = error_message,
            stringsAsFactors = FALSE
          ),
          failure_path
        )
        message("[", stage, "] 样本失败：", sample_id, " | ", error_message)
        list(ok = FALSE, result = NULL, error = error_message)
      }
    )
  }
  sample_workers <- if (.Platform$OS.type == "windows") {
    1L
  } else {
    min(config$cnv$sample_workers, length(sample_ids))
  }
  message(
    "[", stage, "] 样本级并行：", sample_workers,
    "个样本；每个样本inferCNV=", config$cnv$infercnv_threads,
    "线程，CopyKAT=", config$cnv$copykat_cores, "核。"
  )
  sample_runs <- if (sample_workers > 1L) {
    parallel::mclapply(
      sample_ids,
      run_cnv_sample_safely,
      mc.cores = sample_workers,
      mc.preschedule = FALSE,
      mc.set.seed = FALSE,
      mc.cleanup = TRUE
    )
  } else {
    lapply(sample_ids, run_cnv_sample_safely)
  }
  names(sample_runs) <- sample_ids
  sample_ok <- vapply(
    sample_runs,
    function(x) is.list(x) && isTRUE(x$ok),
    logical(1)
  )
  sample_status <- data.frame(
    sample_id = sample_ids,
    status = vapply(
      sample_runs,
      function(x) {
        if (!is.list(x) || !isTRUE(x$ok)) return("failed")
        as.character(x$result$sample_execution_status %||% "completed")
      },
      character(1)
    ),
    checkpoint_reused = vapply(
      sample_runs,
      function(x) {
        if (!is.list(x) || !isTRUE(x$ok)) return(NA)
        isTRUE(x$result$checkpoint_reused)
      },
      logical(1)
    ),
    error_message = vapply(
      sample_runs,
      function(x) {
        if (!is.list(x)) return("worker_returned_invalid_result")
        as.character(x$error %||% "")
      },
      character(1)
    ),
    stringsAsFactors = FALSE
  )
  sample_status_path <- file.path(cnv_dir, "F1_CNV_sample_execution_status.tsv")
  f1_write_tsv(sample_status, sample_status_path)
  if (any(!sample_ok)) {
    failed <- paste(sample_ids[!sample_ok], collapse = ", ")
    stop(
      "以下样本的CNV计算失败：", failed,
      "。已完成样本的检查点保留；详见：", sample_status_path
    )
  }

  sample_results <- lapply(sample_runs, `[[`, "result")
  for (result in sample_results) {
    f1_append_log(config, stage, result$log_note)
  }
  infer_scores_all <- do.call(rbind, lapply(sample_results, `[[`, "infer_scores"))
  copykat_calls_all <- do.call(rbind, lapply(sample_results, `[[`, "copykat_calls"))
  copykat_baseline_audit_all <- do.call(
    rbind,
    lapply(sample_results, `[[`, "copykat_baseline_audit")
  )
  reference_summary_all <- do.call(rbind, lapply(sample_results, `[[`, "reference_summary"))
  reference_manifest_all <- do.call(rbind, lapply(sample_results, `[[`, "reference_manifest"))
  subcluster_membership_all <- infer_scores_all[, c(
    "cell_id_final", "sample_id", "epithelial_cluster", "infercnv_parent_group",
    "infercnv_subcluster_local", "infercnv_subcluster_id",
    "infercnv_evaluation_status"
  )]
  f1_write_tsv(infer_scores_all, file.path(cnv_dir, "infercnv_cell_scores_all.tsv"))
  f1_write_tsv(copykat_calls_all, file.path(cnv_dir, "copykat_cell_calls_all.tsv"))
  f1_write_tsv(copykat_baseline_audit_all, copykat_baseline_audit_path)
  f1_write_tsv(reference_summary_all, file.path(cnv_dir, "infercnv_reference_summary.tsv"))
  f1_write_tsv(subcluster_membership_all, subcluster_membership_path)
  f1_write_tsv(reference_manifest_all, reference_manifest_path)
}

f1_require_columns(
  infer_scores_all,
  c(
    "cell_id_final", "sample_id", "epithelial_cluster", "infercnv_parent_group",
    "infercnv_subcluster_local", "infercnv_subcluster_id",
    "infercnv_cell_burden", "reference_background_P95", "above_reference_P95",
    "infercnv_evaluation_status"
  ),
  "infercnv_cell_scores_all.tsv"
)
f1_require_columns(
  copykat_calls_all,
  c(
    "cell_id_final", "sample_id", "copykat_primary_raw_call", "copykat_call",
    "copykat_primary_baseline_status", "copykat_primary_evidence_usable",
    "copykat_B_immune_call", "copykat_B_immune_status",
    "copykat_C_holdout_call", "copykat_C_holdout_status"
  ),
  "copykat_cell_calls_all.tsv"
)
f1_require_columns(
  copykat_baseline_audit_all,
  c(
    "sample_id", "predicted_diploid_cells",
    "predicted_diploid_candidate_epithelial_fraction",
    "copykat_primary_baseline_status", "copykat_primary_evidence_usable"
  ),
  "copykat_self_estimated_diploid_composition_all.tsv"
)
f1_require_columns(
  reference_summary_all,
  c(
    "sample_id", "reference_cells", "reference_source",
    "infercnv_reference_policy", "infercnv_evaluation_status",
    "copykat_input_scope", "copykat_external_reference"
  ),
  "infercnv_reference_summary.tsv"
)
f1_require_columns(
  subcluster_membership_all,
  c(
    "cell_id_final", "sample_id", "epithelial_cluster",
    "infercnv_subcluster_id", "infercnv_evaluation_status"
  ),
  "infercnv_subcluster_membership_all.tsv"
)
f1_require_columns(
  reference_manifest_all,
  c(
    "target_sample_id", "method", "reference_cell_id", "reference_sample_id",
    "reference_patient_id", "reference_group_analysis",
    "reference_cell_type_major", "annotation_confidence", "reference_source"
  ),
  "cnv_reference_cell_manifest.tsv"
)
copykat_manifest <- reference_manifest_all[
  grepl("^copykat_[BC]_", reference_manifest_all$method),
  ,
  drop = FALSE
]
if (
  nrow(copykat_manifest) &&
    any(copykat_manifest$target_sample_id != copykat_manifest$reference_sample_id)
) {
  stop("CopyKAT已知正常细胞中出现外部样本，违反current-sample-only规则。")
}
copykat_c_manifest <- copykat_manifest[
  grepl("^copykat_C_", copykat_manifest$method),
  ,
  drop = FALSE
]
if (
  nrow(copykat_c_manifest) &&
    any(copykat_c_manifest$reference_group_analysis != "Normal_Gastric")
) {
  stop("CopyKAT C臂known normal只能来自当前Normal_Gastric样本。")
}
if (!setequal(candidate_cells, infer_scores_all$cell_id_final)) {
  stop("inferCNV cell score没有完整覆盖获准的上皮候选细胞。")
}
if (!setequal(candidate_cells, copykat_calls_all$cell_id_final)) {
  stop("CopyKAT call表没有完整覆盖获准的上皮候选细胞。")
}
if (
  anyDuplicated(subcluster_membership_all$cell_id_final) ||
    !setequal(candidate_cells, subcluster_membership_all$cell_id_final)
) {
  stop("inferCNV subcluster membership必须一对一完整覆盖获准的上皮候选细胞。")
}
f1_write_tsv(
  subcluster_membership_all,
  file.path(config$paths$malignancy_dir, "infercnv_subcluster_membership.tsv")
)

infer_subcluster_summary <- do.call(rbind, lapply(
  split(
    infer_scores_all,
    interaction(
      infer_scores_all$sample_id,
      infer_scores_all$infercnv_subcluster_id,
      drop = TRUE
    )
  ),
  function(x) {
    evaluable <- identical(unique(x$infercnv_evaluation_status), "evaluable")
    cluster_counts <- sort(table(as.character(x$epithelial_cluster)), decreasing = TRUE)
    dominant_cluster <- names(cluster_counts)[[1]]
    data.frame(
      sample_id = x$sample_id[[1]],
      infercnv_subcluster_id = x$infercnv_subcluster_id[[1]],
      infercnv_subcluster_local = x$infercnv_subcluster_local[[1]],
      epithelial_cluster = dominant_cluster,
      dominant_epithelial_cluster_fraction =
        as.numeric(cluster_counts[[1]]) / sum(cluster_counts),
      epithelial_cluster_count = length(cluster_counts),
      epithelial_cluster_composition = paste0(
        names(cluster_counts), "=", as.integer(cluster_counts), collapse = "|"
      ),
      cell_count = nrow(x),
      infercnv_burden_median = if (evaluable) stats::median(x$infercnv_cell_burden, na.rm = TRUE) else NA_real_,
      reference_background_P95 = if (evaluable) unique(x$reference_background_P95)[[1]] else NA_real_,
      fraction_above_reference_P95 = if (evaluable) {
        mean(f1_as_logical(x$above_reference_P95, "above_reference_P95"), na.rm = TRUE)
      } else {
        NA_real_
      },
      infercnv_evaluation_status = unique(x$infercnv_evaluation_status)[[1]],
      subcluster_size_note = if (nrow(x) == 1L) {
        "singleton_not_evaluable_on_its_own"
      } else {
        "evaluate_from_heatmap_without_target_subcluster_count"
      },
      broad_segment_support = if (!evaluable || nrow(x) == 1L) {
        "not_evaluable"
      } else {
        "pending_researcher_heatmap_review"
      },
      stringsAsFactors = FALSE
    )
  }
))
reference_context_fields <- c(
  "sample_id", "target_patient_id", "reference_cells", "reference_source",
  "reference_source_samples", "reference_source_patients", "reference_lineages"
)
f1_require_columns(
  reference_summary_all,
  reference_context_fields,
  "infercnv_reference_summary.tsv"
)
infer_subcluster_summary <- merge(
  infer_subcluster_summary,
  reference_summary_all[, reference_context_fields],
  by = "sample_id",
  all.x = TRUE,
  sort = FALSE
)
f1_write_tsv(
  infer_subcluster_summary,
  file.path(config$paths$malignancy_dir, "infercnv_subcluster_summary.tsv")
)
f1_write_tsv(copykat_calls_all, file.path(config$paths$malignancy_dir, "copykat_cell_calls.tsv"))
f1_write_tsv(
  copykat_baseline_audit_all,
  file.path(
    config$paths$malignancy_dir,
    "copykat_self_estimated_diploid_composition.tsv"
  )
)

# 三臂都按候选上皮细胞统计，分母明确为该组织组中进入F1.6的候选上皮细胞。
copykat_group <- as.character(
  epithelial$group_analysis[match(copykat_calls_all$cell_id_final, colnames(epithelial))]
)
copykat_arm_columns <- c(
  A_self_estimated_raw = "copykat_primary_raw_call",
  B_same_sample_immune = "copykat_B_immune_call",
  C_normal_gastric_two_fold_holdout = "copykat_C_holdout_call"
)
copykat_three_arm_summary <- do.call(rbind, lapply(
  names(copykat_arm_columns),
  function(arm_name) {
    call_column <- copykat_arm_columns[[arm_name]]
    do.call(rbind, lapply(
      c("Normal_Gastric", "Primary_Tumor", "Peritoneal_Metastasis"),
      function(group_name) {
        keep <- copykat_group == group_name
        calls <- as.character(copykat_calls_all[[call_column]][keep])
        n_cells <- length(calls)
        n_aneuploid <- sum(calls == "aneuploid", na.rm = TRUE)
        n_diploid <- sum(calls == "diploid", na.rm = TRUE)
        n_uncalled <- n_cells - n_aneuploid - n_diploid
        n_called <- n_aneuploid + n_diploid
        data.frame(
          copykat_arm = arm_name,
          group_analysis = group_name,
          denominator_candidate_epithelial_cells = n_cells,
          aneuploid_cells = n_aneuploid,
          diploid_cells = n_diploid,
          uncalled_cells = n_uncalled,
          aneuploid_fraction_all_candidate_cells =
            if (n_cells) n_aneuploid / n_cells else NA_real_,
          aneuploid_fraction_called_cells =
            if (n_called) n_aneuploid / n_called else NA_real_,
          interpretation = if (
            identical(arm_name, "C_normal_gastric_two_fold_holdout") &&
              !identical(group_name, "Normal_Gastric")
          ) {
            "not_applicable_C_arm_is_normal_gastric_only"
          } else {
            "descriptive_method_comparison_not_independent_validation"
          },
          stringsAsFactors = FALSE
        )
      }
    ))
  }
))
f1_write_tsv(
  copykat_three_arm_summary,
  file.path(config$paths$malignancy_dir, "copykat_three_arm_group_comparison.tsv")
)

# 检查正常胃分泌/消化程序是否与CopyKAT A臂的aneuploid分类系统性相关。
# 使用raw counts中的该panel占总UMI比例；统计单位为样本内两类细胞的中位数差。
secretory_genes <- intersect(config$cnv$copykat_secretory_genes, rownames(full_counts))
normal_copykat_rows <- copykat_calls_all[
  copykat_group == "Normal_Gastric" &
    copykat_calls_all$copykat_primary_raw_call %in% c("aneuploid", "diploid"),
  ,
  drop = FALSE
]
if (length(secretory_genes) && nrow(normal_copykat_rows)) {
  secretory_cells <- normal_copykat_rows$cell_id_final
  secretory_fraction <- Matrix::colSums(
    full_counts[secretory_genes, secretory_cells, drop = FALSE]
  ) / pmax(Matrix::colSums(full_counts[, secretory_cells, drop = FALSE]), 1)
  secretory_cell_source <- data.frame(
    cell_id_final = secretory_cells,
    sample_id = normal_copykat_rows$sample_id,
    copykat_primary_raw_call = normal_copykat_rows$copykat_primary_raw_call,
    secretory_gene_UMI_fraction = as.numeric(secretory_fraction[secretory_cells]),
    stringsAsFactors = FALSE
  )
  secretory_by_sample <- do.call(rbind, lapply(
    split(secretory_cell_source, secretory_cell_source$sample_id),
    function(x) {
      aneuploid_values <- x$secretory_gene_UMI_fraction[
        x$copykat_primary_raw_call == "aneuploid"
      ]
      diploid_values <- x$secretory_gene_UMI_fraction[
        x$copykat_primary_raw_call == "diploid"
      ]
      data.frame(
        sample_id = x$sample_id[[1]],
        aneuploid_cells = length(aneuploid_values),
        diploid_cells = length(diploid_values),
        median_secretory_fraction_aneuploid = if (length(aneuploid_values)) {
          stats::median(aneuploid_values)
        } else {
          NA_real_
        },
        median_secretory_fraction_diploid = if (length(diploid_values)) {
          stats::median(diploid_values)
        } else {
          NA_real_
        },
        stringsAsFactors = FALSE
      )
    }
  ))
  secretory_by_sample$aneuploid_minus_diploid_median <-
    secretory_by_sample$median_secretory_fraction_aneuploid -
      secretory_by_sample$median_secretory_fraction_diploid
  evaluable_difference <- is.finite(secretory_by_sample$aneuploid_minus_diploid_median)
  paired_sample_count <- sum(evaluable_difference)
  paired_p <- if (paired_sample_count >= 3L) {
    suppressWarnings(stats::wilcox.test(
      secretory_by_sample$aneuploid_minus_diploid_median[evaluable_difference],
      mu = 0,
      exact = FALSE
    )$p.value)
  } else {
    NA_real_
  }
  secretory_by_sample$panel_genes_present <- paste(secretory_genes, collapse = "|")
  secretory_by_sample$paired_evaluable_samples <- paired_sample_count
  secretory_by_sample$paired_wilcoxon_p <- paired_p
  secretory_by_sample$interpretation <-
    "method_behavior_check_only_not_malignancy_evidence"
} else {
  secretory_by_sample <- data.frame(
    sample_id = "not_evaluable",
    aneuploid_cells = 0L,
    diploid_cells = 0L,
    median_secretory_fraction_aneuploid = NA_real_,
    median_secretory_fraction_diploid = NA_real_,
    aneuploid_minus_diploid_median = NA_real_,
    panel_genes_present = if (length(secretory_genes)) {
      paste(secretory_genes, collapse = "|")
    } else {
      "none"
    },
    paired_evaluable_samples = 0L,
    paired_wilcoxon_p = NA_real_,
    interpretation = "not_evaluable_missing_genes_or_called_normal_cells",
    stringsAsFactors = FALSE
  )
}
f1_write_tsv(
  secretory_by_sample,
  file.path(config$paths$malignancy_dir, "copykat_secretory_spike_test.tsv")
)

copykat_with_subcluster <- merge(
  copykat_calls_all,
  subcluster_membership_all[, c(
    "cell_id_final", "epithelial_cluster", "infercnv_subcluster_id"
  )],
  by = "cell_id_final",
  all.x = TRUE,
  sort = FALSE
)
if (anyNA(copykat_with_subcluster$infercnv_subcluster_id)) {
  stop("CopyKAT候选细胞无法完整映射inferCNV subcluster。")
}
copykat_subcluster <- do.call(rbind, lapply(
  split(
    copykat_with_subcluster,
    interaction(
      copykat_with_subcluster$sample_id,
      copykat_with_subcluster$infercnv_subcluster_id,
      drop = TRUE
    )
  ),
  function(x) {
    data.frame(
      sample_id = x$sample_id[[1]],
      infercnv_subcluster_id = x$infercnv_subcluster_id[[1]],
      copykat_aneuploid = sum(x$copykat_call == "aneuploid"),
      copykat_diploid = sum(x$copykat_call == "diploid"),
      copykat_uncalled = sum(x$copykat_call == "uncalled"),
      stringsAsFactors = FALSE
    )
  }
))

# 这两个“程序支持”字段只作单方CNV证据时的生物学上下文。
# 谱系marker说明细胞像哪类胃上皮，不自动说明良恶性；审核表同时给出
# 候选inferCNV subcluster和同谱系Normal_Gastric的marker检出比例，供人工比较。
epithelial_context <- f1_read_tsv(config$paths$epithelial_review_template)
f1_require_columns(
  epithelial_context,
  c("epithelial_cluster", "group_cell_counts", "top_markers"),
  "F1_epithelial_cluster_review_template.tsv"
)
epithelial_context$epithelial_cluster <- as.character(epithelial_context$epithelial_cluster)
epithelial_review$epithelial_cluster <- as.character(epithelial_review$epithelial_cluster)
epithelial_context <- merge(
  epithelial_context[, c("epithelial_cluster", "group_cell_counts", "top_markers")],
  epithelial_review[, c("epithelial_cluster", "epithelial_subtype")],
  by = "epithelial_cluster",
  all.x = TRUE,
  sort = FALSE
)

marker_panel <- f1_read_tsv(config$paths$marker_panel)
f1_require_columns(
  marker_panel,
  c("cell_type", "positive_markers", "supporting_markers"),
  "cell_type_marker_panel.tsv"
)
lineage_panel <- marker_panel[grepl("_epithelial$", marker_panel$cell_type), , drop = FALSE]
lineage_marker_map <- setNames(
  lapply(seq_len(nrow(lineage_panel)), function(i) f1_lineage_markers(lineage_panel[i, , drop = FALSE])),
  lineage_panel$cell_type
)

review_template <- merge(
  infer_subcluster_summary,
  copykat_subcluster,
  by = c("sample_id", "infercnv_subcluster_id"),
  all = TRUE,
  sort = TRUE
)
review_template <- merge(
  review_template,
  epithelial_context,
  by = "epithelial_cluster",
  all.x = TRUE,
  sort = FALSE
)

program_context <- do.call(rbind, lapply(seq_len(nrow(review_template)), function(i) {
  sample_id <- as.character(review_template$sample_id[[i]])
  subcluster_id <- as.character(review_template$infercnv_subcluster_id[[i]])
  subtype <- as.character(review_template$epithelial_subtype[[i]])
  if (is.na(subtype)) subtype <- ""
  genes <- lineage_marker_map[[subtype]]
  sample_cluster_cells <- subcluster_membership_all$cell_id_final[
    as.character(subcluster_membership_all$sample_id) %in% sample_id &
      as.character(subcluster_membership_all$infercnv_subcluster_id) %in% subcluster_id
  ]
  same_lineage_normal_cells <- candidate_cells[
    as.character(
      epithelial$group_analysis[match(candidate_cells, colnames(epithelial))]
    ) %in% "Normal_Gastric" &
      as.character(
        epithelial$epithelial_subtype[match(candidate_cells, colnames(epithelial))]
      ) %in% subtype
  ]
  same_lineage_normal_samples <- unique(as.character(
    epithelial$sample_id[match(same_lineage_normal_cells, colnames(epithelial))]
  ))
  missing_genes <- setdiff(genes, rownames(full_counts))
  reference_status <- if (!length(genes)) {
    "not_evaluable_unrecognized_epithelial_subtype"
  } else if (!length(same_lineage_normal_cells)) {
    "not_evaluable_no_same_lineage_normal_reference"
  } else {
    "available"
  }
  interpretation <- if (identical(subtype, "Intestinal_like_epithelial")) {
    "intestinal_like_is_metaplasia_or_tumor_ambiguous"
  } else {
    "lineage_identity_does_not_prove_non_malignancy"
  }
  data.frame(
    canonical_lineage_markers = if (length(genes)) paste(genes, collapse = ",") else "not_mapped",
    lineage_markers_missing_from_matrix =
      if (length(missing_genes)) paste(missing_genes, collapse = ",") else "none",
    candidate_lineage_marker_detection_pct =
      f1_marker_detection_summary(full_counts, sample_cluster_cells, genes),
    same_lineage_normal_reference_cells = length(same_lineage_normal_cells),
    same_lineage_normal_reference_samples = length(same_lineage_normal_samples),
    normal_reference_lineage_marker_detection_pct =
      f1_marker_detection_summary(full_counts, same_lineage_normal_cells, genes),
    same_lineage_normal_reference_status = reference_status,
    lineage_interpretation = interpretation,
    stringsAsFactors = FALSE
  )
}))
review_template <- cbind(review_template, program_context)

ambient_summary <- f1_read_tsv(config$paths$ambient_summary)
f1_require_columns(
  ambient_summary,
  c(
    "sample_id", "status", "contamination_median", "contamination_Q3",
    "contamination_max", "corrected_counts_path"
  ),
  "ambient_rna_summary_by_sample.tsv"
)
ambient_context <- data.frame(
  sample_id = as.character(ambient_summary$sample_id),
  DecontX_status = ambient_summary$status,
  DecontX_contamination_median = ambient_summary$contamination_median,
  DecontX_contamination_Q3 = ambient_summary$contamination_Q3,
  DecontX_contamination_max = ambient_summary$contamination_max,
  DecontX_corrected_counts_path = ambient_summary$corrected_counts_path,
  stringsAsFactors = FALSE
)
review_template <- merge(
  review_template,
  ambient_context,
  by = "sample_id",
  all.x = TRUE,
  sort = FALSE
)

review_template$infercnv_broad_segment_support <- ifelse(
  review_template$infercnv_evaluation_status == "evaluable" &
    review_template$cell_count > 1L,
  "",
  "not_evaluable"
)
review_template$heatmap_broad_segment_description <- ""
review_template$chromosome_regions <- ""
review_template$normal_gastric_epithelial_comparison <- ""
review_template$hmm_segment_support_description <- ""
review_template$run_corrected_infercnv_sensitivity <- FALSE
review_template$corrected_infercnv_broad_segment_support <- "not_run"
review_template$corrected_sensitivity_reason <- ""
review_template$tumor_program_support <- ""
review_template$normal_program_support <- ""
review_template$review_note <- ""
f1_write_tsv(review_template, config$paths$malignancy_review_template)

if (!file.exists(config$paths$malignancy_review_approved)) {
  f1_append_log(config, stage, paste0("CNV结果已生成，等待热图和marker联合审核：", config$paths$malignancy_review_template))
  stop(
    "F1.6已生成inferCNV/CopyKAT摘要和恶性判定模板。请审核热图及marker后建立：",
    config$paths$malignancy_review_approved
  )
}

approved <- f1_read_tsv(config$paths$malignancy_review_approved)
review_fields <- c(
  "sample_id", "infercnv_subcluster_id", "epithelial_cluster",
  "infercnv_broad_segment_support",
  "heatmap_broad_segment_description", "chromosome_regions",
  "normal_gastric_epithelial_comparison", "hmm_segment_support_description",
  "run_corrected_infercnv_sensitivity", "corrected_infercnv_broad_segment_support",
  "corrected_sensitivity_reason",
  "tumor_program_support", "normal_program_support", "review_note"
)
f1_require_columns(approved, review_fields, "F1_malignancy_cluster_review_approved.tsv")
approved$key <- paste(approved$sample_id, approved$infercnv_subcluster_id, sep = "__")
template_keys <- paste(
  review_template$sample_id,
  review_template$infercnv_subcluster_id,
  sep = "__"
)
if (anyDuplicated(approved$key) || !setequal(approved$key, template_keys)) {
  stop("恶性审核批准表必须与当前sample_id × inferCNV subcluster一对一对应。")
}
template_index <- match(approved$key, template_keys)
if (any(
  as.character(approved$epithelial_cluster) !=
    as.character(review_template$epithelial_cluster[template_index])
)) {
  stop("恶性审核批准表不得修改inferCNV subcluster对应的父级epithelial cluster。")
}
allowed_infer <- c("strong", "weak", "absent", "not_evaluable")
if (any(!approved$infercnv_broad_segment_support %in% allowed_infer)) {
  stop("infercnv_broad_segment_support只允许strong、weak、absent或not_evaluable。")
}
evaluable_review <- approved$infercnv_broad_segment_support %in% c("strong", "weak", "absent")
required_heatmap_notes <- c(
  "heatmap_broad_segment_description",
  "chromosome_regions",
  "normal_gastric_epithelial_comparison",
  "hmm_segment_support_description"
)
for (field in required_heatmap_notes) {
  values <- trimws(as.character(approved[[field]]))
  if (any(evaluable_review & (is.na(values) | !nzchar(values)))) {
    stop("可评价subcluster必须填写", field, "，不能只凭自动标签判恶性。")
  }
}
approved$run_corrected_infercnv_sensitivity <- f1_as_logical(
  approved$run_corrected_infercnv_sensitivity,
  "run_corrected_infercnv_sensitivity"
)
reference_status <- reference_summary_all$infercnv_evaluation_status[
  match(approved$sample_id, reference_summary_all$sample_id)
]
if (anyNA(reference_status)) stop("恶性审核表中存在无法匹配reference状态的sample_id。")
if (any(
  reference_status != "evaluable" &
    approved$infercnv_broad_segment_support != "not_evaluable"
)) {
  stop("reference不足的样本必须把infercnv_broad_segment_support保留为not_evaluable。")
}
approved$tumor_program_support <- f1_as_logical(approved$tumor_program_support, "tumor_program_support")
approved$normal_program_support <- f1_as_logical(approved$normal_program_support, "normal_program_support")

# corrected分支只在审核者已根据ambient摘要和marker泄漏证据明确触发时运行。
# 开关按sample一致；CopyKAT仍使用raw counts。
corrected_request <- lapply(
  split(approved$run_corrected_infercnv_sensitivity, approved$sample_id),
  unique
)
if (any(lengths(corrected_request) != 1L)) {
  stop("同一sample_id的run_corrected_infercnv_sensitivity必须一致。")
}
corrected_samples <- names(corrected_request)[vapply(
  corrected_request,
  function(x) isTRUE(x[[1]]),
  logical(1)
)]
if (length(corrected_samples)) {
  requested_rows <- approved$sample_id %in% corrected_samples
  sensitivity_reason <- as.character(approved$corrected_sensitivity_reason[requested_rows])
  if (any(is.na(sensitivity_reason) | !nzchar(trimws(sensitivity_reason)))) {
    stop("运行corrected inferCNV的样本必须填写污染证据和触发理由。")
  }
  requested_status <- reference_summary_all$infercnv_evaluation_status[
    match(corrected_samples, reference_summary_all$sample_id)
  ]
  if (any(requested_status != "evaluable")) {
    stop("reference不足的样本不能运行corrected inferCNV敏感性。")
  }

  ambient_summary <- f1_read_tsv(config$paths$ambient_summary)
  f1_require_columns(
    ambient_summary,
    c("sample_id", "corrected_counts_path"),
    "ambient_rna_summary_by_sample.tsv"
  )
  corrected_scores_path <- file.path(cnv_dir, "infercnv_corrected_cell_scores_all.tsv")
  corrected_scores <- if (file.exists(corrected_scores_path)) {
    f1_read_tsv(corrected_scores_path)
  } else {
    NULL
  }
  if (!is.null(corrected_scores)) {
    corrected_provenance_fields <- c(
      "raw_infercnv_subcluster_id", "corrected_infercnv_subcluster_id",
      "infercnv_observation_group", "infercnv_HMM", "infercnv_HMM_type",
      "infercnv_analysis_mode", "infercnv_partition_method",
      "infercnv_k_nn", "infercnv_leiden_resolution",
      "infercnv_reference_policy", "infercnv_reference_source"
    )
    corrected_provenance_ok <- isTRUE(
      all(corrected_provenance_fields %in% colnames(corrected_scores)) &&
        all(
          corrected_scores$infercnv_observation_group ==
            config$cnv$infercnv_observation_group
        ) &&
        all(
          tolower(as.character(corrected_scores$infercnv_HMM)) ==
            tolower(as.character(config$cnv$infercnv_hmm))
        ) &&
        all(corrected_scores$infercnv_HMM_type == config$cnv$infercnv_hmm_type) &&
        all(corrected_scores$infercnv_analysis_mode == config$cnv$infercnv_analysis_mode) &&
        all(
          corrected_scores$infercnv_partition_method ==
            config$cnv$infercnv_tumor_subcluster_partition_method
        ) &&
        all(as.integer(corrected_scores$infercnv_k_nn) == config$cnv$infercnv_k_nn) &&
        all(
          as.character(corrected_scores$infercnv_leiden_resolution) ==
            as.character(config$cnv$infercnv_leiden_resolution)
        ) &&
        all(corrected_scores$infercnv_reference_policy == infercnv_reference_policy) &&
        all(
          corrected_scores$infercnv_reference_source ==
            reference_summary_all$reference_source[
              match(corrected_scores$sample_id, reference_summary_all$sample_id)
            ]
        )
    )
    if (!corrected_provenance_ok) {
      corrected_scores <- NULL
      f1_append_log(
        config,
        stage,
        "既有corrected inferCNV结果与当前subcluster/reference参数不一致，禁止复用。"
      )
    }
  }
  completed_samples <- if (is.null(corrected_scores)) character() else unique(corrected_scores$sample_id)
  new_rows <- list()
  for (sample_id in setdiff(corrected_samples, completed_samples)) {
    obs_cells <- candidate_cells[
      epithelial$sample_id[match(candidate_cells, colnames(epithelial))] %in% sample_id
    ]
    obs_cluster <- setNames(
      as.character(epithelial$epithelial_cluster_id[match(obs_cells, colnames(epithelial))]),
      obs_cells
    )
    reference <- choose_reference_cells(
      sample_id,
      all_cells,
      config$cnv$minimum_reference_cells,
      config$cnv$maximum_reference_cells
    )
    reference_cells <- setdiff(intersect(reference$cells, colnames(full_counts)), obs_cells)
    if (
      !identical(reference$status, "evaluable") ||
        length(reference_cells) < config$cnv$minimum_reference_cells
    ) {
      stop(sample_id, "当前reference不足，不能运行corrected inferCNV敏感性。")
    }
    input_cells <- c(obs_cells, reference_cells)
    corrected_counts <- load_corrected_infercnv_counts(
      input_cells,
      all_cells,
      matched_gene_order$gene_symbol,
      ambient_summary
    )
    corrected_dir <- file.path(cnv_dir, sample_id, "infercnv_decontX_corrected")
    corrected_run <- run_infercnv_matrix(
      counts = corrected_counts,
      input_cells = input_cells,
      obs_cells = obs_cells,
      obs_cluster = obs_cluster,
      reference_cells = reference_cells,
      sample_id = sample_id,
      gene_order = matched_gene_order,
      out_dir = corrected_dir,
      file_prefix = paste0(sample_id, "_infercnv_decontX_corrected"),
      config = config
    )
    raw_membership_index <- match(
      corrected_run$cell_id_final,
      infer_scores_all$cell_id_final
    )
    if (anyNA(raw_membership_index)) {
      stop(sample_id, " corrected inferCNV细胞无法映射raw subcluster。")
    }
    corrected_run$corrected_infercnv_subcluster_id <-
      corrected_run$infercnv_subcluster_id
    corrected_run$raw_infercnv_subcluster_id <-
      infer_scores_all$infercnv_subcluster_id[raw_membership_index]
    corrected_run$infercnv_observation_group <- config$cnv$infercnv_observation_group
    corrected_run$infercnv_HMM <- config$cnv$infercnv_hmm
    corrected_run$infercnv_HMM_type <- config$cnv$infercnv_hmm_type
    corrected_run$infercnv_analysis_mode <- config$cnv$infercnv_analysis_mode
    corrected_run$infercnv_partition_method <-
      config$cnv$infercnv_tumor_subcluster_partition_method
    corrected_run$infercnv_k_nn <- config$cnv$infercnv_k_nn
    corrected_run$infercnv_leiden_resolution <-
      as.character(config$cnv$infercnv_leiden_resolution)
    corrected_run$infercnv_reference_policy <- infercnv_reference_policy
    corrected_run$infercnv_reference_source <- reference$source
    new_rows[[sample_id]] <- corrected_run
    rm(corrected_counts)
    gc(verbose = FALSE)
  }
  if (length(new_rows)) {
    corrected_scores <- rbind(
      if (is.null(corrected_scores)) NULL else corrected_scores[
        !corrected_scores$sample_id %in% names(new_rows),
        ,
        drop = FALSE
      ],
      do.call(rbind, new_rows)
    )
    f1_write_tsv(corrected_scores, corrected_scores_path)
  }
}

corrected_support <- as.character(approved$corrected_infercnv_broad_segment_support)
requested_rows <- approved$run_corrected_infercnv_sensitivity
if (any(requested_rows & !corrected_support %in% allowed_infer)) {
  approved_index <- match(
    paste(
      review_template$sample_id,
      review_template$infercnv_subcluster_id,
      sep = "__"
    ),
    approved$key
  )
  for (field in review_fields) review_template[[field]] <- approved[[field]][approved_index]
  f1_write_tsv(review_template, config$paths$malignancy_review_template)
  stop("corrected inferCNV已生成；请审核corrected热图并填写其大片段支持等级后重跑F1.6。")
}
if (any(!requested_rows & corrected_support != "not_run")) {
  stop("未运行corrected敏感性的行必须把corrected_infercnv_broad_segment_support写为not_run。")
}
comparison <- data.frame(
  sample_id = approved$sample_id,
  infercnv_subcluster_id = approved$infercnv_subcluster_id,
  epithelial_cluster = approved$epithelial_cluster,
  corrected_sensitivity_run = requested_rows,
  raw_broad_segment_support = approved$infercnv_broad_segment_support,
  corrected_broad_segment_support = corrected_support,
  material_change = requested_rows &
    approved$infercnv_broad_segment_support != corrected_support,
  stringsAsFactors = FALSE
)
comparison_path <- file.path(
  config$paths$malignancy_dir,
  "infercnv_ambient_sensitivity_comparison.tsv"
)
f1_write_tsv(comparison, comparison_path)
if (any(comparison$material_change)) {
  stop(
    "corrected inferCNV改变了大片段支持等级，因此可能改变恶性标签或06a纳入。",
    "已按方案暂停并写出：", comparison_path
  )
}

infer_index <- match(colnames(epithelial), infer_scores_all$cell_id_final)
epithelial$infercnv_subcluster_id <-
  infer_scores_all$infercnv_subcluster_id[infer_index]
epithelial$infercnv_subcluster_local <-
  infer_scores_all$infercnv_subcluster_local[infer_index]
epi_key <- paste(
  epithelial$sample_id,
  epithelial$infercnv_subcluster_id,
  sep = "__"
)
review_index <- match(epi_key, approved$key)
copy_index <- match(colnames(epithelial), copykat_calls_all$cell_id_final)
epithelial$infercnv_broad_segment_support <- approved$infercnv_broad_segment_support[review_index]
epithelial$corrected_infercnv_broad_segment_support <-
  approved$corrected_infercnv_broad_segment_support[review_index]
epithelial$tumor_program_support <- approved$tumor_program_support[review_index]
epithelial$normal_program_support <- approved$normal_program_support[review_index]
epithelial$copykat_primary_raw_call <- copykat_calls_all$copykat_primary_raw_call[copy_index]
epithelial$copykat_call <- copykat_calls_all$copykat_call[copy_index]
epithelial$copykat_primary_baseline_status <-
  copykat_calls_all$copykat_primary_baseline_status[copy_index]
epithelial$copykat_primary_evidence_usable <-
  copykat_calls_all$copykat_primary_evidence_usable[copy_index]
epithelial$copykat_B_immune_call <- copykat_calls_all$copykat_B_immune_call[copy_index]
epithelial$copykat_C_holdout_call <- copykat_calls_all$copykat_C_holdout_call[copy_index]
epithelial$infercnv_cell_burden <- infer_scores_all$infercnv_cell_burden[match(colnames(epithelial), infer_scores_all$cell_id_final)]

cluster_approved <- epithelial$epithelial_cluster_id %in% approved_clusters
infer_strong <- epithelial$infercnv_broad_segment_support == "strong"
infer_negative_or_weak <- epithelial$infercnv_broad_segment_support %in% c("weak", "absent")
infer_not_strong <- epithelial$infercnv_broad_segment_support %in%
  c("weak", "absent", "not_evaluable")
aneuploid <- epithelial$copykat_call == "aneuploid"
diploid <- epithelial$copykat_call == "diploid"
uncalled <- is.na(epithelial$copykat_call) | epithelial$copykat_call == "uncalled"
tumor_source <- epithelial$group_analysis %in% c("Primary_Tumor", "Peritoneal_Metastasis")
tumor_context <- tumor_source & epithelial$tumor_program_support
normal_context <- epithelial$normal_program_support

label <- rep("epithelial_uncertain", ncol(epithelial))
label[!cluster_approved] <- "exclude_contamination_or_doublet"
label[
  cluster_approved & infer_negative_or_weak & diploid &
    normal_context & !epithelial$tumor_program_support
] <- "non_malignant_epithelial"
label[
  cluster_approved & aneuploid & infer_not_strong &
    tumor_context & !normal_context
] <- "malignant_probable_copykat"
label[
  cluster_approved & infer_strong & (diploid | uncalled) &
    tumor_context & !normal_context
] <- "malignant_probable_infercnv"
  label[
    cluster_approved & infer_strong & aneuploid &
      tumor_context & !normal_context
  ] <- "malignant_high_confidence"
  epithelial$malignancy_label <- label
  epithelial$include_in_06a <- label %in% c(
    "malignant_high_confidence", "malignant_probable_infercnv"
  )
  epithelial$include_in_06b <- label == "malignant_high_confidence"
  if (any(epithelial$include_in_06a & !tumor_source)) {
    stop("正常胃或其他非肿瘤来源上皮不得进入06a/06b。")
  }
  if (any(
    epithelial$include_in_06a &
      !epithelial$infercnv_broad_segment_support %in% c("strong", "weak")
  )) {
    stop("06a每个细胞必须至少具有weak inferCNV热图大片段支持。")
  }
  if (any(epithelial$malignancy_label == "malignant_probable_copykat" & epithelial$include_in_06a)) {
    stop("CopyKAT单独支持的细胞不得进入06a。")
  }

calling <- data.frame(
  cell_id_final = colnames(epithelial),
  sample_id = epithelial$sample_id,
  patient_id = epithelial$patient_id,
  group_analysis = epithelial$group_analysis,
  epithelial_cluster = epithelial$epithelial_cluster_id,
  infercnv_subcluster_id = epithelial$infercnv_subcluster_id,
  infercnv_cell_burden = epithelial$infercnv_cell_burden,
  infercnv_broad_segment_support = epithelial$infercnv_broad_segment_support,
  corrected_infercnv_broad_segment_support =
    epithelial$corrected_infercnv_broad_segment_support,
  copykat_primary_raw_call = epithelial$copykat_primary_raw_call,
  copykat_call = epithelial$copykat_call,
  copykat_primary_baseline_status = epithelial$copykat_primary_baseline_status,
  copykat_primary_evidence_usable = epithelial$copykat_primary_evidence_usable,
  copykat_B_immune_call = epithelial$copykat_B_immune_call,
  copykat_C_holdout_call = epithelial$copykat_C_holdout_call,
  tumor_program_support = epithelial$tumor_program_support,
  normal_program_support = epithelial$normal_program_support,
  malignancy_label = epithelial$malignancy_label,
  include_in_06a = epithelial$include_in_06a,
  include_in_06b = epithelial$include_in_06b,
  stringsAsFactors = FALSE
)
f1_write_tsv(calling, file.path(config$paths$malignancy_dir, "malignant_cell_calling_summary.tsv"))
f1_write_tsv(epithelial_review, file.path(config$paths$annotation_dir, "epithelial_contamination_review.tsv"))
f1_save_rds_atomic(epithelial, config$paths$object_05, compress = FALSE)

n_main <- sum(epithelial$include_in_06a)
n_high <- sum(epithelial$include_in_06b)
if (n_main == 0L || n_high == 0L) {
  stop(
    "联合判定后06a或06b没有细胞（06a=", n_main, ", 06b=", n_high,
    "）。05对象和摘要已保存，请检查reference与审核证据，不得放宽规则凑数。"
  )
}
main_object <- subset(epithelial, cells = colnames(epithelial)[epithelial$include_in_06a])
high_object <- subset(epithelial, cells = colnames(epithelial)[epithelial$include_in_06b])
f1_save_rds_atomic(main_object, config$paths$object_06a, compress = FALSE)
f1_save_rds_atomic(high_object, config$paths$object_06b, compress = FALSE)

# 汇总全部公开细胞的关键metadata，方便后续按原始细胞追溯排除原因。
initial <- readRDS(config$paths$object_01)
final_meta <- initial[[]]
full_meta <- all_cells[[]]
for (field in c("seurat_clusters", "cell_type_major", "cell_type_minor", "cell_state", "annotation_confidence")) {
  final_meta[[field]] <- full_meta[[field]][match(rownames(final_meta), rownames(full_meta))]
}
for (field in c(
  "epithelial_cluster_id", "infercnv_subcluster_id",
  "infercnv_subcluster_local", "copykat_primary_raw_call", "copykat_call",
  "copykat_primary_baseline_status", "copykat_primary_evidence_usable",
  "copykat_B_immune_call", "copykat_C_holdout_call", "infercnv_cell_burden",
  "infercnv_broad_segment_support", "corrected_infercnv_broad_segment_support",
  "tumor_program_support", "normal_program_support",
  "malignancy_label", "include_in_06a", "include_in_06b"
)) {
  final_meta[[field]] <- epithelial[[field, drop = TRUE]][match(rownames(final_meta), colnames(epithelial))]
}
final_meta$cell_id_final <- rownames(final_meta)
f1_write_tsv(final_meta, file.path(config$paths$malignancy_dir, "cell_metadata_final.tsv"))

label_counts <- table(epithelial$malignancy_label)
n_infercnv_not_evaluable <- sum(reference_summary_all$infercnv_evaluation_status != "evaluable")
status <- "PASS_WITH_NOTED_LIMITATIONS"
report <- c(
  "# F1 Final Report", "",
  paste0("Generated at: ", f1_now()), "",
  "## Main Result", "",
  paste0("- F1 status: ", status),
  paste0("- QC/doublet-filtered all cells: ", ncol(all_cells)),
  paste0("- Reclustered epithelial cells: ", ncol(epithelial)),
  paste0("- 06a malignant main cells: ", n_main),
  paste0("- 06b high-confidence malignant cells: ", n_high),
  paste0("- inferCNV-not-evaluable samples: ", n_infercnv_not_evaluable),
  paste0("- Corrected inferCNV sensitivity samples: ", length(corrected_samples)),
  paste0("- Malignancy labels: ", paste0(names(label_counts), "=", as.integer(label_counts), collapse = "; ")), "",
  "## Methods", "",
  "- Fixed per-sample min.cells=3 feature space and five frozen cell thresholds were used.",
  "- scDblFinder was the deletion rule; DoubletFinder was sensitivity-only.",
  "- Per-sample SCTransform v2 used vars.to.regress=NULL; Harmony corrected sample_id only.",
  "- UCell/MLMOD and prognosis were not used in F1 decisions.",
  paste0(
    "- inferCNV used analysis_mode=subclusters with Leiden, k_nn=",
    config$cnv$infercnv_k_nn,
    ", resolution=", config$cnv$infercnv_leiden_resolution,
    ", one observations group per sample and formal i6 HMM; heatmap broad segments remained primary."
  ),
  paste0(
    "- inferCNV references followed ", infercnv_reference_policy,
    "; CopyKAT A/B/C inputs and known normals remained current-sample-only."
  ),
  "- CopyKAT A was self-estimated baseline; B used same-sample immune known normals; C was Normal_Gastric two-fold holdout sensitivity.",
  "- CopyKAT-only probable cells were retained in object 05 but excluded from 06a and 06b.",
  "- inferCNV and CopyKAT used RNA raw counts and were interpreted jointly with approved biological review.",
  "- Epithelial lineage markers were compared descriptively with same-lineage Normal_Gastric cells; no universal cancer-marker score was used.",
  "- Corrected inferCNV was conditional sensitivity only; raw counts remained primary.", "",
  "## Limitations", "",
  "- Raw/empty droplets are unavailable; cell calling, emptyDrops, SoupX and CellBender were not evaluable.",
  "- inferCNV is now deprecated upstream; the exact package version and full outputs were retained for reproducibility.",
  "- inferCNV and CopyKAT are both RNA-derived and do not constitute independent DNA validation.",
  "- inferCNV Leiden subclusters are CNV-pattern groups, not proven tumor clones.",
  "- Paired or pooled Normal_Gastric inferCNV references remain susceptible to cross-sample technical differences.",
  "- CopyKAT diploid calls cannot exclude near-diploid malignant cells.",
  "- The public matrix and repeated mt_percent<=20 rule limit conclusions about cells above that boundary.", "",
  "## F2 Entry", "",
  "- 06a is the F2 main object; 06b is the high-confidence sensitivity object.",
  "- Enter F2 only after independent script/result review and user approval."
)
writeLines(report, file.path(config$paths$malignancy_dir, "F1_final_report.md"), useBytes = TRUE)
f1_write_parameter_versions(config)
f1_save_session_info(config, "F1_06_malignancy_inference")
f1_append_log(config, stage, sprintf("完成F1：06a=%d，06b=%d；状态=%s", n_main, n_high, status))
message("F1.6完成：", config$paths$object_06a)
