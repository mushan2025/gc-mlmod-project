# F1.4 主要细胞类型注释与DecontX评估 ---------------------------------------
#
# 生物学目的：依据cluster marker、冻结marker panel和样本覆盖，为每个cluster
# 赋予可解释的主要谱系；批准粗谱系后，再按样本使用已经审核和注释的
# Seurat cluster作为DecontX模型分组，估计retained-cell环境RNA污染。
# 粗谱系只用于结果汇总，最终标签必须经研究者审核，DecontX不反向改写标签。
#
# 第一次正式运行会生成注释模板并停止；填好批准表后再次运行才执行DecontX并
# 保存正式注释对象。raw counts始终是主矩阵，corrected counts只供敏感性检查。

file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), winslash = "/"))
source(file.path(script_dir, "F1_config.R"), encoding = "UTF-8")
source(file.path(script_dir, "F1_utils.R"), encoding = "UTF-8")

args <- f1_parse_args()
config <- f1_build_config(args$project_root)
stage <- "F1.4"
marker_path <- file.path(config$paths$annotation_dir, "cluster_marker_genes.tsv")
outputs <- c(
  config$paths$annotation_template,
  marker_path,
  config$paths$object_03,
  file.path(config$paths$annotation_dir, "celltype_annotation_summary.tsv"),
  file.path(config$paths$annotation_dir, "cluster_quality_summary.tsv"),
  config$paths$ambient_cell_estimates,
  config$paths$ambient_summary,
  config$paths$ambient_lineage_summary,
  config$paths$ambient_cluster_summary,
  config$paths$ambient_corrected_count_summary
)

if (!args$execute) {
  f1_stage_dry_run(
    stage,
    paste0(
      "完成主要谱系人工注释后，按样本以已审核Seurat cluster作为模型分组",
      "运行DecontX；粗谱系只用于汇总解释"
    ),
    c(config$paths$object_03a, config$paths$marker_panel, config$paths$annotation_approved),
    outputs,
    config$packages[[stage]]
  )
  quit(save = "no", status = 0)
}

f1_require_packages(config$packages[[stage]], stage)
f1_check_f0_ready(config)
if (!file.exists(config$paths$object_03a)) stop("请先完成F1.3：", config$paths$object_03a)
f1_prepare_directories(config)
set.seed(config$seed)
f1_append_log(config, stage, "开始cluster marker计算、主要谱系注释与注释后DecontX评估")

run_decontx_after_annotation <- function(counts, model_clusters, sample_id, config) {
  reliable_model_clusters <- unique(model_clusters)
  if (length(reliable_model_clusters) < config$ambient$minimum_reliable_lineages) {
    return(list(
      status = "not_evaluable_fewer_than_two_model_clusters",
      contamination = setNames(rep(NA_real_, ncol(counts)), colnames(counts)),
      corrected = NULL,
      reliable_model_clusters = reliable_model_clusters,
      error = NA_character_
    ))
  }

  tryCatch({
    set.seed(config$seed)
    sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = counts))
    result <- celda::decontX(
      sce,
      z = model_clusters,
      background = NULL,
      seed = config$seed
    )
    cell_data <- as.data.frame(SummarizedExperiment::colData(result))
    if (!"decontX_contamination" %in% colnames(cell_data)) {
      stop("DecontX结果缺少decontX_contamination。")
    }
    if (!setequal(rownames(cell_data), colnames(counts))) {
      stop("DecontX返回的cell ID与输入不一致。")
    }
    corrected <- celda::decontXcounts(result)
    rownames(corrected) <- rownames(counts)
    colnames(corrected) <- colnames(counts)
    list(
      status = "completed_with_researcher_approved_seurat_clusters",
      contamination = setNames(
        as.numeric(cell_data[colnames(counts), "decontX_contamination"]),
        colnames(counts)
      ),
      corrected = corrected,
      reliable_model_clusters = reliable_model_clusters,
      error = NA_character_
    )
  }, error = function(e) {
    list(
      status = "not_evaluable_model_failure",
      contamination = setNames(rep(NA_real_, ncol(counts)), colnames(counts)),
      corrected = NULL,
      reliable_model_clusters = reliable_model_clusters,
      error = conditionMessage(e)
    )
  })
}

object <- readRDS(config$paths$object_03a)
f1_require_columns(
  object[[]],
  c(
    "sample_id", "patient_id", "group_analysis", "seurat_clusters",
    "nCount_RNA", "nFeature_RNA", "mt_percent", "HB_percent",
    "DoubletFinder_class", "scDblFinder_class"
  ),
  "F1.3 clustered object metadata"
)
object <- f1_join_assay(object, "RNA")
SeuratObject::DefaultAssay(object) <- "RNA"
SeuratObject::Idents(object) <- "seurat_clusters"
clusters <- sort(unique(as.character(object$seurat_clusters)))

# 已有marker表且覆盖当前全部cluster时直接复用，避免人工审核后重复长时间计算。
reuse_markers <- FALSE
if (file.exists(marker_path)) {
  candidate <- f1_read_tsv(marker_path)
  if (all(c("cluster", "gene") %in% colnames(candidate)) && setequal(unique(as.character(candidate$cluster)), clusters)) {
    markers <- candidate
    reuse_markers <- TRUE
  }
}
if (!reuse_markers) {
  markers <- Seurat::FindAllMarkers(
    object,
    assay = "RNA",
    only.pos = config$markers$only_positive,
    min.pct = config$markers$min_pct,
    logfc.threshold = config$markers$logfc_threshold,
    test.use = "wilcox",
    verbose = TRUE
  )
  if (!"gene" %in% colnames(markers)) markers$gene <- rownames(markers)
  f1_write_tsv(markers, marker_path)
}

top <- f1_top_markers(markers, n = config$markers$top_n_per_cluster)
composition <- f1_cluster_composition(object[[]])
template <- composition
template <- merge(template, top, by = "cluster", all.x = TRUE, sort = TRUE)
template$cell_type_major <- NA_character_
template$cell_type_minor <- NA_character_
template$cell_state <- NA_character_
template$annotation_confidence <- NA_character_
template$annotation_reason <- NA_character_
template$annotation_review_status <- "pending_full_F1.4_review"
template$downstream_handling_before_full_approval <- NA_character_
template <- f1_apply_f1_annotation_draft(template)
template <- template[
  match(f1_sort_cluster_ids(template$cluster), as.character(template$cluster)),
  ,
  drop = FALSE
]
f1_write_tsv(template, config$paths$annotation_template)

panel <- f1_read_tsv(config$paths$marker_panel)
f1_require_columns(panel, c("cell_type", "positive_markers"), "cell_type_marker_panel.tsv")
broad_rows <- panel$cell_type %in% c(
  "Epithelial", "T_cell", "NK_cell", "B_cell", "Plasma_cell", "Myeloid",
  "Mast_cell", "Fibroblast", "Endothelial", "Pericyte"
)
broad_features <- unique(unlist(strsplit(panel$positive_markers[broad_rows], ",", fixed = TRUE)))
broad_features <- intersect(trimws(broad_features), rownames(object))
if (length(broad_features)) {
  dot <- Seurat::DotPlot(object, features = broad_features, assay = "RNA", group.by = "seurat_clusters") +
    ggplot2::coord_flip() +
    ggplot2::labs(title = "Broad-lineage marker review", x = NULL, y = "cluster") +
    ggplot2::theme_bw(base_size = 9)
  ggplot2::ggsave(
    file.path(config$paths$annotation_dir, "F1_broad_lineage_marker_dotplot.pdf"),
    dot,
    width = 12,
    height = 10,
    limitsize = FALSE
  )
}

if (!file.exists(config$paths$annotation_approved)) {
  f1_append_log(config, stage, paste0("已生成注释模板，等待研究者批准：", config$paths$annotation_template))
  stop(
    "F1.4已生成marker和注释模板。请审核后建立：",
    config$paths$annotation_approved,
    "；保留模板字段并为每个cluster填写major/minor/state/confidence/reason。"
  )
}

approved <- f1_read_tsv(config$paths$annotation_approved)
required_fields <- c(
  "cluster", "cell_type_major", "cell_type_minor", "cell_state",
  "annotation_confidence", "annotation_reason"
)
f1_require_columns(approved, required_fields, "F1_cluster_annotation_approved.tsv")
if (anyDuplicated(as.character(approved$cluster)) || !setequal(as.character(approved$cluster), clusters)) {
  stop("批准注释表必须与当前cluster一对一完整对应。")
}
text_fields <- setdiff(required_fields, "cluster")
incomplete <- vapply(text_fields, function(field) {
  any(is.na(approved[[field]]) | !nzchar(trimws(as.character(approved[[field]]))))
}, logical(1))
if (any(incomplete)) {
  stop(
    "批准注释表仍有未填写字段：",
    paste(text_fields[incomplete], collapse = ", ")
  )
}
allowed_major <- c(
  "Epithelial", "T/NK", "B/Plasma", "Myeloid", "Fibroblast/CAF",
  "Endothelial/Pericyte", "Mast", "Mesothelial",
  "Uncertain", "Mixed_or_doublet_suspect"
)
if (any(!approved$cell_type_major %in% allowed_major)) {
  stop("批准表含未登记的cell_type_major：", paste(unique(approved$cell_type_major[!approved$cell_type_major %in% allowed_major]), collapse = ", "))
}
if (any(!approved$annotation_confidence %in% c("high", "medium", "low"))) {
  stop("annotation_confidence只允许high、medium或low。")
}
approved_annotations <- approved[, required_fields, drop = FALSE]

map_index <- match(as.character(object$seurat_clusters), as.character(approved_annotations$cluster))
for (field in setdiff(required_fields, "cluster")) {
  object[[field]] <- as.character(approved_annotations[[field]][map_index])
}

template_context <- template[
  ,
  setdiff(
    colnames(template),
    c(
      "cell_type_major", "cell_type_minor", "cell_state",
      "annotation_confidence", "annotation_reason"
    )
  ),
  drop = FALSE
]
annotation_summary <- merge(
  template_context,
  approved_annotations,
  by = "cluster",
  sort = TRUE
)
annotation_summary$annotation_review_status <- "researcher_approved_annotation"
annotation_summary <- annotation_summary[
  match(
    f1_sort_cluster_ids(annotation_summary$cluster),
    as.character(annotation_summary$cluster)
  ),
  ,
  drop = FALSE
]
f1_write_tsv(annotation_summary, file.path(config$paths$annotation_dir, "celltype_annotation_summary.tsv"))

meta <- object[[]]
cluster_qc <- f1_cluster_qc_metrics(meta)
cluster_quality <- merge(composition, top, by = "cluster", all.x = TRUE, sort = TRUE)
cluster_quality <- merge(cluster_quality, cluster_qc, by = "cluster", all.x = TRUE, sort = TRUE)
cluster_quality <- merge(
  cluster_quality,
  approved_annotations,
  by = "cluster",
  all.x = TRUE,
  sort = TRUE
)
cluster_quality$review_status <- "researcher_approved_annotation"
cluster_quality <- cluster_quality[
  match(f1_sort_cluster_ids(cluster_quality$cluster), as.character(cluster_quality$cluster)),
  ,
  drop = FALSE
]
f1_write_tsv(cluster_quality, file.path(config$paths$annotation_dir, "cluster_quality_summary.tsv"))

p1 <- Seurat::DimPlot(object, reduction = "umap", group.by = "cell_type_major", label = TRUE, repel = TRUE) +
  ggplot2::ggtitle("F1 major cell lineages")
p2 <- Seurat::DimPlot(object, reduction = "umap", group.by = "sample_id") +
  ggplot2::ggtitle("F1 major cell lineages: sample source")
ggplot2::ggsave(
  file.path(config$paths$annotation_dir, "F1_all_cells_annotated_umap.pdf"),
  patchwork::wrap_plots(p1, p2, ncol = 2),
  width = 14,
  height = 6,
  limitsize = FALSE
)

# DecontX依赖细胞群标签。8个粗谱系会把普通B细胞和多种浆细胞状态合并，
# 在本数据中造成明显的伪高污染，因此模型z固定使用24个已审核Seurat cluster。
# 粗谱系只用于汇总；不使用恶性标签、MLMOD或任何下游结果。
raw_counts <- SeuratObject::LayerData(object, assay = "RNA", layer = "counts")
f1_assert_integer_counts(raw_counts, "F1.4 DecontX raw RNA counts")
sample_ids <- sort(unique(as.character(object$sample_id)))
dir.create(config$paths$decontx_corrected_dir, recursive = TRUE, showWarnings = FALSE)

ambient_cell_rows <- vector("list", length(sample_ids))
ambient_summary_rows <- vector("list", length(sample_ids))
names(ambient_cell_rows) <- sample_ids
names(ambient_summary_rows) <- sample_ids
ambient_contamination <- setNames(rep(NA_real_, ncol(object)), colnames(object))
ambient_status <- setNames(rep(NA_character_, ncol(object)), colnames(object))

for (sample_id in sample_ids) {
  sample_cells <- colnames(object)[as.character(object$sample_id) == sample_id]
  coarse_labels <- setNames(
    as.character(object$cell_type_major[match(sample_cells, colnames(object))]),
    sample_cells
  )
  model_clusters <- setNames(
    as.character(object$seurat_clusters[match(sample_cells, colnames(object))]),
    sample_cells
  )
  if (anyNA(coarse_labels) || any(!nzchar(coarse_labels))) {
    stop(sample_id, "存在缺失粗谱系标签，不能运行DecontX。")
  }
  if (anyNA(model_clusters) || any(!nzchar(model_clusters))) {
    stop(sample_id, "存在缺失Seurat cluster标签，不能运行DecontX。")
  }
  if (!all(unique(model_clusters) %in% as.character(approved_annotations$cluster))) {
    stop(sample_id, "包含未通过研究者注释审核的Seurat cluster。")
  }

  sample_counts <- raw_counts[, sample_cells, drop = FALSE]
  sample_counts <- sample_counts[Matrix::rowSums(sample_counts) > 0, , drop = FALSE]
  result <- run_decontx_after_annotation(
    counts = sample_counts,
    model_clusters = model_clusters,
    sample_id = sample_id,
    config = config
  )
  ambient_contamination[names(result$contamination)] <- result$contamination
  ambient_status[sample_cells] <- result$status

  corrected_path <- NA_character_
  if (!is.null(result$corrected)) {
    corrected_path <- file.path(
      config$paths$decontx_corrected_dir,
      paste0(sample_id, "_decontX_corrected_counts.rds")
    )
    f1_save_rds_atomic(result$corrected, corrected_path, compress = FALSE)
  }

  finite_contamination <- result$contamination[is.finite(result$contamination)]
  contamination_stats <- if (length(finite_contamination)) {
    c(
      f1_summary_stats(finite_contamination, "contamination"),
      contamination_P90 = as.numeric(stats::quantile(
        finite_contamination,
        probs = 0.90,
        na.rm = TRUE,
        names = FALSE,
        type = 7
      )),
      contamination_P95 = as.numeric(stats::quantile(
        finite_contamination,
        probs = 0.95,
        na.rm = TRUE,
        names = FALSE,
        type = 7
      ))
    )
  } else {
    setNames(
      as.list(rep(NA_real_, 7)),
      c(
        "contamination_min", "contamination_Q1", "contamination_median",
        "contamination_Q3", "contamination_max", "contamination_P90",
        "contamination_P95"
      )
    )
  }
  coarse_counts <- table(coarse_labels)
  model_cluster_counts <- table(model_clusters)
  raw_total_counts <- sum(sample_counts)
  corrected_values <- if (!is.null(result$corrected)) {
    if (inherits(result$corrected, "sparseMatrix")) {
      result$corrected@x
    } else {
      as.vector(result$corrected)
    }
  } else {
    numeric()
  }
  corrected_total_counts <- if (length(corrected_values)) {
    sum(corrected_values)
  } else {
    NA_real_
  }
  ambient_summary_rows[[sample_id]] <- cbind(
    data.frame(
      sample_id = sample_id,
      input_cells = length(sample_cells),
      input_genes = nrow(sample_counts),
      input_matrix = "fixed_QC_scDblFinder_singlet_raw_integer_counts",
      method = "DecontX",
      celda_version = as.character(utils::packageVersion("celda")),
      status = result$status,
      cluster_label_source = "researcher_approved_seurat_cluster_partition",
      model_cluster_counts = paste0(
        names(model_cluster_counts), "=", as.integer(model_cluster_counts),
        collapse = "|"
      ),
      coarse_label_counts = paste0(
        names(coarse_counts), "=", as.integer(coarse_counts), collapse = "|"
      ),
      reliable_model_clusters = paste(
        f1_sort_cluster_ids(result$reliable_model_clusters),
        collapse = "|"
      ),
      background_input = "not_available_public_filtered_cells",
      corrected_counts_path = corrected_path,
      corrected_matrix_class = if (!is.null(result$corrected)) {
        paste(class(result$corrected), collapse = "|")
      } else {
        NA_character_
      },
      raw_total_counts = as.numeric(raw_total_counts),
      corrected_total_counts = as.numeric(corrected_total_counts),
      removed_count_fraction = if (
        is.finite(corrected_total_counts) && raw_total_counts > 0
      ) {
        as.numeric(
          (raw_total_counts - corrected_total_counts) / raw_total_counts
        )
      } else {
        NA_real_
      },
      corrected_noninteger_fraction = if (length(corrected_values)) {
        mean(abs(corrected_values - round(corrected_values)) > 1e-8)
      } else {
        NA_real_
      },
      raw_counts_remain_main = TRUE,
      deletion_based_on_contamination = FALSE,
      interpretation = "retained_cell_ambient_estimate_and_conditional_sensitivity_only",
      error = if (is.na(result$error)) {
        NA_character_
      } else {
        gsub("[\r\n\t]+", " ", result$error)
      },
      seed = config$seed,
      stringsAsFactors = FALSE
    ),
    as.data.frame(contamination_stats, check.names = FALSE)
  )
  ambient_cell_rows[[sample_id]] <- data.frame(
    cell_id_final = sample_cells,
    sample_id = sample_id,
    seurat_cluster = as.character(object$seurat_clusters[match(sample_cells, colnames(object))]),
    decontX_model_cluster = unname(model_clusters),
    coarse_lineage_label = unname(coarse_labels),
    cell_type_minor = as.character(
      object$cell_type_minor[match(sample_cells, colnames(object))]
    ),
    cell_state = as.character(
      object$cell_state[match(sample_cells, colnames(object))]
    ),
    model_label_source = "researcher_approved_seurat_cluster_partition",
    decontX_status = result$status,
    retained_cell_ambient_contamination_estimate =
      as.numeric(result$contamination[sample_cells]),
    stringsAsFactors = FALSE
  )

  rm(sample_counts, result, corrected_values)
  gc(verbose = FALSE)
}

ambient_cells <- do.call(rbind, ambient_cell_rows)
ambient_summary <- do.call(rbind, ambient_summary_rows)
object <- Seurat::AddMetaData(
  object,
  metadata = data.frame(
    retained_cell_ambient_contamination_estimate =
      as.numeric(ambient_contamination[colnames(object)]),
    decontX_evaluation_status = as.character(ambient_status[colnames(object)]),
    row.names = colnames(object),
    check.names = FALSE
  )
)
f1_write_tsv(ambient_cells, config$paths$ambient_cell_estimates)
f1_write_tsv(ambient_summary, config$paths$ambient_summary)

ambient_dt <- data.table::as.data.table(ambient_cells)
lineage_summary <- ambient_dt[
  is.finite(retained_cell_ambient_contamination_estimate),
  .(
    n_cells = .N,
    n_samples = data.table::uniqueN(sample_id),
    contamination_mean = mean(
      retained_cell_ambient_contamination_estimate
    ),
    contamination_median = stats::median(
      retained_cell_ambient_contamination_estimate
    ),
    contamination_P90 = as.numeric(stats::quantile(
      retained_cell_ambient_contamination_estimate,
      probs = 0.90,
      names = FALSE,
      type = 7
    )),
    contamination_P95 = as.numeric(stats::quantile(
      retained_cell_ambient_contamination_estimate,
      probs = 0.95,
      names = FALSE,
      type = 7
    )),
    fraction_ge_0_10 = mean(
      retained_cell_ambient_contamination_estimate >= 0.10
    ),
    fraction_ge_0_20 = mean(
      retained_cell_ambient_contamination_estimate >= 0.20
    )
  ),
  by = coarse_lineage_label
][order(-n_cells)]
f1_write_tsv(lineage_summary, config$paths$ambient_lineage_summary)

cluster_summary <- ambient_dt[
  is.finite(retained_cell_ambient_contamination_estimate),
  .(
    n_cells = .N,
    n_samples = data.table::uniqueN(sample_id),
    contamination_mean = mean(
      retained_cell_ambient_contamination_estimate
    ),
    contamination_median = stats::median(
      retained_cell_ambient_contamination_estimate
    ),
    contamination_P90 = as.numeric(stats::quantile(
      retained_cell_ambient_contamination_estimate,
      probs = 0.90,
      names = FALSE,
      type = 7
    )),
    contamination_P95 = as.numeric(stats::quantile(
      retained_cell_ambient_contamination_estimate,
      probs = 0.95,
      names = FALSE,
      type = 7
    )),
    fraction_ge_0_10 = mean(
      retained_cell_ambient_contamination_estimate >= 0.10
    ),
    fraction_ge_0_20 = mean(
      retained_cell_ambient_contamination_estimate >= 0.20
    )
  ),
  by = .(
    seurat_cluster, coarse_lineage_label, cell_type_minor, cell_state
  )
][order(as.integer(seurat_cluster))]
f1_write_tsv(cluster_summary, config$paths$ambient_cluster_summary)

corrected_count_summary <- ambient_summary[
  ,
  c(
    "sample_id", "input_cells", "input_genes", "status",
    "corrected_counts_path", "corrected_matrix_class",
    "raw_total_counts", "corrected_total_counts",
    "removed_count_fraction", "corrected_noninteger_fraction"
  ),
  drop = FALSE
]
f1_write_tsv(
  corrected_count_summary,
  config$paths$ambient_corrected_count_summary
)

plot_data <- ambient_cells[
  is.finite(ambient_cells$retained_cell_ambient_contamination_estimate),
  ,
  drop = FALSE
]
if (nrow(plot_data)) {
  ambient_plot <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = sample_id,
      y = retained_cell_ambient_contamination_estimate
    )
  ) +
    ggplot2::geom_boxplot(outlier.size = 0.2, linewidth = 0.25) +
    ggplot2::coord_flip() +
    ggplot2::labs(
      title = "逐样本DecontX污染估计（已审核cluster作为模型分组）",
      x = "sample_id",
      y = "DecontX contamination"
    ) +
    ggplot2::theme_bw(base_size = 9)
  ggplot2::ggsave(
    file.path(config$paths$qc_dir, "F1_DecontX_contamination_by_sample.pdf"),
    ambient_plot,
    width = 8,
    height = max(8, 0.24 * length(sample_ids)),
    limitsize = FALSE
  )
}

f1_save_rds_atomic(object, config$paths$object_03, compress = FALSE)
f1_save_session_info(config, "F1_04_annotation_and_decontx")
f1_append_log(
  config,
  stage,
  sprintf(
    paste0(
      "完成研究者批准的主要谱系注释（%d个cluster）及逐样本DecontX；",
      "DecontX z使用已审核Seurat cluster，raw counts保持主矩阵"
    ),
    length(clusters)
  )
)
message("F1.4完成：", config$paths$object_03)
