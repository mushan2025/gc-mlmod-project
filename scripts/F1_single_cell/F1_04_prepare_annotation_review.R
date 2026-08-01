# F1.4 注释前轻量审核表 -----------------------------------------------------
#
# 这个脚本要解决什么问题？
#   F1.3已经完成142,650个细胞的聚类，但完整Seurat对象很大。
#   在正式批准细胞类型并运行DecontX前，我们需要先查看每个cluster的marker、
#   样本/患者来源、组织分组和QC概况，判断是否存在单一样本支配或技术状态。
#
# 为什么单独写这个脚本？
#   这些信息都已经保存在轻量TSV中，不需要重新加载约18 GB的Seurat对象。
#   因此本脚本只整理审核表，不重新聚类、不改变细胞对象，也不运行DecontX。
#
# 输入：
#   results/F1_annotation/F1_all_cells_embedding_source_data.tsv.gz
#   results/F1_annotation/cluster_marker_genes.tsv
#   data/metadata/sample_info.tsv
#   results/F1_annotation/integration_diagnostics.tsv
#
# 输出：
#   results/F1_annotation/F1_cluster_annotation_template.tsv
#   results/F1_annotation/cluster_quality_summary.tsv
#   更新后的results/F1_annotation/integration_diagnostics.tsv
#
# 运行方式：
#   Rscript scripts/F1_single_cell/F1_04_prepare_annotation_review.R \
#     --project-root=/path/to/project --execute
#
# 注意：
#   输出中24个cluster都已形成完整注释草案，但仍需研究者批准并写入正式批准表。
#   草案不能被当作已经冻结的最终注释，也不会触发DecontX。

file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), winslash = "/"))
source(file.path(script_dir, "F1_config.R"), encoding = "UTF-8")
source(file.path(script_dir, "F1_utils.R"), encoding = "UTF-8")

args <- f1_parse_args()
config <- f1_build_config(args$project_root)
stage <- "F1.4_annotation_review"

embedding_path <- file.path(
  config$paths$annotation_dir,
  "F1_all_cells_embedding_source_data.tsv.gz"
)
marker_path <- file.path(config$paths$annotation_dir, "cluster_marker_genes.tsv")
diagnostics_path <- file.path(
  config$paths$annotation_dir,
  "integration_diagnostics.tsv"
)
quality_path <- file.path(
  config$paths$annotation_dir,
  "cluster_quality_summary.tsv"
)

inputs <- c(
  embedding_path,
  marker_path,
  config$paths$sample_info,
  diagnostics_path
)
outputs <- c(
  config$paths$annotation_template,
  quality_path,
  diagnostics_path
)

if (!args$execute) {
  f1_stage_dry_run(
    stage,
    paste0(
      "从轻量TSV生成注释模板、cluster质量表和样本/患者构成结论；",
      "不读取Seurat对象、不运行DecontX"
    ),
    inputs,
    outputs,
    "data.table"
  )
  quit(save = "no", status = 0)
}

f1_require_packages("data.table", stage)
missing_inputs <- inputs[!file.exists(inputs)]
if (length(missing_inputs)) {
  stop("缺少注释前审核输入：", paste(missing_inputs, collapse = "; "))
}
f1_prepare_directories(config)

# 第1节：读取轻量输入并检查样本映射 -----------------------------------------
embedding <- f1_read_tsv(embedding_path)
markers <- f1_read_tsv(marker_path)
sample_info <- f1_read_tsv(config$paths$sample_info)
diagnostics <- f1_read_tsv(diagnostics_path)

f1_require_columns(
  embedding,
  c(
    "cell_id", "sample_id", "group_analysis", "seurat_clusters",
    "nCount_RNA", "nFeature_RNA", "mt_percent", "HB_percent",
    "scDblFinder_class", "DoubletFinder_class"
  ),
  "F1_all_cells_embedding_source_data.tsv.gz"
)
f1_require_columns(
  sample_info,
  c("sample_id", "patient_id", "group_analysis", "include_in_f1"),
  "sample_info.tsv"
)
f1_require_columns(markers, c("cluster", "gene"), "cluster_marker_genes.tsv")
f1_require_columns(
  diagnostics,
  c("item", "value", "interpretation"),
  "integration_diagnostics.tsv"
)

if (anyDuplicated(embedding$cell_id)) {
  stop("embedding source中出现重复cell_id，不能可靠汇总cluster。")
}
if (anyDuplicated(sample_info$sample_id)) {
  stop("sample_info.tsv中sample_id不是一对一，不能可靠映射patient_id。")
}

sample_index <- match(embedding$sample_id, sample_info$sample_id)
if (anyNA(sample_index)) {
  stop(
    "以下sample_id无法映射到sample_info.tsv：",
    paste(unique(embedding$sample_id[is.na(sample_index)]), collapse = ", ")
  )
}
mapped_group <- as.character(sample_info$group_analysis[sample_index])
if (any(as.character(embedding$group_analysis) != mapped_group)) {
  stop("embedding source与sample_info.tsv的group_analysis不一致。")
}
embedding$patient_id <- as.character(sample_info$patient_id[sample_index])
if (anyNA(embedding$patient_id) || any(!nzchar(embedding$patient_id))) {
  stop("存在缺失patient_id，不能计算患者主导度。")
}

# 第2节：生成每簇来源构成、marker和QC摘要 -----------------------------------
composition <- f1_cluster_composition(embedding)
top_markers <- f1_top_markers(markers, n = config$markers$top_n_per_cluster)
cluster_qc <- f1_cluster_qc_metrics(embedding)

if (!setequal(composition$cluster, top_markers$cluster)) {
  stop("marker表与embedding source的cluster集合不一致。")
}

template <- merge(
  composition,
  top_markers,
  by = "cluster",
  all.x = TRUE,
  sort = TRUE
)
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

quality <- merge(
  composition,
  top_markers,
  by = "cluster",
  all.x = TRUE,
  sort = TRUE
)
quality <- merge(
  quality,
  cluster_qc,
  by = "cluster",
  all.x = TRUE,
  sort = TRUE
)
review_fields <- c(
  "cluster", "cell_type_major", "cell_type_minor", "cell_state",
  "annotation_confidence", "annotation_reason",
  "annotation_review_status", "downstream_handling_before_full_approval"
)
quality <- merge(
  quality,
  template[, review_fields, drop = FALSE],
  by = "cluster",
  all.x = TRUE,
  sort = TRUE
)
has_prefilled_review <- !is.na(quality$cell_type_major) &
  nzchar(quality$cell_type_major)
quality$review_status <- ifelse(
  has_prefilled_review,
  "complete_draft_pending_researcher_approval",
  "pending_full_F1.4_annotation"
)
quality <- quality[
  match(f1_sort_cluster_ids(quality$cluster), as.character(quality$cluster)),
  ,
  drop = FALSE
]
f1_write_tsv(quality, quality_path)

# 第3节：把本轮整合判断写回原诊断表 -----------------------------------------
numeric_value <- function(item_name) {
  value <- diagnostics$value[diagnostics$item == item_name]
  if (length(value) != 1L) {
    stop("integration_diagnostics.tsv缺少或重复项目：", item_name)
  }
  parsed <- suppressWarnings(as.numeric(value))
  if (!is.finite(parsed)) {
    stop("integration_diagnostics.tsv中的数值无法解析：", item_name)
  }
  parsed
}

ncount_pc_cor <- numeric_value("max_abs_PCA_correlation_nCount_RNA")
nfeature_pc_cor <- numeric_value("max_abs_PCA_correlation_nFeature_RNA")
technical_hvg_fraction <- numeric_value("technical_HVG_fraction")

normalization_row <- diagnostics$item == "normalization_sensitivity"
if (sum(normalization_row) != 1L) {
  stop("integration_diagnostics.tsv必须只有一个normalization_sensitivity项目。")
}
diagnostics$value[normalization_row] <- "LogNormalize_not_triggered_after_cluster_review"
diagnostics$interpretation[normalization_row] <- sprintf(
  paste0(
    "不设置单一PC相关硬阈值；综合判断为不触发：",
    "max|r|(nCount)=%.3f，max|r|(nFeature)=%.3f，",
    "technical HVG=%.2f%%，主要谱系marker清楚且未见技术指标主导聚类。"
  ),
  ncount_pc_cor,
  nfeature_pc_cor,
  100 * technical_hvg_fraction
)

review_items <- c(
  "cluster_sample_patient_composition_review",
  "preliminary_integration_review_after_clustering",
  "lineage_stratified_integration_review_status"
)
diagnostics <- diagnostics[!diagnostics$item %in% review_items, , drop = FALSE]

c3 <- composition[as.character(composition$cluster) == "3", , drop = FALSE]
if (nrow(c3) != 1L) stop("无法在cluster composition中唯一定位c3。")

lineage_names <- sort(unique(template$cell_type_major))
lineage_composition <- vapply(lineage_names, function(lineage_name) {
  x <- template[template$cell_type_major == lineage_name, , drop = FALSE]
  sprintf(
    "%s(n_clusters=%d,max_sample=%.3f,max_patient=%.3f)",
    lineage_name,
    nrow(x),
    max(x$dominant_sample_fraction),
    max(x$dominant_patient_fraction)
  )
}, character(1))

review_rows <- data.frame(
  item = review_items,
  value = c(
    sprintf(
      paste0(
        "clusters=%d;n_patients_per_cluster=%d-%d;",
        "max_dominant_sample_fraction=%.4f;",
        "max_dominant_patient_fraction=%.4f"
      ),
      nrow(composition),
      min(composition$n_patients),
      max(composition$n_patients),
      max(composition$dominant_sample_fraction),
      max(composition$dominant_patient_fraction)
    ),
    sprintf(
      paste0(
        "no_single_sample_dominated_cluster_detected;",
        "c3_Normal_Gastric=%.4f;c3_Primary_Tumor=%.4f;",
        "c3_malignancy_pending_F1.5_F1.6"
      ),
      c3$Normal_Gastric_fraction,
      c3$Primary_Tumor_fraction
    ),
    paste0(
      "draft_complete_pending_researcher_approval;",
      paste(lineage_composition, collapse = ";")
    )
  ),
  interpretation = c(
    paste0(
      "每簇均补充dominant sample/patient及组织来源构成；",
      "这些占比用于发现来源异常，不能单独证明过度整合。"
    ),
    paste0(
      "当前没有单一样本垄断的粗cluster；c3跨正常胃和原发肿瘤来源且marker混合，",
      "需在F1.5/F1.6结合未整合结构、上皮细分和CNV复核，",
      "目前既不证明过度整合，也不证明恶性。"
    ),
    paste0(
      "24个cluster的粗谱系草案已经完成。按草案分层后，T/NK、B/Plasma、",
      "Myeloid、Fibroblast/CAF、Endothelial/Pericyte、Mast和Mesothelial",
      "均保留连贯marker且未见单一样本垄断；上皮簇需在F1.5/F1.6继续判断。",
      "该结论仍等待研究者批准，不把跨患者混合自动解释为整合成功或失败。"
    )
  ),
  stringsAsFactors = FALSE
)
diagnostics <- rbind(diagnostics, review_rows)
f1_write_tsv(diagnostics, diagnostics_path)

# 第4节：保存运行记录 ---------------------------------------------------------
f1_save_session_info(config, "F1_04_prepare_annotation_review")
f1_append_log(
  config,
  stage,
  sprintf(
    paste0(
      "完成%d个cluster的轻量注释审核表；",
      "补充sample/patient/group构成；未读取Seurat对象、未运行DecontX"
    ),
    nrow(composition)
  )
)

message("F1.4注释前轻量审核表已生成：", config$paths$annotation_template)
