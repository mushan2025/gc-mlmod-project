# F1.3 SCTransform、Harmony与全细胞聚类 ------------------------------------
#
# 生物学目的：在不回归线粒体比例的前提下稳定不同测序深度的表达方差，建立可用于
# 主要细胞谱系注释的全细胞图谱。Harmony只校正sample_id的低维坐标，不改表达矩阵。
#
# 主要输入：02_all_cells_qc_filtered.rds。
# 主要输出：03a_all_cells_sct_harmony_clustered.rds、整合诊断和UMAP图。

file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), winslash = "/"))
source(file.path(script_dir, "F1_config.R"), encoding = "UTF-8")
source(file.path(script_dir, "F1_utils.R"), encoding = "UTF-8")

args <- f1_parse_args()
config <- f1_build_config(args$project_root)
stage <- "F1.3"
outputs <- c(
  config$paths$object_03a,
  file.path(config$paths$annotation_dir, "integration_diagnostics.tsv"),
  file.path(config$paths$annotation_dir, "resolution_cluster_counts.tsv"),
  file.path(config$paths$annotation_dir, "F1_all_cells_embedding_source_data.tsv.gz")
)

if (!args$execute) {
  f1_stage_dry_run(
    stage,
    "逐样本SCTransform v2，不回归变量；SCT PCA后按sample_id运行Harmony和Leiden",
    c(config$paths$object_02),
    outputs,
    config$packages[[stage]]
  )
  quit(save = "no", status = 0)
}

f1_require_packages(config$packages[[stage]], stage)
f1_check_f0_ready(config)
if (!file.exists(config$paths$object_02)) stop("请先完成F1.2：", config$paths$object_02)
f1_prepare_directories(config)
set.seed(config$seed)
f1_append_log(config, stage, "开始SCTransform v2、PCA、Harmony和Leiden全细胞聚类")

object <- readRDS(config$paths$object_02)
forbidden <- grep("MLMOD|UCell|prognos|survival", colnames(object[[]]), ignore.case = TRUE, value = TRUE)
if (length(forbidden)) {
  stop("F1对象中出现不应参与F1的结果字段：", paste(forbidden, collapse = ", "))
}
f1_require_columns(
  object[[]],
  c("sample_id", "group_analysis", "nCount_RNA", "nFeature_RNA", "mt_percent", "HB_percent"),
  "F1.2 object metadata"
)

object <- f1_run_sct_harmony(object, config, stage)

hvg <- SeuratObject::VariableFeatures(object[["SCT"]])
technical_hvg <- grepl("^MT-|^RPL|^RPS", toupper(hvg)) |
  toupper(hvg) %in% config$qc$globin_panel
pca <- Seurat::Embeddings(object, reduction = "pca")
pca_input_features_n <- nrow(Seurat::Loadings(object, reduction = "pca"))
qc_fields <- intersect(c("nCount_RNA", "nFeature_RNA", "mt_percent", "HB_percent"), colnames(object[[]]))
qc_cor <- vapply(qc_fields, function(field) {
  max(abs(stats::cor(pca, object[[field, drop = TRUE]], use = "pairwise.complete.obs")))
}, numeric(1))

diagnostics <- data.frame(
  item = c(
    "normalization", "SCTransform_vst_flavor", "SCTransform_vars_to_regress",
    "variable_features_n", "technical_HVG_fraction", "PCA_npcs", "PCA_input_features_n",
    "Harmony_group_by",
    "main_dims", "Leiden_algorithm", "default_resolution",
    paste0("max_abs_PCA_correlation_", names(qc_cor)),
    "unintegrated_reference", "main_embedding", "normalization_sensitivity"
  ),
  value = c(
    "per_sample_SCTransform_v2", config$sct$vst_flavor, "NULL",
    length(hvg), mean(technical_hvg), config$sct$pca_npcs, pca_input_features_n,
    config$sct$harmony_group,
    paste(range(config$sct$main_dims), collapse = ":"), config$sct$leiden_algorithm,
    config$sct$default_resolution, as.character(qc_cor),
    "SCT_PCA_and_umap.sct_pca_retained", "Harmony_sample_id_UMAP_Leiden",
    "LogNormalize_not_triggered_no_specific_concern"
  ),
  interpretation = c(
    "主归一化", "固定v2方差稳定化", "不回归mt/nCount/细胞周期",
    "SCT候选高变基因数", "只作技术基因富集提示，不自动删基因", "PCA维数",
    "各样本SCT共同保留且实际进入PCA的已缩放高变基因数", "仅校正样本低维坐标",
    "预登记主维数", "Seurat algorithm 4", "marker审核前默认值",
    rep("用于发现QC变量是否明显主导PC；不据此反向改QC阈值", length(qc_cor)),
    "用于识别Harmony是否抹除患者特异肿瘤结构", "主要细胞大类聚类空间",
    "只有明确深度驱动或审稿质疑时才运行LogNormalize敏感性"
  ),
  stringsAsFactors = FALSE
)
f1_write_tsv(diagnostics, file.path(config$paths$annotation_dir, "integration_diagnostics.tsv"))

resolution_columns <- vapply(config$sct$resolutions, function(resolution) {
  paste0("SCT_harmony_res.", format(resolution, trim = TRUE, scientific = FALSE))
}, character(1))
resolution_counts <- do.call(rbind, lapply(resolution_columns, function(column) {
  tab <- table(object[[column, drop = TRUE]])
  data.frame(
    resolution_field = column,
    cluster = names(tab),
    cell_count = as.integer(tab),
    stringsAsFactors = FALSE
  )
}))
f1_write_tsv(resolution_counts, file.path(config$paths$annotation_dir, "resolution_cluster_counts.tsv"))

umap_harmony <- Seurat::Embeddings(object, reduction = "umap")
umap_sct_pca <- Seurat::Embeddings(object, reduction = "umap.sct_pca")
embedding_source <- data.frame(
  cell_id = colnames(object),
  sample_id = object$sample_id,
  group_analysis = object$group_analysis,
  seurat_clusters = object$seurat_clusters,
  harmony_UMAP_1 = umap_harmony[, 1],
  harmony_UMAP_2 = umap_harmony[, 2],
  unintegrated_SCT_PCA_UMAP_1 = umap_sct_pca[, 1],
  unintegrated_SCT_PCA_UMAP_2 = umap_sct_pca[, 2],
  nCount_RNA = object$nCount_RNA,
  nFeature_RNA = object$nFeature_RNA,
  mt_percent = object$mt_percent,
  HB_percent = object$HB_percent,
  stringsAsFactors = FALSE
)
data.table::fwrite(
  embedding_source,
  file.path(config$paths$annotation_dir, "F1_all_cells_embedding_source_data.tsv.gz"),
  sep = "\t",
  quote = FALSE,
  compress = "gzip"
)

p_cluster <- Seurat::DimPlot(object, reduction = "umap", group.by = "seurat_clusters", label = TRUE, repel = TRUE) +
  ggplot2::ggtitle("Harmony UMAP: Leiden clusters (resolution 0.6)")
p_sample <- Seurat::DimPlot(object, reduction = "umap", group.by = "sample_id") +
  ggplot2::ggtitle("Harmony UMAP: sample")
p_group <- Seurat::DimPlot(object, reduction = "umap", group.by = "group_analysis") +
  ggplot2::ggtitle("Harmony UMAP: tissue group")
p_unintegrated <- Seurat::DimPlot(object, reduction = "umap.sct_pca", group.by = "sample_id") +
  ggplot2::ggtitle("Unintegrated SCT-PCA UMAP: sample")
ggplot2::ggsave(
  file.path(config$paths$annotation_dir, "F1_all_cells_embedding_review.pdf"),
  patchwork::wrap_plots(p_cluster, p_sample, p_group, p_unintegrated, ncol = 2),
  width = 14,
  height = 11,
  limitsize = FALSE
)

f1_save_rds_atomic(object, config$paths$object_03a, compress = FALSE)
f1_write_parameter_versions(config)
f1_save_session_info(config, "F1_03_sct_harmony_cluster")
f1_append_log(
  config,
  stage,
  sprintf("完成：%d细胞；SCTransform v2 vars.to.regress=NULL；Harmony仅sample_id；默认resolution=%s", ncol(object), config$sct$default_resolution)
)
message("F1.3完成：", config$paths$object_03a)
