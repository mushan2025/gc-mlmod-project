# F1.5 上皮细胞提取与二次聚类 ----------------------------------------------
#
# 生物学目的：把可靠上皮细胞从全细胞图谱中提取出来，重新建立只反映上皮内部
# 异质性的SCTransform/Harmony图谱。上皮亚型与“是否恶性”在这里仍保持分开。
#
# 主要输入：03_all_cells_integrated_annotated.rds。
# 主要输出：04_epithelial_reclustered.rds、上皮marker和污染复核模板。

file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), winslash = "/"))
source(file.path(script_dir, "F1_config.R"), encoding = "UTF-8")
source(file.path(script_dir, "F1_utils.R"), encoding = "UTF-8")

args <- f1_parse_args()
config <- f1_build_config(args$project_root)
stage <- "F1.5"
marker_path <- file.path(config$paths$annotation_dir, "epithelial_cluster_marker_genes.tsv")
outputs <- c(
  config$paths$object_04,
  file.path(config$paths$annotation_dir, "epithelial_input_selection.tsv"),
  file.path(config$paths$annotation_dir, "epithelial_sample_input_selection.tsv"),
  file.path(config$paths$annotation_dir, "epithelial_excluded_cells.tsv"),
  file.path(config$paths$annotation_dir, "epithelial_recluster_summary.tsv"),
  config$paths$epithelial_review_template
)

if (!args$execute) {
  f1_stage_dry_run(
    stage,
    "提取高/中置信上皮细胞并用相同SCTransform v2设置二次聚类",
    c(config$paths$object_03),
    outputs,
    config$packages[[stage]]
  )
  quit(save = "no", status = 0)
}

f1_require_packages(config$packages[[stage]], stage)
f1_check_f0_ready(config)
if (!file.exists(config$paths$object_03)) stop("请先完成F1.4批准注释：", config$paths$object_03)
f1_prepare_directories(config)
set.seed(config$seed)
f1_append_log(config, stage, "开始提取高/中置信上皮细胞并二次聚类")

all_cells <- readRDS(config$paths$object_03)
required_meta <- c(
  "cell_type_major", "annotation_confidence", "sample_id", "patient_id",
  "group_analysis", "seurat_clusters"
)
f1_require_columns(all_cells[[]], required_meta, "F1.4 annotated object metadata")
eligible <- all_cells$cell_type_major == "Epithelial" &
  all_cells$annotation_confidence %in% c("high", "medium")
if (!any(eligible)) stop("F1.4批准注释中没有高/中置信Epithelial细胞。")
initial_epithelial_cells <- colnames(all_cells)[eligible]

all_cells <- f1_join_assay(all_cells, "RNA")
raw_initial <- SeuratObject::LayerData(
  all_cells,
  assay = "RNA",
  layer = "counts"
)[, initial_epithelial_cells, drop = FALSE]

# 逐样本SCTransform需要每个样本自身提供足够的表达特征。正式首跑发现
# sample37仅有7个候选上皮细胞；按SCTransform固定的min_cells=5计算时只有
# 283个可用基因，会把所有样本共同PCA特征压到500以下。研究者于2026-07-30
# 批准：不足500个SCT可用基因的样本不进入F1.5上皮重聚类，但细胞仍保留在
# F1.4全细胞对象并写入明确排除记录。
minimum_detected_features_for_sample_sct <- 500L
minimum_cells_per_gene_for_sct <- config$sct$minimum_cells_per_gene
sample_ids <- sort(unique(as.character(
  all_cells$sample_id[match(initial_epithelial_cells, colnames(all_cells))]
)))
sample_selection <- do.call(rbind, lapply(sample_ids, function(sample_id) {
  sample_cells <- initial_epithelial_cells[
    as.character(all_cells$sample_id[match(initial_epithelial_cells, colnames(all_cells))]) ==
      sample_id
  ]
  sct_usable_features <- sum(
    Matrix::rowSums(raw_initial[, sample_cells, drop = FALSE] > 0) >=
      minimum_cells_per_gene_for_sct
  )
  include_sample <-
    sct_usable_features >= minimum_detected_features_for_sample_sct
  sample_groups <- unique(as.character(
    all_cells$group_analysis[match(sample_cells, colnames(all_cells))]
  ))
  data.frame(
    sample_id = sample_id,
    group_analysis = paste(sort(sample_groups), collapse = "|"),
    initial_epithelial_cells = length(sample_cells),
    sct_min_cells_per_gene = minimum_cells_per_gene_for_sct,
    sct_usable_features = sct_usable_features,
    minimum_sct_usable_features_required =
      minimum_detected_features_for_sample_sct,
    include_in_epithelial_recluster = include_sample,
    exclusion_reason = if (include_sample) {
      ""
    } else {
      "fewer_than_500_SCT_usable_features_cannot_support_independent_sample_SCT_model"
    },
    stringsAsFactors = FALSE
  )
}))
f1_write_tsv(
  sample_selection,
  file.path(config$paths$annotation_dir, "epithelial_sample_input_selection.tsv")
)

included_samples <- sample_selection$sample_id[
  sample_selection$include_in_epithelial_recluster
]
epithelial_cells <- initial_epithelial_cells[
  as.character(all_cells$sample_id[match(initial_epithelial_cells, colnames(all_cells))]) %in%
    included_samples
]
excluded_cells <- setdiff(initial_epithelial_cells, epithelial_cells)
excluded_meta <- all_cells[[]][excluded_cells, , drop = FALSE]
excluded_table <- data.frame(
  cell_id_final = excluded_cells,
  sample_id = as.character(excluded_meta$sample_id),
  group_analysis = as.character(excluded_meta$group_analysis),
  source_all_cell_cluster = as.character(excluded_meta$seurat_clusters),
  exclusion_reason = rep(
    "sample_has_fewer_than_500_SCT_usable_features_for_independent_SCT_model",
    length(excluded_cells)
  ),
  retained_in_F1_4_all_cell_object = rep(TRUE, length(excluded_cells)),
  stringsAsFactors = FALSE
)
f1_write_tsv(
  excluded_table,
  file.path(config$paths$annotation_dir, "epithelial_excluded_cells.tsv")
)
if (!length(epithelial_cells)) {
  stop("应用逐样本SCT最低表达特征规则后没有可用于F1.5的上皮细胞。")
}

input_selection <- do.call(rbind, lapply(
  sort(unique(as.character(all_cells$seurat_clusters))),
  function(cluster_id) {
    cells <- colnames(all_cells)[as.character(all_cells$seurat_clusters) == cluster_id]
    initial_selected <- intersect(cells, initial_epithelial_cells)
    selected <- intersect(cells, epithelial_cells)
    data.frame(
      source_all_cell_cluster = cluster_id,
      source_cell_type_major =
        unique(as.character(all_cells$cell_type_major[match(cells, colnames(all_cells))]))[[1]],
      source_annotation_confidence = paste(
        sort(unique(as.character(
          all_cells$annotation_confidence[match(cells, colnames(all_cells))]
        ))),
        collapse = "|"
      ),
      total_cells = length(cells),
      initial_eligible_epithelial_cells = length(initial_selected),
      included_epithelial_cells = length(selected),
      technical_excluded_cells = length(initial_selected) - length(selected),
      inclusion_rule = paste0(
        "cell_type_major=Epithelial_and_confidence_high_or_medium;",
        "sample_SCT_usable_features>=500_with_min_cells_per_gene=5"
      ),
      stringsAsFactors = FALSE
    )
  }
))
f1_write_tsv(
  input_selection,
  file.path(config$paths$annotation_dir, "epithelial_input_selection.tsv")
)

raw <- raw_initial[, epithelial_cells, drop = FALSE]
rm(raw_initial)
f1_assert_integer_counts(raw, "F1.5 epithelial RNA counts")
meta <- all_cells[[]][epithelial_cells, , drop = FALSE]
# 保留全细胞图谱中的原始cluster，避免二次聚类覆盖seurat_clusters后失去来源。
meta$source_all_cell_cluster <- as.character(meta$seurat_clusters)
epithelial <- Seurat::CreateSeuratObject(
  counts = raw,
  assay = "RNA",
  project = "GSE183904_epithelial",
  meta.data = meta,
  min.cells = 0,
  min.features = 0
)
rm(raw)
epithelial <- f1_run_sct_harmony(epithelial, config, stage)
epithelial$epithelial_cluster_id <- as.character(epithelial$seurat_clusters)
SeuratObject::Idents(epithelial) <- "epithelial_cluster_id"

markers <- Seurat::FindAllMarkers(
  epithelial,
  assay = "RNA",
  only.pos = config$markers$only_positive,
  min.pct = config$markers$min_pct,
  logfc.threshold = config$markers$logfc_threshold,
  test.use = "wilcox",
  verbose = TRUE
)
if (!"gene" %in% colnames(markers)) markers$gene <- rownames(markers)
f1_write_tsv(markers, marker_path)
top <- f1_top_markers(markers, n = config$markers$top_n_per_cluster)

clusters <- sort(unique(epithelial$epithelial_cluster_id))
review_template <- do.call(rbind, lapply(clusters, function(cluster_id) {
  cells <- colnames(epithelial)[epithelial$epithelial_cluster_id == cluster_id]
  groups <- table(epithelial$group_analysis[match(cells, colnames(epithelial))])
  top_row <- top[top$cluster == cluster_id, , drop = FALSE]
  data.frame(
    epithelial_cluster = cluster_id,
    cell_count = length(cells),
    sample_count = length(unique(epithelial$sample_id[match(cells, colnames(epithelial))])),
    group_cell_counts = paste0(names(groups), "=", as.integer(groups), collapse = "|"),
    top_markers = top_row$top_markers %||% "",
    include_in_malignancy = "",
    contamination_status = "",
    epithelial_subtype = "",
    review_reason = "",
    stringsAsFactors = FALSE
  )
}))
f1_write_tsv(review_template, config$paths$epithelial_review_template)

summary <- data.frame(
  item = c(
    "initial_epithelial_cells", "technical_excluded_cells",
    "technical_excluded_samples", "SCT_min_cells_per_gene",
    "minimum_SCT_usable_features_per_sample",
    "input_epithelial_cells", "input_samples", "normalization", "vars_to_regress",
    "Harmony_group_by", "main_dims", "default_resolution", "epithelial_clusters",
    "malignancy_status"
  ),
  value = c(
    length(initial_epithelial_cells), length(excluded_cells),
    paste(setdiff(sample_ids, included_samples), collapse = "|"),
    minimum_cells_per_gene_for_sct,
    minimum_detected_features_for_sample_sct,
    ncol(epithelial), length(unique(epithelial$sample_id)),
    "per_sample_SCTransform_v2", "NULL",
    config$sct$harmony_group, paste(range(config$sct$main_dims), collapse = ":"),
    config$sct$default_resolution, length(clusters), "not_assigned_in_F1.5"
  ),
  stringsAsFactors = FALSE
)
f1_write_tsv(summary, file.path(config$paths$annotation_dir, "epithelial_recluster_summary.tsv"))

p1 <- Seurat::DimPlot(epithelial, reduction = "umap", group.by = "epithelial_cluster_id", label = TRUE, repel = TRUE) +
  ggplot2::ggtitle("Epithelial-cell reclustering")
p2 <- Seurat::DimPlot(epithelial, reduction = "umap", group.by = "group_analysis") +
  ggplot2::ggtitle("Epithelial cells: tissue group")
p3 <- Seurat::DimPlot(epithelial, reduction = "umap", group.by = "sample_id") +
  ggplot2::ggtitle("Epithelial cells: sample")
ggplot2::ggsave(
  file.path(config$paths$annotation_dir, "F1_epithelial_recluster_umap.pdf"),
  patchwork::wrap_plots(p1, p2, p3, ncol = 2),
  width = 14,
  height = 10,
  limitsize = FALSE
)

f1_save_rds_atomic(epithelial, config$paths$object_04, compress = FALSE)
f1_save_session_info(config, "F1_05_epithelial_recluster")
f1_append_log(
  config,
  stage,
  sprintf("完成：%d个高/中置信上皮细胞，%d个二次聚类；尚未判定恶性", ncol(epithelial), length(clusters))
)
message("F1.5完成：", config$paths$object_04)
message("F1.6前请审核：", config$paths$epithelial_review_template)
