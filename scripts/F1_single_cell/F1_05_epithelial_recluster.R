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
epithelial_cells <- colnames(all_cells)[eligible]

input_selection <- do.call(rbind, lapply(sort(unique(as.character(all_cells$seurat_clusters))), function(cluster_id) {
  cells <- colnames(all_cells)[as.character(all_cells$seurat_clusters) == cluster_id]
  selected <- intersect(cells, epithelial_cells)
  data.frame(
    source_all_cell_cluster = cluster_id,
    source_cell_type_major = unique(as.character(all_cells$cell_type_major[match(cells, colnames(all_cells))]))[[1]],
    source_annotation_confidence = paste(sort(unique(as.character(all_cells$annotation_confidence[match(cells, colnames(all_cells))]))), collapse = "|"),
    total_cells = length(cells),
    included_epithelial_cells = length(selected),
    inclusion_rule = "cell_type_major=Epithelial_and_confidence_high_or_medium",
    stringsAsFactors = FALSE
  )
}))
f1_write_tsv(input_selection, file.path(config$paths$annotation_dir, "epithelial_input_selection.tsv"))

all_cells <- f1_join_assay(all_cells, "RNA")
raw <- SeuratObject::LayerData(all_cells, assay = "RNA", layer = "counts")[, epithelial_cells, drop = FALSE]
f1_assert_integer_counts(raw, "F1.5 epithelial RNA counts")
meta <- all_cells[[]][epithelial_cells, , drop = FALSE]
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
    "input_epithelial_cells", "input_samples", "normalization", "vars_to_regress",
    "Harmony_group_by", "main_dims", "default_resolution", "epithelial_clusters",
    "malignancy_status"
  ),
  value = c(
    ncol(epithelial), length(unique(epithelial$sample_id)), "per_sample_SCTransform_v2", "NULL",
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
