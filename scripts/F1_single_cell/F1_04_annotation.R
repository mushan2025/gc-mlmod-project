# F1.4 主要细胞类型注释 ------------------------------------------------------
#
# 生物学目的：依据cluster marker、冻结marker panel和样本覆盖，为每个cluster
# 赋予可解释的主要谱系。最终标签必须经研究者审核，脚本不从单个marker自动下结论。
#
# 第一次正式运行会生成注释模板并停止；填好批准表后再次运行才保存正式注释对象。

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
  file.path(config$paths$annotation_dir, "celltype_annotation_summary.tsv")
)

if (!args$execute) {
  f1_stage_dry_run(
    stage,
    "计算RNA marker并通过批准表完成主要谱系人工注释",
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
f1_append_log(config, stage, "开始cluster marker计算与主要谱系注释准备")

object <- readRDS(config$paths$object_03a)
f1_require_columns(
  object[[]],
  c("sample_id", "group_analysis", "seurat_clusters", "nCount_RNA", "nFeature_RNA", "mt_percent", "HB_percent"),
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
cluster_counts <- as.data.frame(table(cluster = as.character(object$seurat_clusters)), stringsAsFactors = FALSE)
colnames(cluster_counts)[2] <- "cell_count"
sample_counts <- aggregate(
  sample_id ~ seurat_clusters,
  data = unique(object[[]][, c("sample_id", "seurat_clusters")]),
  FUN = length
)
colnames(sample_counts) <- c("cluster", "sample_count")
template <- merge(cluster_counts, sample_counts, by = "cluster", all.x = TRUE, sort = TRUE)
template <- merge(template, top, by = "cluster", all.x = TRUE, sort = TRUE)
template$cell_type_major <- ""
template$cell_type_minor <- ""
template$cell_state <- ""
template$annotation_confidence <- ""
template$annotation_reason <- ""
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
    ggplot2::labs(title = "主要谱系marker（供人工注释审核）", x = NULL, y = "cluster") +
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
allowed_major <- c(
  "Epithelial", "T/NK", "B/Plasma", "Myeloid", "Fibroblast/CAF",
  "Endothelial/Pericyte", "Mast", "Uncertain", "Mixed_or_doublet_suspect"
)
if (any(!approved$cell_type_major %in% allowed_major)) {
  stop("批准表含未登记的cell_type_major：", paste(unique(approved$cell_type_major[!approved$cell_type_major %in% allowed_major]), collapse = ", "))
}
if (any(!approved$annotation_confidence %in% c("high", "medium", "low"))) {
  stop("annotation_confidence只允许high、medium或low。")
}

map_index <- match(as.character(object$seurat_clusters), as.character(approved$cluster))
for (field in setdiff(required_fields, "cluster")) {
  object[[field]] <- as.character(approved[[field]][map_index])
}

annotation_summary <- merge(template[, c("cluster", "cell_count", "sample_count", "top_markers")], approved, by = "cluster", sort = TRUE)
f1_write_tsv(annotation_summary, file.path(config$paths$annotation_dir, "celltype_annotation_summary.tsv"))

meta <- object[[]]
cluster_quality <- do.call(rbind, lapply(clusters, function(cluster_id) {
  x <- meta[as.character(meta$seurat_clusters) == cluster_id, , drop = FALSE]
  marker_row <- top[top$cluster == cluster_id, , drop = FALSE]
  data.frame(
    cluster = cluster_id,
    cell_count = nrow(x),
    sample_count = length(unique(x$sample_id)),
    top_markers = marker_row$top_markers %||% "",
    median_nCount_RNA = stats::median(x$nCount_RNA, na.rm = TRUE),
    median_nFeature_RNA = stats::median(x$nFeature_RNA, na.rm = TRUE),
    median_mt_percent = stats::median(x$mt_percent, na.rm = TRUE),
    median_HB_percent = stats::median(x$HB_percent, na.rm = TRUE),
    DoubletFinder_only_fraction = mean(
      tolower(x$DoubletFinder_class) == "doublet" & tolower(x$scDblFinder_class) == "singlet",
      na.rm = TRUE
    ),
    cell_type_major = unique(x$cell_type_major)[[1]],
    annotation_confidence = unique(x$annotation_confidence)[[1]],
    review_status = "researcher_approved_annotation",
    stringsAsFactors = FALSE
  )
}))
f1_write_tsv(cluster_quality, file.path(config$paths$annotation_dir, "cluster_quality_summary.tsv"))

p1 <- Seurat::DimPlot(object, reduction = "umap", group.by = "cell_type_major", label = TRUE, repel = TRUE) +
  ggplot2::ggtitle("F1主要细胞类型")
p2 <- Seurat::DimPlot(object, reduction = "umap", group.by = "sample_id") +
  ggplot2::ggtitle("F1主要细胞类型图谱：样本来源")
ggplot2::ggsave(
  file.path(config$paths$annotation_dir, "F1_all_cells_annotated_umap.pdf"),
  patchwork::wrap_plots(p1, p2, ncol = 2),
  width = 14,
  height = 6,
  limitsize = FALSE
)

f1_save_rds_atomic(object, config$paths$object_03, compress = FALSE)
f1_save_session_info(config, "F1_04_annotation")
f1_append_log(config, stage, sprintf("完成研究者批准的主要谱系注释：%d个cluster", length(clusters)))
message("F1.4完成：", config$paths$object_03)
