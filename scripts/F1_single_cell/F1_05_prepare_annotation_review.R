# F1.5 上皮cluster注释审核图 ---------------------------------------------
#
# 该图只帮助研究者回答三个问题：
# 1. 每个cluster更接近哪一种胃上皮分化程序；
# 2. 是否主要表现为增殖、应激、炎症等状态；
# 3. 是否存在免疫、基质或内皮污染。
# 图中任何单个marker都不能直接决定良恶性；恶性身份留到F1.6。

file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), winslash = "/"))
source(file.path(script_dir, "F1_config.R"), encoding = "UTF-8")
source(file.path(script_dir, "F1_utils.R"), encoding = "UTF-8")

args <- f1_parse_args()
config <- f1_build_config(args$project_root)
stage <- "F1.5 annotation review"
output_path <- file.path(
  config$paths$annotation_dir,
  "F1_epithelial_marker_review_dotplot.pdf"
)
feature_path <- file.path(
  config$paths$annotation_dir,
  "F1_epithelial_marker_review_features.tsv"
)
doublet_path <- file.path(
  config$paths$annotation_dir,
  "epithelial_cluster_doublet_sensitivity.tsv"
)

if (!args$execute) {
  f1_stage_dry_run(
    stage,
    "绘制上皮亚型、状态和污染marker审核图",
    c(config$paths$object_04),
    c(output_path, feature_path, doublet_path),
    c("Seurat", "SeuratObject", "ggplot2", "data.table")
  )
  quit(save = "no", status = 0)
}

f1_require_packages(
  c("Seurat", "SeuratObject", "ggplot2", "data.table"),
  stage
)
epithelial <- readRDS(config$paths$object_04)
if (!"epithelial_cluster_id" %in% colnames(epithelial[[]])) {
  stop("F1.5对象缺少epithelial_cluster_id。")
}
SeuratObject::DefaultAssay(epithelial) <- "RNA"
SeuratObject::Idents(epithelial) <- "epithelial_cluster_id"
meta <- epithelial[[]]
f1_require_columns(
  meta,
  c(
    "epithelial_cluster_id", "scDblFinder_score", "scDblFinder_class",
    "DoubletFinder_class"
  ),
  "F1.5 epithelial metadata"
)

feature_groups <- list(
  Broad_epithelial = c("EPCAM", "KRT8", "KRT18", "KRT19", "CDH1"),
  Pit_surface_mucous = c("MUC5AC", "TFF1", "GKN1", "GKN2"),
  Mucous_neck = c("MUC6", "TFF2"),
  Chief = c("LIPF", "PGA3", "PGA4", "PGA5", "PGC"),
  Parietal = c("ATP4A", "ATP4B", "GIF", "CBLIF"),
  Enteroendocrine = c("CHGA", "CHGB", "TPH1", "GHRL", "NEUROD1"),
  Intestinal_or_tumor_like = c(
    "REG4", "TFF3", "CEACAM5", "CEACAM6", "MMP7",
    "S100P", "OLFM4", "SOX9", "KRT17", "CLDN4"
  ),
  Cell_state = c(
    "MKI67", "TOP2A", "MT1G", "MT2A", "HSPA1A",
    "FOS", "CXCL8", "HLA-DRA"
  ),
  Non_epithelial_conflict = c(
    "PTPRC", "TYROBP", "LST1", "CD3D", "MS4A1",
    "CD79A", "MZB1", "JCHAIN", "COL1A1", "PECAM1"
  )
)
available <- rownames(epithelial[["RNA"]])
feature_groups <- lapply(feature_groups, intersect, y = available)
feature_groups <- feature_groups[lengths(feature_groups) > 0L]

feature_table <- do.call(rbind, lapply(names(feature_groups), function(group) {
  data.frame(
    marker_group = group,
    gene = feature_groups[[group]],
    stringsAsFactors = FALSE
  )
}))
f1_write_tsv(feature_table, feature_path)

# scDblFinder阳性细胞已经在F1.2删除；这里仅查看保留细胞中的残余分数，
# 以及未用于删除的DoubletFinder敏感性标签是否集中在某些新cluster。
cluster_ids <- f1_sort_cluster_ids(meta$epithelial_cluster_id)
doublet_summary <- do.call(rbind, lapply(cluster_ids, function(cluster_id) {
  keep <- as.character(meta$epithelial_cluster_id) == cluster_id
  df_class <- as.character(meta$DoubletFinder_class[keep])
  df_evaluable <- !is.na(df_class) & nzchar(df_class)
  data.frame(
    epithelial_cluster = cluster_id,
    n_cells = sum(keep),
    scDblFinder_score_median =
      stats::median(as.numeric(meta$scDblFinder_score[keep]), na.rm = TRUE),
    scDblFinder_score_p90 = as.numeric(stats::quantile(
      as.numeric(meta$scDblFinder_score[keep]),
      0.90,
      na.rm = TRUE,
      names = FALSE
    )),
    retained_scDblFinder_doublet_fraction = mean(
      as.character(meta$scDblFinder_class[keep]) == "doublet",
      na.rm = TRUE
    ),
    DoubletFinder_evaluable_cells = sum(df_evaluable),
    DoubletFinder_doublet_fraction = if (any(df_evaluable)) {
      mean(tolower(df_class[df_evaluable]) == "doublet")
    } else {
      NA_real_
    },
    stringsAsFactors = FALSE
  )
}))
f1_write_tsv(doublet_summary, doublet_path)

plot <- Seurat::DotPlot(
  epithelial,
  features = feature_groups,
  assay = "RNA",
  group.by = "epithelial_cluster_id",
  dot.scale = 6
) +
  Seurat::RotatedAxis() +
  ggplot2::labs(
    title = "F1.5 epithelial cluster marker review",
    x = NULL,
    y = "Epithelial cluster"
  ) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(size = 7),
    axis.text.y = ggplot2::element_text(size = 8),
    strip.text.x = ggplot2::element_text(size = 7)
  )

ggplot2::ggsave(
  output_path,
  plot,
  width = 21,
  height = 9,
  limitsize = FALSE
)
f1_save_session_info(config, "F1_05_prepare_annotation_review")
f1_append_log(
  config,
  stage,
  sprintf(
    "完成：%d个marker、%d个上皮cluster",
    nrow(feature_table),
    length(unique(epithelial$epithelial_cluster_id))
  )
)
message("F1.5上皮marker审核图完成：", output_path)
