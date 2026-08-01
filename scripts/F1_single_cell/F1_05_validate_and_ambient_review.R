# F1.5 输出验证与上皮ambient RNA复核 --------------------------------------
#
# 生物学目的：
# 1. 确认上皮二次聚类对象完整、raw counts未被替换，且排除记录与实际对象一致；
# 2. 按新的上皮cluster检查ambient污染分布；
# 3. 比较raw与DecontX corrected中的上皮亚型panel和前列marker，判断污染校正
#    是否足以改变F1.5注释。corrected矩阵仍只作敏感性分析。

file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), winslash = "/"))
source(file.path(script_dir, "F1_config.R"), encoding = "UTF-8")
source(file.path(script_dir, "F1_utils.R"), encoding = "UTF-8")

args <- f1_parse_args()
config <- f1_build_config(args$project_root)
stage <- "F1.5 validation"

validation_path <- file.path(
  config$paths$annotation_dir,
  "F1_5_output_validation.tsv"
)
composition_path <- file.path(
  config$paths$annotation_dir,
  "epithelial_cluster_composition.tsv"
)
resolution_path <- file.path(
  config$paths$annotation_dir,
  "epithelial_resolution_cluster_counts.tsv"
)
source_composition_path <- file.path(
  config$paths$annotation_dir,
  "epithelial_cluster_source_composition.tsv"
)
ambient_cluster_path <- file.path(
  config$paths$annotation_dir,
  "epithelial_ambient_summary_by_cluster.tsv"
)
ambient_sample_cluster_path <- file.path(
  config$paths$annotation_dir,
  "epithelial_ambient_summary_by_sample_cluster.tsv"
)
panel_path <- file.path(
  config$paths$annotation_dir,
  "epithelial_panel_raw_corrected_by_cluster.tsv"
)
impact_path <- file.path(
  config$paths$annotation_dir,
  "epithelial_ambient_annotation_impact_summary.tsv"
)

if (!args$execute) {
  f1_stage_dry_run(
    stage,
    "验证F1.5对象并复核上皮cluster的raw/corrected marker稳定性",
    c(
      config$paths$object_04,
      config$paths$ambient_cell_estimates,
      config$paths$decontx_corrected_dir,
      file.path(config$paths$annotation_dir, "epithelial_cluster_marker_genes.tsv")
    ),
    c(
      validation_path, composition_path, resolution_path,
      source_composition_path,
      ambient_cluster_path,
      ambient_sample_cluster_path, panel_path, impact_path
    ),
    c("Seurat", "SeuratObject", "Matrix", "data.table")
  )
  quit(save = "no", status = 0)
}

f1_require_packages(
  c("Seurat", "SeuratObject", "Matrix", "data.table"),
  stage
)
if (!file.exists(config$paths$object_04)) {
  stop("缺少F1.5正式对象：", config$paths$object_04)
}

epithelial <- readRDS(config$paths$object_04)
meta <- epithelial[[]]
f1_require_columns(
  meta,
  c(
    "sample_id", "patient_id", "group_analysis", "seurat_clusters",
    "epithelial_cluster_id", "source_all_cell_cluster"
  ),
  "F1.5 epithelial metadata"
)

sample_selection <- f1_read_tsv(file.path(
  config$paths$annotation_dir,
  "epithelial_sample_input_selection.tsv"
))
excluded_cells <- f1_read_tsv(file.path(
  config$paths$annotation_dir,
  "epithelial_excluded_cells.tsv"
))
f1_require_columns(
  sample_selection,
  c(
    "sample_id", "initial_epithelial_cells", "include_in_epithelial_recluster",
    "sct_usable_features"
  ),
  "epithelial_sample_input_selection.tsv"
)
sample_selection$include_in_epithelial_recluster <- f1_as_logical(
  sample_selection$include_in_epithelial_recluster,
  "include_in_epithelial_recluster"
)

raw <- SeuratObject::LayerData(epithelial, assay = "RNA", layer = "counts")
f1_assert_integer_counts(raw, "F1.5 epithelial RNA raw counts")

required_reductions <- c("pca", "umap.sct_pca", "harmony", "umap")
available_reductions <- names(epithelial@reductions)
reduction_ok <- all(required_reductions %in% available_reductions)
embedding_ok <- reduction_ok && all(vapply(required_reductions, function(name) {
  embedding <- Seurat::Embeddings(epithelial, reduction = name)
  all(is.finite(embedding)) && all(apply(embedding, 2, stats::sd) > 0)
}, logical(1)))

resolution_fields <- vapply(config$sct$resolutions, function(resolution) {
  paste0(
    "SCT_harmony_res.",
    format(resolution, trim = TRUE, scientific = FALSE)
  )
}, character(1))
expected_cells <- sum(
  sample_selection$initial_epithelial_cells[
    sample_selection$include_in_epithelial_recluster
  ]
)
expected_samples <- sum(sample_selection$include_in_epithelial_recluster)
excluded_not_in_object <- if (nrow(excluded_cells)) {
  !any(excluded_cells$cell_id_final %in% colnames(epithelial))
} else {
  TRUE
}
pca_features <- nrow(Seurat::Loadings(epithelial, reduction = "pca"))
sct_models <- length(methods::slot(epithelial[["SCT"]], "SCTModel.list"))
forbidden_fields <- grep(
  "MLMOD|prognos|survival",
  colnames(meta),
  value = TRUE,
  ignore.case = TRUE
)

validation <- data.frame(
  check = c(
    "cell_count_matches_approved_input",
    "sample_count_matches_approved_input",
    "excluded_cells_absent_from_F1_5",
    "source_all_cell_cluster_complete",
    "RNA_raw_counts_nonnegative_integer",
    "required_reductions_present",
    "required_embeddings_finite_nonblank",
    "required_resolution_fields_present",
    "PCA_features_at_least_500",
    "SCT_models_match_included_samples",
    "no_MLMOD_or_prognosis_fields"
  ),
  status = c(
    ncol(epithelial) == expected_cells,
    length(unique(epithelial$sample_id)) == expected_samples,
    excluded_not_in_object,
    all(!is.na(meta$source_all_cell_cluster)) &&
      all(nzchar(as.character(meta$source_all_cell_cluster))),
    TRUE,
    reduction_ok,
    embedding_ok,
    all(resolution_fields %in% colnames(meta)),
    pca_features >= 500L,
    sct_models == expected_samples,
    length(forbidden_fields) == 0L
  ),
  observed = c(
    paste0(ncol(epithelial), "/", expected_cells),
    paste0(length(unique(epithelial$sample_id)), "/", expected_samples),
    if (nrow(excluded_cells)) nrow(excluded_cells) else 0L,
    paste(sort(unique(as.character(meta$source_all_cell_cluster))), collapse = "|"),
    paste0(nrow(raw), " features; ", ncol(raw), " cells"),
    paste(intersect(required_reductions, available_reductions), collapse = "|"),
    embedding_ok,
    paste(intersect(resolution_fields, colnames(meta)), collapse = "|"),
    pca_features,
    paste0(sct_models, "/", expected_samples),
    paste(forbidden_fields, collapse = "|")
  ),
  stringsAsFactors = FALSE
)
validation$status <- ifelse(validation$status, "PASS", "FAIL")
f1_write_tsv(validation, validation_path)
if (any(validation$status == "FAIL")) {
  stop(
    "F1.5对象验证失败：",
    paste(validation$check[validation$status == "FAIL"], collapse = ", ")
  )
}

composition_meta <- meta
composition_meta$seurat_clusters <- as.character(
  composition_meta$epithelial_cluster_id
)
composition <- f1_cluster_composition(composition_meta)
f1_write_tsv(composition, composition_path)

resolution_summary <- data.frame(
  resolution = config$sct$resolutions,
  metadata_field = resolution_fields,
  n_clusters = vapply(resolution_fields, function(field) {
    length(unique(as.character(meta[[field]])))
  }, integer(1)),
  is_main = config$sct$resolutions == config$sct$default_resolution,
  stringsAsFactors = FALSE
)
f1_write_tsv(resolution_summary, resolution_path)

source_composition <- as.data.frame(table(
  epithelial_cluster = as.character(meta$epithelial_cluster_id),
  source_all_cell_cluster = as.character(meta$source_all_cell_cluster)
), stringsAsFactors = FALSE)
source_composition <- source_composition[
  source_composition$Freq > 0,
  ,
  drop = FALSE
]
source_totals <- tapply(
  source_composition$Freq,
  source_composition$epithelial_cluster,
  sum
)
source_composition$fraction_within_epithelial_cluster <-
  source_composition$Freq /
  source_totals[source_composition$epithelial_cluster]
colnames(source_composition)[colnames(source_composition) == "Freq"] <-
  "cell_count"
f1_write_tsv(source_composition, source_composition_path)

# 把F1.4已有的逐细胞污染估计映射到新的上皮cluster。
ambient <- f1_read_tsv(config$paths$ambient_cell_estimates)
f1_require_columns(
  ambient,
  c("cell_id_final", "retained_cell_ambient_contamination_estimate"),
  "ambient_rna_cell_estimates.tsv"
)
ambient_index <- match(colnames(epithelial), ambient$cell_id_final)
if (anyNA(ambient_index)) {
  stop("F1.5上皮细胞无法完整匹配F1.4 ambient污染估计。")
}
meta$ambient_estimate <- as.numeric(
  ambient$retained_cell_ambient_contamination_estimate[ambient_index]
)
if (any(!is.finite(meta$ambient_estimate))) {
  stop("F1.5上皮对象存在非有限ambient污染估计。")
}

summarize_ambient <- function(x) {
  data.frame(
    n_cells = length(x),
    contamination_median = stats::median(x),
    contamination_p90 = as.numeric(stats::quantile(x, 0.90, names = FALSE)),
    fraction_ge_0_10 = mean(x >= 0.10),
    fraction_ge_0_20 = mean(x >= 0.20),
    stringsAsFactors = FALSE
  )
}

cluster_ids <- f1_sort_cluster_ids(meta$epithelial_cluster_id)
ambient_cluster <- do.call(rbind, lapply(cluster_ids, function(cluster_id) {
  x <- meta$ambient_estimate[
    as.character(meta$epithelial_cluster_id) == cluster_id
  ]
  cbind(
    data.frame(epithelial_cluster = cluster_id, stringsAsFactors = FALSE),
    summarize_ambient(x)
  )
}))
f1_write_tsv(ambient_cluster, ambient_cluster_path)

sample_cluster_keys <- unique(data.frame(
  sample_id = as.character(meta$sample_id),
  epithelial_cluster = as.character(meta$epithelial_cluster_id),
  stringsAsFactors = FALSE
))
ambient_sample_cluster <- do.call(rbind, lapply(
  seq_len(nrow(sample_cluster_keys)),
  function(i) {
    sample_id <- sample_cluster_keys$sample_id[[i]]
    cluster_id <- sample_cluster_keys$epithelial_cluster[[i]]
    keep <- as.character(meta$sample_id) == sample_id &
      as.character(meta$epithelial_cluster_id) == cluster_id
    cbind(
      sample_cluster_keys[i, , drop = FALSE],
      summarize_ambient(meta$ambient_estimate[keep])
    )
  }
))
f1_write_tsv(ambient_sample_cluster, ambient_sample_cluster_path)

# 选择冻结的上皮亚型panel和每个新cluster的前30个marker。
marker_panel <- f1_read_tsv(config$paths$marker_panel)
f1_require_columns(
  marker_panel,
  c("cell_type", "positive_markers", "supporting_markers"),
  "cell_type_marker_panel.tsv"
)
epithelial_panel_names <- c(
  "Epithelial", "Pit_mucous_epithelial", "Mucous_neck_epithelial",
  "Chief_epithelial", "Parietal_epithelial",
  "Intestinal_like_epithelial", "Enteroendocrine_epithelial"
)
marker_panel <- marker_panel[
  marker_panel$cell_type %in% epithelial_panel_names,
  ,
  drop = FALSE
]
split_genes <- function(x) {
  x <- as.character(x)
  if (is.na(x) || !nzchar(x) || identical(x, "NA")) return(character())
  trimws(unlist(strsplit(x, ",", fixed = TRUE)))
}
panel_genes <- setNames(lapply(seq_len(nrow(marker_panel)), function(i) {
  unique(c(
    split_genes(marker_panel$positive_markers[[i]]),
    split_genes(marker_panel$supporting_markers[[i]])
  ))
}), marker_panel$cell_type)

markers <- f1_read_tsv(file.path(
  config$paths$annotation_dir,
  "epithelial_cluster_marker_genes.tsv"
))
f1_require_columns(
  markers,
  c("cluster", "gene", "p_val_adj"),
  "epithelial_cluster_marker_genes.tsv"
)
top <- f1_top_markers(markers, n = 30L)
top_genes <- setNames(
  lapply(top$top_markers, split_genes),
  as.character(top$cluster)
)
review_genes <- intersect(
  unique(c(unlist(panel_genes), unlist(top_genes))),
  rownames(raw)
)
if (!length(review_genes)) {
  stop("上皮panel和前列marker均无法匹配F1.5 RNA counts。")
}

panel_accumulator <- expand.grid(
  epithelial_cluster = cluster_ids,
  panel = names(panel_genes),
  stringsAsFactors = FALSE
)
panel_accumulator$raw_panel_counts <- 0
panel_accumulator$corrected_panel_counts <- 0
cluster_accumulator <- data.frame(
  epithelial_cluster = cluster_ids,
  raw_total_counts = 0,
  corrected_total_counts = 0,
  raw_top_marker_counts = 0,
  corrected_top_marker_counts = 0,
  raw_top_marker_detected_pairs = 0,
  corrected_top_marker_detected_pairs = 0,
  stringsAsFactors = FALSE
)

for (sample_id in sort(unique(as.character(meta$sample_id)))) {
  corrected_path <- file.path(
    config$paths$decontx_corrected_dir,
    paste0(sample_id, "_decontX_corrected_counts.rds")
  )
  if (!file.exists(corrected_path)) {
    stop("缺少DecontX corrected矩阵：", corrected_path)
  }
  corrected <- readRDS(corrected_path)
  sample_cells <- colnames(epithelial)[as.character(meta$sample_id) == sample_id]
  if (!all(sample_cells %in% colnames(corrected))) {
    stop(sample_id, "的corrected矩阵缺少F1.5上皮细胞。")
  }
  available_genes <- intersect(review_genes, rownames(corrected))
  raw_review <- raw[available_genes, sample_cells, drop = FALSE]
  corrected_review <- corrected[available_genes, sample_cells, drop = FALSE]
  raw_cell_totals <- Matrix::colSums(raw[, sample_cells, drop = FALSE])
  corrected_cell_totals <- Matrix::colSums(
    corrected[, sample_cells, drop = FALSE]
  )

  for (cluster_id in unique(as.character(
    meta$epithelial_cluster_id[match(sample_cells, rownames(meta))]
  ))) {
    cluster_cells <- sample_cells[
      as.character(meta$epithelial_cluster_id[match(sample_cells, rownames(meta))]) ==
        cluster_id
    ]
    cluster_row <- match(cluster_id, cluster_accumulator$epithelial_cluster)
    cluster_accumulator$raw_total_counts[[cluster_row]] <-
      cluster_accumulator$raw_total_counts[[cluster_row]] +
      sum(raw_cell_totals[cluster_cells])
    cluster_accumulator$corrected_total_counts[[cluster_row]] <-
      cluster_accumulator$corrected_total_counts[[cluster_row]] +
      sum(corrected_cell_totals[cluster_cells])

    cluster_top <- intersect(top_genes[[cluster_id]], available_genes)
    if (length(cluster_top)) {
      raw_top <- raw_review[cluster_top, cluster_cells, drop = FALSE]
      corrected_top <-
        corrected_review[cluster_top, cluster_cells, drop = FALSE]
      cluster_accumulator$raw_top_marker_counts[[cluster_row]] <-
        cluster_accumulator$raw_top_marker_counts[[cluster_row]] + sum(raw_top)
      cluster_accumulator$corrected_top_marker_counts[[cluster_row]] <-
        cluster_accumulator$corrected_top_marker_counts[[cluster_row]] +
        sum(corrected_top)
      cluster_accumulator$raw_top_marker_detected_pairs[[cluster_row]] <-
        cluster_accumulator$raw_top_marker_detected_pairs[[cluster_row]] +
        sum(raw_top > 0)
      cluster_accumulator$corrected_top_marker_detected_pairs[[cluster_row]] <-
        cluster_accumulator$corrected_top_marker_detected_pairs[[cluster_row]] +
        sum(corrected_top > 0)
    }

    for (panel_name in names(panel_genes)) {
      genes <- intersect(panel_genes[[panel_name]], available_genes)
      if (!length(genes)) next
      panel_row <- panel_accumulator$epithelial_cluster == cluster_id &
        panel_accumulator$panel == panel_name
      panel_accumulator$raw_panel_counts[panel_row] <-
        panel_accumulator$raw_panel_counts[panel_row] +
        sum(raw_review[genes, cluster_cells, drop = FALSE])
      panel_accumulator$corrected_panel_counts[panel_row] <-
        panel_accumulator$corrected_panel_counts[panel_row] +
        sum(corrected_review[genes, cluster_cells, drop = FALSE])
    }
  }
  rm(corrected, raw_review, corrected_review)
}

panel_accumulator <- merge(
  panel_accumulator,
  cluster_accumulator[
    ,
    c("epithelial_cluster", "raw_total_counts", "corrected_total_counts")
  ],
  by = "epithelial_cluster",
  all.x = TRUE,
  sort = FALSE
)
panel_accumulator$raw_panel_fraction <-
  panel_accumulator$raw_panel_counts / panel_accumulator$raw_total_counts
panel_accumulator$corrected_panel_fraction <-
  panel_accumulator$corrected_panel_counts /
  panel_accumulator$corrected_total_counts
panel_accumulator$panel_count_retention <- ifelse(
  panel_accumulator$raw_panel_counts > 0,
  panel_accumulator$corrected_panel_counts /
    panel_accumulator$raw_panel_counts,
  NA_real_
)
f1_write_tsv(panel_accumulator, panel_path)

choose_winner <- function(x, value_col) {
  x <- x[order(-x[[value_col]], x$panel), , drop = FALSE]
  as.character(x$panel[[1]])
}
winner_rows <- do.call(rbind, lapply(cluster_ids, function(cluster_id) {
  x <- panel_accumulator[
    panel_accumulator$epithelial_cluster == cluster_id,
    ,
    drop = FALSE
  ]
  data.frame(
    epithelial_cluster = cluster_id,
    raw_panel_winner = choose_winner(x, "raw_panel_fraction"),
    corrected_panel_winner = choose_winner(x, "corrected_panel_fraction"),
    stringsAsFactors = FALSE
  )
}))
winner_rows$panel_winner_unchanged <-
  winner_rows$raw_panel_winner == winner_rows$corrected_panel_winner

cluster_accumulator$top_marker_count_retention <- ifelse(
  cluster_accumulator$raw_top_marker_counts > 0,
  cluster_accumulator$corrected_top_marker_counts /
    cluster_accumulator$raw_top_marker_counts,
  NA_real_
)
cluster_accumulator$top_marker_detection_retention <- ifelse(
  cluster_accumulator$raw_top_marker_detected_pairs > 0,
  cluster_accumulator$corrected_top_marker_detected_pairs /
    cluster_accumulator$raw_top_marker_detected_pairs,
  NA_real_
)

impact <- merge(
  ambient_cluster,
  winner_rows,
  by = "epithelial_cluster",
  all = TRUE,
  sort = FALSE
)
impact <- merge(
  impact,
  cluster_accumulator[
    ,
    c(
      "epithelial_cluster", "top_marker_count_retention",
      "top_marker_detection_retention"
    )
  ],
  by = "epithelial_cluster",
  all = TRUE,
  sort = FALSE
)
impact$targeted_review_flag <-
  !impact$panel_winner_unchanged |
  impact$top_marker_count_retention < 0.80 |
  impact$contamination_median >= 0.10 |
  impact$contamination_p90 >= 0.50
impact$interpretation <- ifelse(
  !impact$panel_winner_unchanged |
    impact$top_marker_count_retention < 0.80,
  "review_raw_corrected_before_subtype_approval",
  ifelse(
    impact$targeted_review_flag,
    "marker_direction_stable_but_high_ambient_tail_review",
    "raw_marker_direction_supported"
  )
)
f1_write_tsv(impact, impact_path)

f1_save_session_info(config, "F1_05_validate_and_ambient_review")
f1_append_log(
  config,
  stage,
  sprintf(
    paste0(
      "完成：%d个上皮cluster；panel winner稳定%d/%d；",
      "前列marker count retention最低%.3f"
    ),
    nrow(impact),
    sum(impact$panel_winner_unchanged),
    nrow(impact),
    min(impact$top_marker_count_retention, na.rm = TRUE)
  )
)
message("F1.5输出验证与ambient复核完成：", impact_path)
