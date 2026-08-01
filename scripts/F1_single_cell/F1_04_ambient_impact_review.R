# F1.4 ambient RNA影响复核 -------------------------------------------------
#
# 生物学目的：
#   DecontX分数本身不能决定删细胞或换矩阵。本脚本比较raw与corrected中
#   经典谱系marker及原始cluster top markers是否仍保留，用少量定量指标判断
#   retained-cell ambient估计是否足以改变粗谱系注释或上皮候选范围。
#
# 输出只用于F1.4/F1.5的敏感性判断。raw counts继续作为主矩阵；
# corrected浮点值不用于DESeq2，也不反向修改已批准注释。

file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), winslash = "/"))
source(file.path(script_dir, "F1_config.R"), encoding = "UTF-8")
source(file.path(script_dir, "F1_utils.R"), encoding = "UTF-8")

args <- f1_parse_args()
config <- f1_build_config(args$project_root)
f1_require_packages(
  c("SeuratObject", "Matrix", "data.table"),
  "F1.4 ambient impact review"
)

panel_table <- f1_read_tsv(config$paths$marker_panel)
approved <- f1_read_tsv(config$paths$annotation_approved)
ambient_cells <- f1_read_tsv(config$paths$ambient_cell_estimates)
f1_require_columns(
  panel_table,
  c("cell_type", "positive_markers"),
  "cell_type_marker_panel.tsv"
)
f1_require_columns(
  approved,
  c(
    "cluster", "cell_type_major", "cell_type_minor", "cell_state",
    "top_markers"
  ),
  "F1_cluster_annotation_approved.tsv"
)

parse_markers <- function(cell_types) {
  values <- panel_table$positive_markers[
    panel_table$cell_type %in% cell_types
  ]
  values <- values[!is.na(values) & nzchar(values)]
  unique(trimws(unlist(strsplit(values, ",", fixed = TRUE))))
}

# Mesothelial没有写回冻结panel；这里只复用c22审核时已登记的marker作影响复核。
broad_panels <- list(
  "Epithelial" = parse_markers("Epithelial"),
  "T/NK" = parse_markers(c("T_cell", "NK_cell")),
  "B/Plasma" = parse_markers(c("B_cell", "Plasma_cell")),
  "Myeloid" = parse_markers(
    c("Myeloid", "Monocyte", "Macrophage", "Dendritic_cell")
  ),
  "Fibroblast/CAF" = parse_markers(c("Fibroblast", "CAF")),
  "Endothelial/Pericyte" = parse_markers(
    c("Endothelial", "Pericyte", "Smooth_muscle_cell")
  ),
  "Mast" = parse_markers("Mast_cell"),
  "Mesothelial" = c(
    "UPK3B", "LRRN4", "MSLN", "WT1", "PRG4", "CALB2", "ITLN1"
  )
)
stopifnot(
  setequal(names(broad_panels), unique(approved$cell_type_major)),
  all(lengths(broad_panels) > 0)
)

approved$cluster <- as.character(approved$cluster)
top_markers <- setNames(
  lapply(approved$top_markers, function(x) {
    unique(trimws(strsplit(x, ",", fixed = TRUE)[[1]]))
  }),
  approved$cluster
)

object <- readRDS(config$paths$object_03)
object <- f1_join_assay(object, "RNA")
raw_counts <- SeuratObject::LayerData(object, assay = "RNA", layer = "counts")
stopifnot(
  ncol(object) == nrow(ambient_cells),
  setequal(colnames(object), ambient_cells$cell_id_final)
)
raw_cell_totals <- Matrix::colSums(raw_counts)
review_genes <- unique(c(
  unlist(broad_panels, use.names = FALSE),
  unlist(top_markers, use.names = FALSE)
))
review_genes <- intersect(review_genes, rownames(raw_counts))

marker_metrics <- function(matrix, genes, cells, total_counts) {
  present <- intersect(genes, rownames(matrix))
  marker_counts <- if (length(present) && length(cells)) {
    matrix[present, cells, drop = FALSE]
  } else {
    NULL
  }
  detected <- if (is.null(marker_counts)) 0 else sum(marker_counts > 0)
  count_sum <- if (is.null(marker_counts)) 0 else sum(marker_counts)
  data.frame(
    marker_genes_expected = length(genes),
    marker_genes_present = length(present),
    marker_cell_pairs = length(genes) * length(cells),
    marker_detected_pairs = as.numeric(detected),
    marker_count_sum = as.numeric(count_sum),
    marker_count_fraction = if (total_counts > 0) {
      as.numeric(count_sum / total_counts)
    } else {
      NA_real_
    },
    stringsAsFactors = FALSE
  )
}

sample_ids <- sort(unique(as.character(object$sample_id)))
panel_rows <- list()
top_rows <- list()
panel_index <- 0L
top_index <- 0L
ambient_dt <- data.table::as.data.table(ambient_cells)

for (sid in sample_ids) {
  message("复核 ", sid, "（", match(sid, sample_ids), "/", length(sample_ids), "）")
  sample_cells <- colnames(object)[as.character(object$sample_id) == sid]
  corrected_path <- file.path(
    config$paths$decontx_corrected_dir,
    paste0(sid, "_decontX_corrected_counts.rds")
  )
  if (!file.exists(corrected_path)) {
    stop("缺少", sid, "的DecontX corrected矩阵：", corrected_path)
  }
  corrected <- readRDS(corrected_path)
  if (!setequal(colnames(corrected), sample_cells)) {
    stop(sid, "的corrected矩阵cell ID与03对象不一致。")
  }
  corrected <- corrected[, sample_cells, drop = FALSE]
  corrected_cell_totals <- Matrix::colSums(corrected)
  raw_review <- raw_counts[review_genes, sample_cells, drop = FALSE]
  corrected_review <- corrected[
    intersect(review_genes, rownames(corrected)),
    sample_cells,
    drop = FALSE
  ]
  sample_clusters <- as.character(
    object$seurat_clusters[match(sample_cells, colnames(object))]
  )

  for (cluster_id in f1_sort_cluster_ids(unique(sample_clusters))) {
    cluster_cells <- sample_cells[sample_clusters == cluster_id]
    raw_total <- sum(raw_cell_totals[cluster_cells])
    corrected_total <- sum(corrected_cell_totals[cluster_cells])
    annotation_row <- approved[approved$cluster == cluster_id, , drop = FALSE]
    if (nrow(annotation_row) != 1L) {
      stop("cluster ", cluster_id, "无法一对一匹配批准注释。")
    }

    for (panel_name in names(broad_panels)) {
      raw_metric <- marker_metrics(
        raw_review,
        broad_panels[[panel_name]],
        cluster_cells,
        raw_total
      )
      corrected_metric <- marker_metrics(
        corrected_review,
        broad_panels[[panel_name]],
        cluster_cells,
        corrected_total
      )
      panel_index <- panel_index + 1L
      panel_rows[[panel_index]] <- data.frame(
        sample_id = sid,
        cluster = cluster_id,
        cell_type_major = annotation_row$cell_type_major,
        panel_name = panel_name,
        is_expected_panel = panel_name == annotation_row$cell_type_major,
        n_cells = length(cluster_cells),
        raw_total_counts = as.numeric(raw_total),
        corrected_total_counts = as.numeric(corrected_total),
        raw_marker_genes_present = raw_metric$marker_genes_present,
        corrected_marker_genes_present =
          corrected_metric$marker_genes_present,
        marker_genes_expected = raw_metric$marker_genes_expected,
        marker_cell_pairs = raw_metric$marker_cell_pairs,
        raw_marker_detected_pairs = raw_metric$marker_detected_pairs,
        corrected_marker_detected_pairs =
          corrected_metric$marker_detected_pairs,
        raw_marker_count_sum = raw_metric$marker_count_sum,
        corrected_marker_count_sum = corrected_metric$marker_count_sum,
        raw_marker_count_fraction = raw_metric$marker_count_fraction,
        corrected_marker_count_fraction =
          corrected_metric$marker_count_fraction,
        stringsAsFactors = FALSE
      )
    }

    raw_top <- marker_metrics(
      raw_review,
      top_markers[[cluster_id]],
      cluster_cells,
      raw_total
    )
    corrected_top <- marker_metrics(
      corrected_review,
      top_markers[[cluster_id]],
      cluster_cells,
      corrected_total
    )
    cell_ambient <- ambient_dt[
      sample_id == sid &
        as.character(seurat_cluster) == cluster_id,
      retained_cell_ambient_contamination_estimate
    ]
    top_index <- top_index + 1L
    top_rows[[top_index]] <- data.frame(
      sample_id = sid,
      cluster = cluster_id,
      cell_type_major = annotation_row$cell_type_major,
      cell_type_minor = annotation_row$cell_type_minor,
      cell_state = annotation_row$cell_state,
      n_cells = length(cluster_cells),
      contamination_median = stats::median(cell_ambient),
      contamination_P90 = as.numeric(stats::quantile(
        cell_ambient,
        probs = 0.90,
        names = FALSE,
        type = 7
      )),
      raw_total_counts = as.numeric(raw_total),
      corrected_total_counts = as.numeric(corrected_total),
      top_marker_genes_expected = raw_top$marker_genes_expected,
      raw_top_marker_detected_pairs = raw_top$marker_detected_pairs,
      corrected_top_marker_detected_pairs =
        corrected_top$marker_detected_pairs,
      top_marker_cell_pairs = raw_top$marker_cell_pairs,
      raw_top_marker_count_sum = raw_top$marker_count_sum,
      corrected_top_marker_count_sum = corrected_top$marker_count_sum,
      stringsAsFactors = FALSE
    )
  }
  rm(
    corrected, corrected_cell_totals, raw_review, corrected_review
  )
  invisible(gc())
}

panel_sample <- data.table::rbindlist(panel_rows)
top_sample <- data.table::rbindlist(top_rows)

panel_cluster <- panel_sample[
  ,
  .(
    cell_type_major = unique(cell_type_major),
    n_cells = sum(n_cells),
    marker_genes_expected = unique(marker_genes_expected),
    raw_marker_detection = sum(raw_marker_detected_pairs) /
      sum(marker_cell_pairs),
    corrected_marker_detection = sum(corrected_marker_detected_pairs) /
      sum(marker_cell_pairs),
    raw_marker_count_fraction = sum(raw_marker_count_sum) /
      sum(raw_total_counts),
    corrected_marker_count_fraction = sum(corrected_marker_count_sum) /
      sum(corrected_total_counts),
    marker_count_retention = if (sum(raw_marker_count_sum) > 0) {
      sum(corrected_marker_count_sum) / sum(raw_marker_count_sum)
    } else {
      NA_real_
    }
  ),
  by = .(cluster, panel_name, is_expected_panel)
]
panel_cluster[
  ,
  marker_detection_retention := ifelse(
    raw_marker_detection > 0,
    corrected_marker_detection / raw_marker_detection,
    NA_real_
  )
]

top_cluster <- top_sample[
  ,
  .(
    cell_type_major = unique(cell_type_major),
    cell_type_minor = unique(cell_type_minor),
    cell_state = unique(cell_state),
    n_cells = sum(n_cells),
    n_samples = data.table::uniqueN(sample_id),
    raw_top_marker_detection = sum(raw_top_marker_detected_pairs) /
      sum(top_marker_cell_pairs),
    corrected_top_marker_detection =
      sum(corrected_top_marker_detected_pairs) / sum(top_marker_cell_pairs),
    top_marker_count_retention = if (sum(raw_top_marker_count_sum) > 0) {
      sum(corrected_top_marker_count_sum) / sum(raw_top_marker_count_sum)
    } else {
      NA_real_
    },
    total_count_retention = sum(corrected_total_counts) /
      sum(raw_total_counts)
  ),
  by = cluster
]
ambient_cluster <- ambient_dt[
  ,
  .(
    contamination_median = stats::median(
      retained_cell_ambient_contamination_estimate
    ),
    contamination_P90 = as.numeric(stats::quantile(
      retained_cell_ambient_contamination_estimate,
      probs = 0.90,
      names = FALSE,
      type = 7
    ))
  ),
  by = .(cluster = as.character(seurat_cluster))
]
top_cluster <- merge(
  top_cluster,
  ambient_cluster,
  by = "cluster",
  all.x = TRUE
)
top_cluster[
  ,
  top_marker_detection_retention := ifelse(
    raw_top_marker_detection > 0,
    corrected_top_marker_detection / raw_top_marker_detection,
    NA_real_
  )
]

winner_table <- panel_cluster[
  ,
  .(
    raw_panel_winner = panel_name[which.max(raw_marker_detection)],
    corrected_panel_winner =
      panel_name[which.max(corrected_marker_detection)]
  ),
  by = cluster
]
expected_panel <- panel_cluster[is_expected_panel == TRUE]
impact_summary <- merge(
  top_cluster,
  expected_panel[
    ,
    .(
      cluster,
      expected_panel = panel_name,
      raw_expected_panel_detection = raw_marker_detection,
      corrected_expected_panel_detection = corrected_marker_detection,
      expected_panel_detection_retention = marker_detection_retention,
      expected_panel_count_retention = marker_count_retention
    )
  ],
  by = "cluster",
  all.x = TRUE
)
impact_summary <- merge(
  impact_summary,
  winner_table,
  by = "cluster",
  all.x = TRUE
)
impact_summary[
  ,
  `:=`(
    raw_expected_panel_is_winner = expected_panel == raw_panel_winner,
    corrected_expected_panel_is_winner =
      expected_panel == corrected_panel_winner,
    panel_winner_changed = raw_panel_winner != corrected_panel_winner
  )
]
impact_summary[, cluster_order := as.integer(cluster)]
data.table::setorder(impact_summary, cluster_order)
impact_summary[, cluster_order := NULL]

top_sample[
  ,
  `:=`(
    top_marker_detection_retention = ifelse(
      raw_top_marker_detected_pairs > 0,
      corrected_top_marker_detected_pairs /
        raw_top_marker_detected_pairs,
      NA_real_
    ),
    top_marker_count_retention = ifelse(
      raw_top_marker_count_sum > 0,
      corrected_top_marker_count_sum / raw_top_marker_count_sum,
      NA_real_
    ),
    total_count_retention = ifelse(
      raw_total_counts > 0,
      corrected_total_counts / raw_total_counts,
      NA_real_
    )
  )
]
data.table::setorder(top_sample, -contamination_median)

f1_write_tsv(
  panel_cluster,
  file.path(
    config$paths$qc_dir,
    "ambient_marker_panel_raw_corrected_by_cluster.tsv"
  )
)
f1_write_tsv(
  top_sample,
  file.path(
    config$paths$qc_dir,
    "ambient_top_marker_retention_by_sample_cluster.tsv"
  )
)
f1_write_tsv(
  impact_summary,
  file.path(
    config$paths$qc_dir,
    "ambient_annotation_impact_summary.tsv"
  )
)

message("F1.4 ambient RNA注释影响复核完成。")
