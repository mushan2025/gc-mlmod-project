# F1.4 DecontX输出验证 ------------------------------------------------------
#
# 目的：
#   确认正式cluster版DecontX已经覆盖全部样本和细胞，校正矩阵合法，
#   且向对象增加注释和污染估计后没有改变RNA raw counts、细胞顺序或降维结果。
# 本脚本只验证技术完整性，不根据污染分数删除细胞或改变注释。

file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), winslash = "/"))
source(file.path(script_dir, "F1_config.R"), encoding = "UTF-8")
source(file.path(script_dir, "F1_utils.R"), encoding = "UTF-8")

args <- f1_parse_args()
config <- f1_build_config(args$project_root)
f1_require_packages(c("SeuratObject", "Matrix", "data.table"), "F1.4 validation")

ambient_cells <- f1_read_tsv(config$paths$ambient_cell_estimates)
ambient_summary <- f1_read_tsv(config$paths$ambient_summary)
corrected_summary <- f1_read_tsv(
  config$paths$ambient_corrected_count_summary
)
expected_status <- "completed_with_researcher_approved_seurat_clusters"

stopifnot(
  nrow(ambient_summary) == 40L,
  nrow(corrected_summary) == 40L,
  nrow(ambient_cells) == 142650L,
  !anyDuplicated(ambient_cells$cell_id_final),
  all(ambient_summary$status == expected_status),
  all(corrected_summary$status == expected_status),
  all(ambient_cells$decontX_status == expected_status),
  all(
    ambient_cells$model_label_source ==
      "researcher_approved_seurat_cluster_partition"
  )
)

contamination <- ambient_cells$retained_cell_ambient_contamination_estimate
stopifnot(
  all(is.finite(contamination)),
  all(contamination >= 0 & contamination <= 1)
)

corrected_paths <- list.files(
  config$paths$decontx_corrected_dir,
  pattern = "_decontX_corrected_counts[.]rds$",
  full.names = TRUE
)
stopifnot(length(corrected_paths) == 40L)

# 两个大对象顺序读取并及时释放，避免在内存中同时保留两份。
object_03a <- readRDS(config$paths$object_03a)
object_03a <- f1_join_assay(object_03a, "RNA")
counts_03a <- SeuratObject::LayerData(
  object_03a,
  assay = "RNA",
  layer = "counts"
)
signature_03a <- f1_counts_signature(counts_03a)
cells_03a <- colnames(object_03a)
reductions_03a <- names(object_03a@reductions)
rm(counts_03a, object_03a)
invisible(gc())

object_03 <- readRDS(config$paths$object_03)
object_03 <- f1_join_assay(object_03, "RNA")
counts_03 <- SeuratObject::LayerData(
  object_03,
  assay = "RNA",
  layer = "counts"
)
signature_03 <- f1_counts_signature(counts_03)

stopifnot(
  identical(unname(signature_03a), unname(signature_03)),
  identical(cells_03a, colnames(object_03)),
  identical(reductions_03a, names(object_03@reductions)),
  all(
    c(
      "cell_type_major", "cell_type_minor", "cell_state",
      "annotation_confidence", "annotation_reason",
      "retained_cell_ambient_contamination_estimate",
      "decontX_evaluation_status"
    ) %in% colnames(object_03[[]])
  )
)

cell_index <- match(colnames(object_03), ambient_cells$cell_id_final)
stopifnot(
  !anyNA(cell_index),
  isTRUE(all.equal(
    as.numeric(object_03$retained_cell_ambient_contamination_estimate),
    as.numeric(contamination[cell_index]),
    tolerance = 1e-12,
    check.attributes = FALSE
  )),
  all(
    as.character(object_03$decontX_evaluation_status) ==
      ambient_cells$decontX_status[cell_index]
  ),
  all(
    as.character(object_03$cell_type_major) ==
      ambient_cells$coarse_lineage_label[cell_index]
  )
)

sample_ids <- sort(unique(as.character(object_03$sample_id)))
calculated_rows <- vector("list", length(sample_ids))
names(calculated_rows) <- sample_ids

for (sample_id in sample_ids) {
  path <- file.path(
    config$paths$decontx_corrected_dir,
    paste0(sample_id, "_decontX_corrected_counts.rds")
  )
  stopifnot(file.exists(path))
  corrected <- readRDS(path)
  sample_cells <- colnames(object_03)[
    as.character(object_03$sample_id) == sample_id
  ]
  stopifnot(
    !anyDuplicated(colnames(corrected)),
    setequal(colnames(corrected), sample_cells),
    all(rownames(corrected) %in% rownames(counts_03))
  )

  corrected_values <- if (inherits(corrected, "sparseMatrix")) {
    corrected@x
  } else {
    as.vector(corrected)
  }
  stopifnot(
    all(is.finite(corrected_values)),
    all(corrected_values >= 0)
  )

  raw_total <- sum(counts_03[, sample_cells, drop = FALSE])
  corrected_total <- sum(corrected_values)
  stopifnot(corrected_total <= raw_total + 1e-6)

  calculated_rows[[sample_id]] <- data.frame(
    sample_id = sample_id,
    input_cells = length(sample_cells),
    raw_total_counts = as.numeric(raw_total),
    corrected_total_counts = as.numeric(corrected_total),
    removed_count_fraction = as.numeric(
      (raw_total - corrected_total) / raw_total
    ),
    corrected_noninteger_fraction = if (length(corrected_values)) {
      mean(abs(corrected_values - round(corrected_values)) > 1e-8)
    } else {
      0
    },
    stringsAsFactors = FALSE
  )
  rm(corrected, corrected_values)
  invisible(gc())
}

calculated <- do.call(rbind, calculated_rows)
reported <- corrected_summary[
  match(calculated$sample_id, corrected_summary$sample_id),
  ,
  drop = FALSE
]
stopifnot(
  identical(calculated$sample_id, reported$sample_id),
  identical(as.integer(calculated$input_cells), as.integer(reported$input_cells)),
  isTRUE(all.equal(
    calculated$raw_total_counts,
    reported$raw_total_counts,
    tolerance = 1e-10,
    check.attributes = FALSE
  )),
  isTRUE(all.equal(
    calculated$corrected_total_counts,
    reported$corrected_total_counts,
    tolerance = 1e-10,
    check.attributes = FALSE
  )),
  isTRUE(all.equal(
    calculated$removed_count_fraction,
    reported$removed_count_fraction,
    tolerance = 1e-10,
    check.attributes = FALSE
  ))
)

validation <- data.frame(
  check = c(
    "samples_completed",
    "cells_with_finite_estimate",
    "corrected_matrices",
    "raw_counts_signature_preserved",
    "cell_order_preserved",
    "reductions_preserved",
    "annotation_and_estimate_metadata_match",
    "corrected_values_nonnegative",
    "corrected_totals_match_summary"
  ),
  status = "PASS",
  value = c(
    nrow(ambient_summary),
    nrow(ambient_cells),
    length(corrected_paths),
    paste(names(signature_03), signature_03, sep = "=", collapse = ";"),
    length(cells_03a),
    paste(reductions_03a, collapse = "|"),
    ncol(object_03),
    nrow(calculated),
    nrow(calculated)
  ),
  stringsAsFactors = FALSE
)
f1_write_tsv(
  validation,
  file.path(config$paths$qc_dir, "F1_DecontX_output_validation.tsv")
)

message("F1.4 DecontX输出验证通过。")
