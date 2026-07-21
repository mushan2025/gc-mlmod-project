# F1.1 导入公开raw count矩阵 -------------------------------------------------
#
# 生物学目的：建立不丢失公开细胞和feature的初始Seurat对象，作为后续QC的
# 可追溯起点。本脚本只导入和核对，不过滤细胞，也不做标准化。
#
# 主要输入：F0正式生成的manifest、sample_info、data_audit及40个csv.gz矩阵。
# 主要输出：01_all_cells_raw_or_initial.rds、import_integrity_check.tsv。

file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), winslash = "/"))
source(file.path(script_dir, "F1_config.R"), encoding = "UTF-8")
source(file.path(script_dir, "F1_utils.R"), encoding = "UTF-8")

args <- f1_parse_args()
config <- f1_build_config(args$project_root)
stage <- "F1.1"
outputs <- c(
  config$paths$object_01,
  file.path(config$paths$qc_dir, "import_integrity_check.tsv")
)

if (!args$execute) {
  f1_stage_dry_run(
    stage,
    "逐样本导入gene-by-cell raw counts并建立完整初始对象",
    unlist(config$paths[c("processed_manifest", "sample_info", "data_audit", "f0_gate", "f0_report")]),
    outputs,
    config$packages[[stage]]
  )
  quit(save = "no", status = 0)
}

f1_require_packages(config$packages[[stage]], stage)
f1_check_f0_ready(config)
f1_prepare_directories(config)
set.seed(config$seed)
f1_append_log(config, stage, "开始导入F0批准的公开raw count矩阵")

manifest <- f1_read_tsv(config$paths$processed_manifest)
sample_info <- f1_read_tsv(config$paths$sample_info)
audit <- f1_read_tsv(config$paths$data_audit)
f1_require_columns(manifest, c("sample_id", "geo_accession", "extracted_path"), "processed_input_manifest.tsv")
f1_require_columns(
  sample_info,
  c("sample_id", "geo_accession", "patient_id", "group_analysis", "include_in_f1"),
  "sample_info.tsv"
)
f1_require_columns(
  audit,
  c(
    "sample_id", "file_name", "matrix_rows_genes", "matrix_cols_cells",
    "matrix_orientation", "matrix_orientation_validation_status", "audit_decision"
  ),
  "data_audit.tsv"
)

include_samples <- sample_info$sample_id[f1_as_logical(sample_info$include_in_f1, "include_in_f1")]
audit <- audit[
  audit$sample_id %in% include_samples &
    audit$audit_decision == "enter_full_F1_independent_reQC",
  , drop = FALSE
]
if (!nrow(audit) || !setequal(include_samples, audit$sample_id)) {
  stop("sample_info拟纳入样本与data_audit允许进入F1的样本不完全一致。")
}
if (anyDuplicated(manifest$sample_id) || anyDuplicated(sample_info$sample_id) || anyDuplicated(audit$sample_id)) {
  stop("F0登记表中sample_id不是一对一，不能安全导入。")
}

objects <- vector("list", length(include_samples))
names(objects) <- include_samples
integrity_rows <- vector("list", length(include_samples))

for (i in seq_along(include_samples)) {
  sample_id <- include_samples[[i]]
  message(sprintf("[%s] 导入 %s (%d/%d)", stage, sample_id, i, length(include_samples)))
  m <- manifest[manifest$sample_id == sample_id, , drop = FALSE]
  s <- sample_info[sample_info$sample_id == sample_id, , drop = FALSE]
  a <- audit[audit$sample_id == sample_id, , drop = FALSE]
  if (nrow(m) != 1L || nrow(s) != 1L || nrow(a) != 1L) {
    stop(sample_id, "在manifest/sample_info/data_audit中不是唯一一行。")
  }
  if (a$matrix_orientation != "gene_by_cell" || a$matrix_orientation_validation_status != "pass_gene_by_cell") {
    stop(sample_id, "未通过gene-by-cell方向验证；脚本不会自动转置。")
  }

  matrix_path <- as.character(m$extracted_path)
  if (!grepl("^[A-Za-z]:[/\\\\]", matrix_path) && !startsWith(matrix_path, "/")) {
    matrix_path <- file.path(config$project_root, matrix_path)
  }
  matrix_path <- normalizePath(matrix_path, winslash = "/", mustWork = TRUE)

  # 一次只把一个样本读成稠密临时表，随后立即转稀疏矩阵并释放内存。
  tab <- data.table::fread(
    matrix_path,
    sep = ",",
    header = TRUE,
    check.names = FALSE,
    data.table = TRUE,
    showProgress = TRUE
  )
  if (ncol(tab) < 2L) stop(sample_id, "矩阵没有细胞列。")
  genes <- as.character(tab[[1]])
  barcodes <- colnames(tab)[-1]

  duplicate_genes <- unique(genes[duplicated(genes)])
  duplicate_barcodes <- unique(barcodes[duplicated(barcodes)])
  if (length(duplicate_genes) || length(duplicate_barcodes)) {
    issue_parts <- list()
    if (length(duplicate_genes)) {
      issue_parts[["gene"]] <- data.frame(
        sample_id = sample_id, id_type = "gene", duplicated_id = duplicate_genes
      )
    }
    if (length(duplicate_barcodes)) {
      issue_parts[["barcode"]] <- data.frame(
        sample_id = sample_id, id_type = "barcode", duplicated_id = duplicate_barcodes
      )
    }
    issue <- do.call(rbind, issue_parts)
    f1_write_tsv(issue, file.path(config$paths$qc_dir, "duplicate_feature_barcode_report.tsv"))
    stop(sample_id, "存在重复gene或barcode；已输出报告，未使用make.unique()。")
  }
  expected_rows <- as.integer(a$matrix_rows_genes)
  expected_cells <- as.integer(a$matrix_cols_cells)
  if (length(genes) != expected_rows || length(barcodes) != expected_cells) {
    stop(
      sample_id, "实际维度与F0不一致：实际", length(genes), " × ", length(barcodes),
      "；F0记录", expected_rows, " × ", expected_cells, "。"
    )
  }

  dense_counts <- as.matrix(tab[, -1, with = FALSE])
  rm(tab)
  rownames(dense_counts) <- genes
  colnames(dense_counts) <- barcodes
  f1_assert_integer_counts(dense_counts, paste0(sample_id, " raw matrix"))
  sparse_counts <- methods::as(Matrix::Matrix(dense_counts, sparse = TRUE), "dgCMatrix")
  rm(dense_counts)
  gc(verbose = FALSE)

  original_barcodes <- colnames(sparse_counts)
  final_ids <- paste0(sample_id, "__", original_barcodes)
  colnames(sparse_counts) <- final_ids
  object <- Seurat::CreateSeuratObject(
    counts = sparse_counts,
    assay = "RNA",
    project = "GSE183904",
    min.cells = 0,
    min.features = 0
  )
  object$cell_id_final <- final_ids
  object$original_barcode <- original_barcodes
  for (field in colnames(s)) {
    object[[field]] <- rep(as.character(s[[field]]), ncol(object))
  }
  integrity_rows[[i]] <- data.frame(
    sample_id = sample_id,
    geo_accession = as.character(s$geo_accession),
    matrix_path = matrix_path,
    expected_features = expected_rows,
    actual_features = nrow(object),
    expected_cells = expected_cells,
    actual_cells = ncol(object),
    matrix_orientation = "gene_by_cell",
    duplicate_gene_count = 0L,
    duplicate_barcode_count = 0L,
    metadata_match = all(object$sample_id == sample_id),
    raw_counts_integer_nonnegative = TRUE,
    import_status = "PASS",
    stringsAsFactors = FALSE
  )
  objects[[i]] <- object
  rm(object, sparse_counts)
  gc(verbose = FALSE)
}

# merge只合并稀疏raw counts和metadata；随后把Seurat v5分层counts接回一个RNA counts层。
initial <- objects[[1]]
if (length(objects) > 1L) {
  initial <- merge(initial, y = objects[-1], merge.data = FALSE)
}
rm(objects)
initial <- f1_join_assay(initial, "RNA")
counts <- SeuratObject::LayerData(initial, assay = "RNA", layer = "counts")
f1_assert_integer_counts(counts, "merged F1.1 RNA counts")
if (anyDuplicated(colnames(initial))) stop("合并后的cell_id_final仍有重复。")
if (ncol(initial) != sum(vapply(integrity_rows, function(x) x$actual_cells, integer(1)))) {
  stop("合并对象细胞数不等于逐样本细胞数之和。")
}

integrity <- do.call(rbind, integrity_rows)
f1_write_tsv(integrity, file.path(config$paths$qc_dir, "import_integrity_check.tsv"))
f1_save_rds_atomic(initial, config$paths$object_01, compress = FALSE)
f1_save_session_info(config, "F1_01_import")
f1_append_log(
  config,
  stage,
  sprintf("完成导入：%d个样本，%d个公开细胞，%d个feature；未过滤细胞", nrow(integrity), ncol(initial), nrow(initial))
)
message("F1.1完成：", config$paths$object_01)
