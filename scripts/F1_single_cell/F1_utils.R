# F1 共用函数 ---------------------------------------------------------------
#
# 这些函数只处理六个脚本反复需要的读写、输入检查和Seurat公共步骤。
# 生物学判断仍保留在对应的F1分节脚本与人工审核表中。

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L || is.na(x[[1]]) || !nzchar(as.character(x[[1]]))) y else x
}

f1_script_dir <- function() {
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) != 1L) {
    stop("无法确定当前R脚本位置；请用Rscript运行，而不是在交互式R中直接粘贴。")
  }
  dirname(normalizePath(sub("^--file=", "", file_arg), winslash = "/", mustWork = TRUE))
}

f1_parse_args <- function(args = commandArgs(trailingOnly = TRUE)) {
  value_after <- function(prefix, default = NULL) {
    hit <- grep(paste0("^", prefix, "="), args, value = TRUE)
    if (!length(hit)) return(default)
    sub(paste0("^", prefix, "="), "", hit[[length(hit)]])
  }

  list(
    execute = "--execute" %in% args,
    approve_cnv_execution = "--approve-cnv-execution" %in% args,
    project_root = value_after("--project-root", getwd()),
    from = value_after("--from", "F1.1"),
    to = value_after("--to", "F1.6")
  )
}

f1_prepare_directories <- function(config) {
  dirs <- unlist(config$paths[c("object_dir", "qc_dir", "annotation_dir", "malignancy_dir", "log_dir")])
  invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))
}

f1_now <- function() format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")

f1_append_log <- function(config, stage, message) {
  f1_prepare_directories(config)
  if (!file.exists(config$paths$analysis_log)) {
    writeLines(c("# F1 Analysis Log", ""), config$paths$analysis_log, useBytes = TRUE)
  }
  cat(
    sprintf("- %s | %s | %s\n", f1_now(), stage, message),
    file = config$paths$analysis_log,
    append = TRUE
  )
}

f1_read_tsv <- function(path, required = TRUE) {
  if (!file.exists(path)) {
    if (required) stop("缺少必需文件：", path)
    return(NULL)
  }
  data.table::fread(path, sep = "\t", header = TRUE, data.table = FALSE, na.strings = c("NA"))
}

f1_write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(x, path, sep = "\t", quote = FALSE, na = "NA")
  invisible(path)
}

f1_require_columns <- function(x, columns, label) {
  missing <- setdiff(columns, colnames(x))
  if (length(missing)) {
    stop(label, "缺少字段：", paste(missing, collapse = ", "))
  }
  invisible(TRUE)
}

f1_as_logical <- function(x, label = "logical field") {
  y <- tolower(trimws(as.character(x)))
  if (any(!y %in% c("true", "false", "1", "0", "yes", "no"))) {
    stop(label, "含有无法解释的布尔值：", paste(unique(y[!y %in% c("true", "false", "1", "0", "yes", "no")]), collapse = ", "))
  }
  y %in% c("true", "1", "yes")
}

f1_stage_dry_run <- function(stage, description, inputs, outputs, packages) {
  cat("\n", stage, "：", description, "\n", sep = "")
  cat("模式：只检查计划，不读取表达矩阵、不写分析结果。\n")
  cat("输入：\n")
  for (x in inputs) cat("  - ", x, if (file.exists(x)) " [存在]" else " [尚不存在]", "\n", sep = "")
  cat("输出：\n")
  for (x in outputs) cat("  - ", x, "\n", sep = "")
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  cat("依赖：", if (length(missing)) paste0("缺少 ", paste(missing, collapse = ", ")) else "均可用", "\n", sep = "")
  invisible(length(missing) == 0L)
}

f1_require_packages <- function(packages, stage) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop(stage, "缺少必需R包：", paste(missing, collapse = ", "),
         "。本脚本不会自动安装包，请先按environment/F1/required_packages.tsv补齐。")
  }
  invisible(TRUE)
}

f1_package_versions <- function(packages) {
  data.frame(
    package = packages,
    version = vapply(packages, function(pkg) {
      if (requireNamespace(pkg, quietly = TRUE)) as.character(utils::packageVersion(pkg)) else "NOT_AVAILABLE"
    }, character(1)),
    stringsAsFactors = FALSE
  )
}

f1_save_session_info <- function(config, stage) {
  f1_prepare_directories(config)
  path <- file.path(config$paths$log_dir, paste0(stage, "_sessionInfo.txt"))
  capture.output(sessionInfo(), file = path)
  invisible(path)
}

f1_save_rds_atomic <- function(object, path, compress = FALSE) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temp <- paste0(path, ".tmp")
  on.exit(unlink(temp), add = TRUE)
  saveRDS(object, temp, compress = compress)
  if (file.exists(path)) unlink(path)
  if (!file.rename(temp, path)) stop("RDS临时文件无法改名为正式文件：", path)
  invisible(path)
}

f1_check_f0_ready <- function(config, stop_on_failure = TRUE) {
  required <- unlist(config$paths[c("processed_manifest", "sample_info", "data_audit", "f0_gate", "f0_report")])
  missing <- required[!file.exists(required)]
  problems <- character()
  if (length(missing)) {
    problems <- c(problems, paste0("缺少正式F0输出：", paste(missing, collapse = "; ")))
  } else {
    gate <- f1_read_tsv(config$paths$f0_gate)
    f1_require_columns(gate, c("pass_fail", "blocking_level"), "F0_gate_checklist.tsv")
    blocking_fail <- gate$blocking_level == "blocking" & gate$pass_fail == "FAIL"
    if (any(blocking_fail)) {
      problems <- c(problems, paste0("F0仍有blocking FAIL：", paste(gate$gate_item[blocking_fail], collapse = ", ")))
    }
    report <- readLines(config$paths$f0_report, warn = FALSE, encoding = "UTF-8")
    gate_line <- grep("F0_scRNA_F1_gate:", report, value = TRUE)
    if (!length(gate_line) || !any(grepl("F0_scRNA_F1_gate:\\s*(PASS|PASS_WITH_NOTED_LIMITATIONS)", gate_line))) {
      problems <- c(problems, "F0_execution_report.md未记录可进入F1的PASS状态")
    }

    sample_info <- f1_read_tsv(config$paths$sample_info)
    audit <- f1_read_tsv(config$paths$data_audit)
    f1_require_columns(sample_info, c("sample_id", "include_in_f1"), "sample_info.tsv")
    f1_require_columns(
      audit,
      c("sample_id", "audit_decision", "matrix_orientation", "matrix_orientation_validation_status"),
      "data_audit.tsv"
    )
    eligible <- sample_info$sample_id[f1_as_logical(sample_info$include_in_f1, "sample_info$include_in_f1")]
    audit_ok <- audit$sample_id[
      audit$audit_decision == "enter_full_F1_independent_reQC" &
        audit$matrix_orientation == "gene_by_cell" &
        audit$matrix_orientation_validation_status == "pass_gene_by_cell"
    ]
    if (!length(eligible) || !all(eligible %in% audit_ok)) {
      problems <- c(problems, "拟纳入F1的样本没有全部通过F0方向与输入边界审计")
    }
  }

  if (length(problems) && stop_on_failure) stop(paste(problems, collapse = "\n"))
  list(ready = !length(problems), problems = problems)
}

f1_assert_integer_counts <- function(counts, label) {
  values <- if (inherits(counts, "sparseMatrix")) counts@x else as.vector(counts)
  if (any(!is.finite(values))) stop(label, "含NA/Inf，不能作为raw counts。")
  if (any(values < 0)) stop(label, "含负值，不能作为raw counts。")
  if (any(abs(values - round(values)) > 1e-8)) stop(label, "含非整数值，不能作为raw integer counts。")
  invisible(TRUE)
}

f1_counts_signature <- function(counts) {
  values <- if (inherits(counts, "sparseMatrix")) counts@x else as.vector(counts)
  c(
    n_features = nrow(counts),
    n_cells = ncol(counts),
    nonzero = if (inherits(counts, "sparseMatrix")) length(counts@x) else sum(counts != 0),
    total_counts = sum(values)
  )
}

f1_join_assay <- function(object, assay) {
  if (!assay %in% names(object@assays)) return(object)
  layers <- SeuratObject::Layers(object[[assay]])
  if (length(layers) > 1L) object <- SeuratObject::JoinLayers(object, assay = assay)
  object
}

f1_get_counts <- function(object, assay = "RNA") {
  object <- f1_join_assay(object, assay)
  counts <- SeuratObject::LayerData(object, assay = assay, layer = "counts")
  list(object = object, counts = counts)
}

f1_summary_stats <- function(x, prefix) {
  q <- stats::quantile(x, probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE, names = FALSE, type = 7)
  out <- as.list(q)
  names(out) <- paste0(prefix, c("_min", "_Q1", "_median", "_Q3", "_max"))
  out
}

f1_collapse_reasons <- function(...) {
  flags <- list(...)
  labels <- names(flags)
  n <- length(flags[[1]])
  vapply(seq_len(n), function(i) {
    failed <- labels[!vapply(flags, function(x) isTRUE(x[[i]]), logical(1))]
    if (length(failed)) paste(failed, collapse = "|") else "none"
  }, character(1))
}

f1_top_markers <- function(markers, cluster_col = "cluster", n = 30L) {
  if (!nrow(markers)) return(data.frame(cluster = character(), top_markers = character()))
  fc_col <- intersect(c("avg_log2FC", "avg_logFC"), colnames(markers))[[1]]
  split_markers <- split(markers, as.character(markers[[cluster_col]]))
  rows <- lapply(names(split_markers), function(cluster_id) {
    x <- split_markers[[cluster_id]]
    x <- x[order(x$p_val_adj, -x[[fc_col]], x$gene), , drop = FALSE]
    data.frame(
      cluster = cluster_id,
      top_markers = paste(utils::head(unique(x$gene), n), collapse = ","),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

f1_run_sct_harmony <- function(object, config, stage_label) {
  set.seed(config$seed)
  object <- f1_join_assay(object, "RNA")
  before <- f1_counts_signature(SeuratObject::LayerData(object, assay = "RNA", layer = "counts"))

  old_future_plan <- future::plan("list")
  old_future_max_size <- getOption("future.globals.maxSize")
  on.exit(future::plan(old_future_plan), add = TRUE)
  on.exit(options(future.globals.maxSize = old_future_max_size), add = TRUE)
  options(
    future.globals.maxSize =
      config$execution$future_globals_max_gb * 1024^3
  )
  if (.Platform$OS.type != "windows" && config$execution$future_workers > 1L) {
    # Linux使用fork并行，主要加速Seurat可由future拆分的步骤；seed仍由各算法固定。
    future::plan(
      future::multicore,
      workers = config$execution$future_workers
    )
  } else {
    future::plan(future::sequential)
  }

  # Seurat v5按sample_id拆分RNA counts层后，SCTransform会为各样本分别拟合模型。
  object[["RNA"]] <- split(object[["RNA"]], f = factor(object$sample_id))
  SeuratObject::DefaultAssay(object) <- "RNA"
  object <- Seurat::SCTransform(
    object = object,
    assay = "RNA",
    new.assay.name = "SCT",
    vst.flavor = config$sct$vst_flavor,
    variable.features.n = config$sct$variable_features_n,
    vars.to.regress = config$sct$vars_to_regress,
    conserve.memory = config$sct$conserve_memory,
    return.only.var.genes = TRUE,
    seed.use = config$seed,
    verbose = TRUE
  )

  SeuratObject::DefaultAssay(object) <- "SCT"
  object <- Seurat::RunPCA(
    object,
    assay = "SCT",
    npcs = config$sct$pca_npcs,
    reduction.name = "pca",
    seed.use = config$seed,
    verbose = TRUE
  )
  object <- Seurat::RunUMAP(
    object,
    reduction = "pca",
    dims = config$sct$main_dims,
    reduction.name = "umap.sct_pca",
    reduction.key = "SCTPCAUMAP_",
    seed.use = config$seed,
    verbose = TRUE
  )
  object <- harmony::RunHarmony(
    object = object,
    group.by.vars = config$sct$harmony_group,
    reduction.use = "pca",
    dims.use = seq_len(config$sct$pca_npcs),
    reduction.save = "harmony",
    verbose = TRUE
  )
  object <- Seurat::RunUMAP(
    object,
    reduction = "harmony",
    dims = config$sct$main_dims,
    reduction.name = "umap",
    reduction.key = "UMAP_",
    seed.use = config$seed,
    verbose = TRUE
  )
  object <- Seurat::FindNeighbors(
    object,
    reduction = "harmony",
    dims = config$sct$main_dims,
    graph.name = c("SCT_harmony_nn", "SCT_harmony_snn"),
    verbose = TRUE
  )
  for (resolution in config$sct$resolutions) {
    cluster_name <- paste0("SCT_harmony_res.", format(resolution, trim = TRUE, scientific = FALSE))
    object <- Seurat::FindClusters(
      object,
      graph.name = "SCT_harmony_snn",
      resolution = resolution,
      algorithm = config$sct$leiden_algorithm,
      random.seed = config$seed,
      cluster.name = cluster_name,
      verbose = TRUE
    )
  }
  selected <- paste0(
    "SCT_harmony_res.",
    format(config$sct$default_resolution, trim = TRUE, scientific = FALSE)
  )
  if (!selected %in% colnames(object[[]])) stop(stage_label, "未生成默认resolution字段：", selected)
  object$seurat_clusters <- factor(object[[selected, drop = TRUE]])
  SeuratObject::Idents(object) <- "seurat_clusters"

  # marker使用RNA的常规log-normalized data层；这不会改变raw counts，也不参与主聚类。
  object <- f1_join_assay(object, "RNA")
  SeuratObject::DefaultAssay(object) <- "RNA"
  object <- Seurat::NormalizeData(
    object,
    assay = "RNA",
    normalization.method = "LogNormalize",
    scale.factor = 10000,
    verbose = TRUE
  )
  after <- f1_counts_signature(SeuratObject::LayerData(object, assay = "RNA", layer = "counts"))
  if (!identical(unname(before), unname(after))) {
    stop(stage_label, "运行SCTransform/Harmony后RNA raw counts签名发生变化，已停止。")
  }
  object
}

f1_write_parameter_versions <- function(config, extra = NULL) {
  packages <- unique(unlist(config$packages, use.names = FALSE))
  versions <- f1_package_versions(packages)
  versions$record_type <- "package_version"
  versions$parameter <- versions$package
  versions$value <- versions$version
  out <- versions[, c("record_type", "parameter", "value")]
  base_parameters <- data.frame(
    record_type = "analysis_parameter",
    parameter = c(
      "seed", "SCTransform_vst_flavor", "SCTransform_variable_features_n",
      "SCTransform_vars_to_regress", "PCA_npcs", "Harmony_group_by",
      "main_dims", "default_resolution", "UCell_F2_input",
      "SCTransform_and_UMAP_future_workers", "future_globals_max_GB",
      "scDblFinder_workers",
      "DecontX_timing", "DecontX_cluster_label_source",
      "DecontX_minimum_reliable_lineages", "DecontX_background",
      "inferCNV_analysis_mode", "inferCNV_internal_subclustering",
      "inferCNV_HMM", "inferCNV_output_format", "inferCNV_threads",
      "CopyKAT_cores",
      "CopyKAT_input_scope", "CopyKAT_known_normal_rule"
    ),
    value = c(
      config$seed, config$sct$vst_flavor, config$sct$variable_features_n,
      "NULL", config$sct$pca_npcs, config$sct$harmony_group,
      paste(range(config$sct$main_dims), collapse = ":"), config$sct$default_resolution,
      "RNA_raw_counts",
      config$execution$future_workers, config$execution$future_globals_max_gb,
      config$doublet$scdblfinder_workers,
      "after_researcher_approved_coarse_lineage_annotation",
      "cell_type_major", config$ambient$minimum_reliable_lineages, "NULL",
      config$cnv$infercnv_analysis_mode,
      identical(config$cnv$infercnv_analysis_mode, "subclusters"),
      config$cnv$infercnv_hmm, config$cnv$infercnv_output_format,
      config$cnv$infercnv_threads, config$cnv$copykat_cores,
      "all_QC_singlet_cells_from_current_sample_only",
      paste0(
        "same_sample_T_NK_then_B_Plasma_if_at_least_",
        config$cnv$minimum_reference_cells,
        "_else_automatic_within_sample"
      )
    ),
    stringsAsFactors = FALSE
  )
  out <- rbind(out, base_parameters)
  if (!is.null(extra)) out <- rbind(out, extra)
  f1_write_tsv(out, file.path(config$paths$annotation_dir, "F1_parameters_and_versions.tsv"))
}
