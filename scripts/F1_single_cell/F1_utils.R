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

f1_sort_cluster_ids <- function(cluster_ids) {
  cluster_ids <- unique(as.character(cluster_ids))
  numeric_ids <- suppressWarnings(as.numeric(cluster_ids))
  if (all(is.finite(numeric_ids))) {
    cluster_ids[order(numeric_ids)]
  } else {
    sort(cluster_ids)
  }
}

f1_dominant_category <- function(values) {
  values <- as.character(values)
  valid <- !is.na(values) & nzchar(values)
  if (!any(valid)) {
    return(list(id = NA_character_, fraction = NA_real_))
  }

  counts <- table(values[valid])
  dominant_ids <- sort(names(counts)[counts == max(counts)])
  dominant_id <- dominant_ids[[1]]
  list(
    id = dominant_id,
    # 分母使用cluster全部细胞，避免有缺失patient_id时把占比人为放大。
    fraction = as.numeric(counts[[dominant_id]]) / length(values)
  )
}

f1_cluster_composition <- function(metadata) {
  f1_require_columns(
    metadata,
    c("seurat_clusters", "sample_id", "patient_id", "group_analysis"),
    "cluster composition metadata"
  )

  cluster_ids <- f1_sort_cluster_ids(metadata$seurat_clusters)
  group_levels <- c(
    "Normal_Gastric", "Primary_Tumor",
    "Peritoneal_Metastasis", "Normal_Peritoneum"
  )

  rows <- lapply(cluster_ids, function(cluster_id) {
    x <- metadata[
      as.character(metadata$seurat_clusters) == cluster_id,
      ,
      drop = FALSE
    ]
    dominant_sample <- f1_dominant_category(x$sample_id)
    dominant_patient <- f1_dominant_category(x$patient_id)
    valid_samples <- !is.na(x$sample_id) & nzchar(as.character(x$sample_id))
    valid_patients <- !is.na(x$patient_id) & nzchar(as.character(x$patient_id))

    out <- data.frame(
      cluster = cluster_id,
      cell_count = nrow(x),
      sample_count = length(unique(as.character(x$sample_id[valid_samples]))),
      n_patients = length(unique(as.character(x$patient_id[valid_patients]))),
      dominant_sample_id = dominant_sample$id,
      dominant_sample_fraction = dominant_sample$fraction,
      dominant_patient_id = dominant_patient$id,
      dominant_patient_fraction = dominant_patient$fraction,
      stringsAsFactors = FALSE
    )

    for (group_name in group_levels) {
      group_count <- sum(as.character(x$group_analysis) == group_name, na.rm = TRUE)
      count_name <- paste0(group_name, "_count")
      fraction_name <- paste0(group_name, "_fraction")
      out[[count_name]] <- group_count
      out[[fraction_name]] <- group_count / nrow(x)
    }
    out
  })

  do.call(rbind, rows)
}

f1_cluster_qc_metrics <- function(metadata) {
  f1_require_columns(
    metadata,
    c(
      "seurat_clusters", "nCount_RNA", "nFeature_RNA", "mt_percent",
      "HB_percent", "DoubletFinder_class", "scDblFinder_class"
    ),
    "cluster QC metadata"
  )

  cluster_ids <- f1_sort_cluster_ids(metadata$seurat_clusters)
  rows <- lapply(cluster_ids, function(cluster_id) {
    x <- metadata[
      as.character(metadata$seurat_clusters) == cluster_id,
      ,
      drop = FALSE
    ]
    data.frame(
      cluster = cluster_id,
      median_nCount_RNA = stats::median(x$nCount_RNA, na.rm = TRUE),
      median_nFeature_RNA = stats::median(x$nFeature_RNA, na.rm = TRUE),
      median_mt_percent = stats::median(x$mt_percent, na.rm = TRUE),
      median_HB_percent = stats::median(x$HB_percent, na.rm = TRUE),
      DoubletFinder_only_fraction = mean(
        tolower(x$DoubletFinder_class) == "doublet" &
          tolower(x$scDblFinder_class) == "singlet",
        na.rm = TRUE
      ),
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, rows)
}

f1_apply_f1_annotation_draft <- function(template) {
  required <- c(
    "cluster", "cell_type_major", "cell_type_minor", "cell_state",
    "annotation_confidence", "annotation_reason",
    "annotation_review_status", "downstream_handling_before_full_approval"
  )
  f1_require_columns(template, required, "F1 annotation review template")

  expected_clusters <- as.character(seq_len(24L))
  if (!setequal(as.character(template$cluster), expected_clusters)) {
    stop(
      "当前注释草案只适用于已审核的resolution 0.6、24-cluster结果；",
      "cluster集合改变后必须重新审核marker，不能沿用旧标签。"
    )
  }

  # major是DecontX和后续对象提取使用的粗谱系；minor/state只负责更细的生物学解释。
  major <- c(
    "1" = "T/NK", "2" = "T/NK", "3" = "Epithelial",
    "4" = "B/Plasma", "5" = "Fibroblast/CAF",
    "6" = "Endothelial/Pericyte", "7" = "T/NK", "8" = "B/Plasma",
    "9" = "Epithelial", "10" = "T/NK", "11" = "Myeloid",
    "12" = "B/Plasma", "13" = "Myeloid", "14" = "B/Plasma",
    "15" = "Endothelial/Pericyte", "16" = "Epithelial",
    "17" = "Mast", "18" = "B/Plasma", "19" = "Epithelial",
    "20" = "Myeloid", "21" = "Epithelial", "22" = "Mesothelial",
    "23" = "T/NK", "24" = "T/NK"
  )
  minor <- c(
    "1" = "CD8_GZMK_cytotoxic_T",
    "2" = "CCR6_IL7R_CD4_like_T",
    "3" = "Gastric_or_intestinal_like_epithelial_malignancy_pending_CNV",
    "4" = "JCHAIN_TNFRSF17_plasma",
    "5" = "LUM_DCN_PDGFRA_fibroblast",
    "6" = "Vascular_endothelial",
    "7" = "GNLY_NK_like_cytotoxic_lymphocyte",
    "8" = "IGHG_high_plasma",
    "9" = "Pit_surface_mucous_like_epithelial",
    "10" = "IFNG_CCL4_CD8_cytotoxic_T",
    "11" = "C1Q_FOLR2_TREM2_macrophage",
    "12" = "Conventional_B_cell",
    "13" = "FCN1_S100A8_A9_inflammatory_monocyte",
    "14" = "IGHA_high_plasma",
    "15" = "Pericyte_smooth_muscle_like",
    "16" = "Mucous_neck_chief_like_epithelial",
    "17" = "Mast_cell",
    "18" = "B_or_plasma_like_stress_high",
    "19" = "Chief_like_epithelial",
    "20" = "CD1C_FCER1A_dendritic_cell",
    "21" = "Enteroendocrine_epithelial",
    "22" = "Mesothelial",
    "23" = "CCR6_KLRB1_CCL20_CD4_like_T",
    "24" = "Cycling_T_cell"
  )
  state <- c(
    "1" = "GZMK_high_cytotoxic_memory_like",
    "2" = "activated_costimulatory_checkpoint_high",
    "3" = "mixed_gastric_and_intestinal_program",
    "4" = "antibody_secreting",
    "5" = "CCL11_CXCL14_stromal_program",
    "6" = "PLVAP_EMCN_high",
    "7" = "GNLY_GZMB_high_cytotoxic",
    "8" = "IGHG_high_antibody_secreting",
    "9" = "MUC5AC_CYP3A5_CLDN4_high",
    "10" = "IFNG_CCL4_TNF_high_activated",
    "11" = "C1Q_FOLR2_TREM2_TAM_like",
    "12" = "B_cell",
    "13" = "inflammatory",
    "14" = "IGHA_high_antibody_secreting",
    "15" = "contractile",
    "16" = "MUC6_TFF2_glandular_secretory",
    "17" = "mast_cell",
    "18" = "stress_high_likely_dissociation_associated",
    "19" = "metallothionein_high",
    "20" = "antigen_presentation_LAMP3_IDO1",
    "21" = "neuroendocrine_secretory",
    "22" = "mesothelial_ECM_program",
    "23" = "Th17_like_CCL20_high",
    "24" = "cycling"
  )
  confidence <- c(
    "1" = "high", "2" = "medium", "3" = "medium", "4" = "high",
    "5" = "high", "6" = "high", "7" = "medium", "8" = "high",
    "9" = "high", "10" = "high", "11" = "high", "12" = "high",
    "13" = "high", "14" = "high", "15" = "high", "16" = "medium",
    "17" = "high", "18" = "medium", "19" = "medium", "20" = "high",
    "21" = "high", "22" = "medium", "23" = "medium", "24" = "medium"
  )
  reason <- c(
    "1" = paste0(
      "CD3D/CD3E/CD3G与CD8A/CD8B成套表达，GZMK/CCL5/NKG7/GZMB提示",
      "细胞毒和效应记忆样程序；无上皮、B细胞或基质谱系冲突。"
    ),
    "2" = paste0(
      "CD3D/CD3E/CD3G支持T细胞，IL7R/CCR6/KLRB1及CTLA4/TNFRSF4/TIGIT",
      "提示活化CD4样程序；FOXP3和IL2RA未形成marker证据，因此不命名为Treg。"
    ),
    "3" = paste0(
      "GKN1/GKN2/TFF1胃表面-小凹程序与S100P/AKR1B10/CEACAM5/REG4/GPX2",
      "肠型或肿瘤候选程序并存；该cluster同时含正常胃和原发肿瘤来源细胞。",
      "组织来源和这些marker均不能单独证明恶性，需在F1.5细分并由F1.6 CNV联合判断。"
    ),
    "4" = paste0(
      "JCHAIN/TNFRSF17/MZB1/DERL3/IGKC及大量免疫球蛋白共同支持浆细胞；",
      "IGHM和IGHA并存，暂不按单一抗体亚型进一步细分。"
    ),
    "5" = paste0(
      "LUM/DCN/PDGFRA/DPT/MFAP4/TCF21为连贯成纤维细胞程序，",
      "CCL11/CXCL14/POSTN提示状态差异；缺少足够FAP/INHBA/CTHRC1证据，",
      "因此不直接命名为CAF亚型。"
    ),
    "6" = paste0(
      "VWF/CDH5/ENG/PLVAP/EMCN/CLDN5/TEK成套支持血管内皮细胞；",
      "即使DoubletFinder-only比例偏高，谱系marker仍连贯且跨样本存在。"
    ),
    "7" = paste0(
      "GNLY/NKG7/KLRD1/PRF1/GZMB强烈支持NK样细胞毒程序；",
      "同时CD3D/E及TRBC在部分细胞中表达，因此保留T/NK粗谱系，",
      "不把整个cluster写成纯NK。"
    ),
    "8" = paste0(
      "IGHG1/IGHG3/IGHG4与MZB1/DERL3/JCHAIN/TNFRSF17/XBP1共同支持",
      "IGHG-high抗体分泌浆细胞。"
    ),
    "9" = paste0(
      "EPCAM/KRT8/KRT18/KRT19/CDH1支持上皮，MUC5AC/TFF1/GKN1/GKN2",
      "支持胃表面-小凹样分化；CYP3A5/CLDN4升高作为状态记录，",
      "不能据此直接判断恶性。"
    ),
    "10" = paste0(
      "CD3D/CD3E/CD3G与CD8A/CD8B支持CD8 T细胞，IFNG/CCL4/TNF/",
      "GZMB/NKG7提示活化细胞毒状态。"
    ),
    "11" = paste0(
      "C1QA/C1QB/C1QC、CD163/MRC1/CSF1R/MS4A7支持巨噬细胞，",
      "FOLR2/TREM2/CCL18/APOC1/SPP1提示肿瘤相关巨噬状态，但不进一步硬分亚型。"
    ),
    "12" = paste0(
      "MS4A1/CD19/CD79A/CD79B/CD22/BANK1成套支持常规B细胞；",
      "缺少MZB1/DERL3主导的浆细胞程序。"
    ),
    "13" = paste0(
      "FCN1/S100A8/S100A9/CD14/LYZ与IL1B/TREM1共同支持炎症性单核细胞，",
      "区别于c11的成熟C1Q巨噬细胞程序。"
    ),
    "14" = paste0(
      "IGHA1/IGHA2与MZB1/DERL3/JCHAIN/TNFRSF17/IGKC共同支持",
      "IGHA-high抗体分泌浆细胞。"
    ),
    "15" = paste0(
      "RGS5/NOTCH3/PDGFRB/CSPG4与MYH11/ACTA2/TAGLN/CNN1共同支持",
      "周细胞-平滑肌样收缩程序；内皮VWF/PECAM1并不主导。"
    ),
    "16" = paste0(
      "EPCAM/KRT8/KRT18/KRT19支持上皮，MUC6/TFF2与PGC/LIPF/PGA3",
      "共同提示黏液颈细胞-主细胞样腺体分泌程序；两种分化程序并存，",
      "暂不强行拆成单一胃上皮亚型。"
    ),
    "17" = paste0(
      "CPA3/TPSAB1/TPSB2/MS4A2/KIT/HDC为成套肥大细胞marker，",
      "谱系清楚且跨样本存在。"
    ),
    "18" = paste0(
      "HSPA6/JUN/FOS/DUSP1/HSPA1A-B提示强应激或解离相关状态；",
      "MZB1/IGKC/DERL3/CD79A仍提供B/浆细胞谱系证据。",
      "应激是cell_state而不是细胞类型，也不能仅凭应激自动删除真实细胞。"
    ),
    "19" = paste0(
      "EPCAM/KRT8/KRT18/KRT19及PGA3/PGA4/PGA5/PGC/LIPF支持",
      "胃主细胞样上皮；MT1G/MT1H/MT1X/MT2A为金属硫蛋白和泛氧化应激应答，",
      "不是MT-CO/MT-ND等线粒体编码基因，不能据此推断线粒体转录本升高。"
    ),
    "20" = paste0(
      "CD1C/FCER1A/CD1E/CLEC10A与HLA-DP/DQ/DR共同支持树突细胞，",
      "LAMP3/IDO1/CD86提示抗原呈递和成熟/活化状态。"
    ),
    "21" = paste0(
      "EPCAM/KRT8/KRT18/KRT19支持上皮背景，CHGA/CHGB/NEUROD1/",
      "INSM1/PCSK1/SYP/GHRL为连贯神经内分泌程序，支持胃肠内分泌上皮。"
    ),
    "22" = paste0(
      "UPK3B/LRRN4/MSLN/WT1/PRG4/CALB2/ITLN1形成连贯间皮细胞程序，",
      "同时出现角蛋白和ECM基因符合间皮特征而非必然双细胞。",
      "DoubletFinder-only为16.3%，但该群跨28位患者且marker内部一致，",
      "因此保留为Mesothelial并维持中等置信，不自动删除。"
    ),
    "23" = paste0(
      "CD3D/CD3E/CD3G支持T细胞，CCR6/KLRB1/IL7R/RORA与CCL20高表达",
      "提示Th17-like CD4样状态；缺少更完整细胞因子证据，故只作like注释。"
    ),
    "24" = paste0(
      "MKI67/TOP2A/CDK1/UBE2C等细胞周期基因主导，CD3D/CD3E/CD3G及TRBC1/2",
      "仍支持T细胞母谱系；因此命名为Cycling T而不是独立细胞类型。"
    )
  )

  decision_names <- list(
    names(major), names(minor), names(state), names(confidence), names(reason)
  )
  if (!all(vapply(decision_names, identical, logical(1), expected_clusters))) {
    stop("F1注释草案的cluster键不一致。")
  }

  handling <- setNames(
    ifelse(
      major == "Epithelial",
      paste0(
        "保留为F1.5上皮候选；恶性身份等待F1.6；",
        "不得仅凭组织来源或亚型marker判定恶性。"
      ),
      paste0(
        "保留在全细胞对象；不进入F1.5上皮主对象；",
        "不得仅凭单一状态或DoubletFinder-only标签自动删除。"
      )
    ),
    names(major)
  )
  handling["3"] <- paste0(
    "保留为F1.5上皮候选；恶性身份等待F1.6；",
    "不得作为inferCNV或CopyKAT正常reference。"
  )
  handling["18"] <- paste0(
    "保留在全细胞对象；当前证据不纳入F1.5上皮主对象；",
    "不得仅因stress-high删除，研究者批准后冻结B/Plasma粗谱系。"
  )
  handling["19"] <- paste0(
    "保留为F1.5上皮候选；记录metallothionein-high状态；",
    "不得把MT1家族或MT2A解释为线粒体编码基因。"
  )
  handling["22"] <- paste0(
    "保留在全细胞对象并作为Mesothelial独立粗谱系；不进入F1.5上皮主对象；",
    "复核DoubletFinder-only富集，但不因其单独删群。"
  )

  decisions <- data.frame(
    cluster = expected_clusters,
    cell_type_major = unname(major),
    cell_type_minor = unname(minor),
    cell_state = unname(state),
    annotation_confidence = unname(confidence),
    annotation_reason = unname(reason),
    annotation_review_status =
      "complete_draft_pending_researcher_approval",
    downstream_handling_before_full_approval = unname(handling),
    stringsAsFactors = FALSE
  )

  matched <- match(decisions$cluster, as.character(template$cluster))
  decision_fields <- setdiff(colnames(decisions), "cluster")
  for (field in decision_fields) {
    template[[field]][matched] <- decisions[[field]]
  }
  template
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
    min_cells = config$sct$minimum_cells_per_gene,
    conserve.memory = config$sct$conserve_memory,
    return.only.var.genes = TRUE,
    seed.use = config$seed,
    verbose = TRUE
  )

  SeuratObject::DefaultAssay(object) <- "SCT"
  # 多样本SCT对象的3000个候选HVG中，只有各样本共同保留在scale.data里的基因
  # 能进入PCA。显式取交集，避免Seurat静默丢弃未缩放基因。
  pca_features <- intersect(
    SeuratObject::VariableFeatures(object[["SCT"]]),
    rownames(SeuratObject::LayerData(object, assay = "SCT", layer = "scale.data"))
  )
  if (length(pca_features) < 500L) {
    stop(stage_label, "可用于PCA的SCT已缩放高变基因少于500个，需检查逐样本SCT结果。")
  }
  object <- Seurat::RunPCA(
    object,
    assay = "SCT",
    features = pca_features,
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
  # Leiden本身已显式固定random.seed；切回顺序future，避免其内部并行封装
  # 对已固定随机数产生误报。耗时主体SCT和UMAP此前仍按服务器并行。
  future::plan(future::sequential)
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
      "SCTransform_vars_to_regress", "SCTransform_min_cells_per_gene",
      "PCA_npcs", "Harmony_group_by",
      "main_dims", "default_resolution", "UCell_F2_input",
      "SCTransform_and_UMAP_future_workers", "future_globals_max_GB",
      "scDblFinder_workers",
      "DecontX_timing", "DecontX_cluster_label_source",
      "DecontX_minimum_reliable_lineages", "DecontX_background",
      "inferCNV_analysis_mode", "inferCNV_internal_subclustering",
      "inferCNV_observation_group", "inferCNV_min_cells_per_gene",
      "inferCNV_window_length", "inferCNV_HMM", "inferCNV_HMM_type",
      "inferCNV_HMM_report_by", "inferCNV_BayesMaxPNormal",
      "inferCNV_reassignCNVs", "inferCNV_output_format",
      "inferCNV_write_expr_matrix", "inferCNV_scipen",
      "inferCNV_bitmap_type", "CNV_sample_workers",
      "inferCNV_threads",
      "CopyKAT_cores", "CopyKAT_internal_min_genes_per_cell",
      "CopyKAT_input_scope", "CopyKAT_arm_A",
      "CopyKAT_arm_B", "CopyKAT_arm_C",
      "CopyKAT_holdout_seed", "CopyKAT_baseline_suspect_rule",
      "CopyKAT_package_data_loading"
    ),
    value = c(
      config$seed, config$sct$vst_flavor, config$sct$variable_features_n,
      "NULL", config$sct$minimum_cells_per_gene,
      config$sct$pca_npcs, config$sct$harmony_group,
      paste(range(config$sct$main_dims), collapse = ":"), config$sct$default_resolution,
      "RNA_raw_counts",
      config$execution$future_workers, config$execution$future_globals_max_gb,
      config$doublet$scdblfinder_workers,
      "after_researcher_approved_cluster_annotation",
      "researcher_approved_seurat_cluster",
      config$ambient$minimum_reliable_lineages, "NULL",
      config$cnv$infercnv_analysis_mode,
      identical(config$cnv$infercnv_analysis_mode, "subclusters"),
      config$cnv$infercnv_observation_group,
      config$cnv$infercnv_min_cells_per_gene,
      config$cnv$infercnv_window_length,
      config$cnv$infercnv_hmm, config$cnv$infercnv_hmm_type,
      config$cnv$infercnv_hmm_report_by,
      config$cnv$infercnv_bayes_max_p_normal,
      config$cnv$infercnv_reassign_cnvs,
      config$cnv$infercnv_output_format,
      config$cnv$infercnv_write_expr_matrix,
      config$cnv$infercnv_scipen,
      config$cnv$infercnv_bitmap_type,
      config$cnv$sample_workers,
      config$cnv$infercnv_threads, config$cnv$copykat_cores,
      config$cnv$copykat_internal_min_genes_per_cell,
      "all_QC_singlet_cells_from_current_sample_only",
      "primary_self_estimated_baseline_without_norm.cell.names",
      paste0("sensitivity_same_sample_immune_if_at_least_", config$cnv$minimum_reference_cells),
      "Normal_Gastric_only_two_fold_held_out_epithelium_plus_same_sample_immune",
      config$cnv$copykat_holdout_seed,
      paste0(
        "self_estimated_diploid_candidate_epithelial_fraction_gt_",
        config$cnv$copykat_baseline_suspect_epithelial_fraction,
        "_is_not_evaluable_baseline_suspect"
      ),
      "library(copykat)_LazyData_full.anno_DNA.hg20_cyclegenes"
    ),
    stringsAsFactors = FALSE
  )
  out <- rbind(out, base_parameters)
  if (!is.null(extra)) out <- rbind(out, extra)
  f1_write_tsv(out, file.path(config$paths$annotation_dir, "F1_parameters_and_versions.tsv"))
}
