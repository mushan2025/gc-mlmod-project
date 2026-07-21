# F1.6 inferCNV、CopyKAT与恶性上皮判定 ------------------------------------
#
# 生物学目的：在可靠上皮细胞内，用大片段CNV模式、CopyKAT非整倍体结果、
# 上皮marker和组织背景联合区分恶性、非恶性与不确定细胞。
#
# 重要边界：inferCNV和CopyKAT都来自RNA表达，只能称互补方法稳健性，不能称
# 独立DNA验证。正式CNV计算需要额外参数 --approve-cnv-execution，防止误触。

file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), winslash = "/"))
source(file.path(script_dir, "F1_config.R"), encoding = "UTF-8")
source(file.path(script_dir, "F1_utils.R"), encoding = "UTF-8")

args <- f1_parse_args()
config <- f1_build_config(args$project_root)
stage <- "F1.6"
outputs <- c(
  config$paths$malignancy_review_template,
  config$paths$object_05,
  config$paths$object_06a,
  config$paths$object_06b,
  file.path(config$paths$malignancy_dir, "malignant_cell_calling_summary.tsv"),
  file.path(config$paths$malignancy_dir, "F1_final_report.md")
)

if (!args$execute) {
  f1_stage_dry_run(
    stage,
    "按样本运行raw-count inferCNV/CopyKAT，经人工审核后形成06a/06b恶性上皮对象",
    c(
      config$paths$object_04, config$paths$object_03,
      config$paths$epithelial_review_approved, config$paths$gene_order
    ),
    outputs,
    config$packages[[stage]]
  )
  cat("CNV安全开关：正式计算还需 --approve-cnv-execution。\n")
  quit(save = "no", status = 0)
}

f1_require_packages(config$packages[[stage]], stage)
f1_check_f0_ready(config)
if (!file.exists(config$paths$object_04)) stop("请先完成F1.5：", config$paths$object_04)
if (!file.exists(config$paths$object_03)) stop("缺少全细胞注释对象：", config$paths$object_03)
if (!file.exists(config$paths$epithelial_review_approved)) {
  stop("F1.6前必须审核上皮cluster并建立：", config$paths$epithelial_review_approved)
}
if (!file.exists(config$paths$gene_order)) stop("缺少inferCNV gene order：", config$paths$gene_order)
f1_prepare_directories(config)
set.seed(config$seed)
f1_append_log(config, stage, "开始准备上皮CNV证据与联合恶性判定")

epithelial <- readRDS(config$paths$object_04)
all_cells <- readRDS(config$paths$object_03)
f1_require_columns(
  epithelial[[]],
  c("sample_id", "patient_id", "group_analysis", "epithelial_cluster_id"),
  "F1.5 epithelial object metadata"
)
f1_require_columns(
  all_cells[[]],
  c("sample_id", "cell_type_major", "annotation_confidence"),
  "F1.4 all-cell object metadata"
)
epithelial_review <- f1_read_tsv(config$paths$epithelial_review_approved)
f1_require_columns(
  epithelial_review,
  c("epithelial_cluster", "include_in_malignancy", "contamination_status", "epithelial_subtype", "review_reason"),
  "F1_epithelial_cluster_review_approved.tsv"
)
clusters <- sort(unique(as.character(epithelial$epithelial_cluster_id)))
if (anyDuplicated(as.character(epithelial_review$epithelial_cluster)) ||
    !setequal(as.character(epithelial_review$epithelial_cluster), clusters)) {
  stop("批准的上皮cluster复核表必须与F1.5当前cluster一对一对应。")
}
epithelial_review$include_in_malignancy <- f1_as_logical(
  epithelial_review$include_in_malignancy,
  "include_in_malignancy"
)
approved_clusters <- as.character(epithelial_review$epithelial_cluster[epithelial_review$include_in_malignancy])
if (!length(approved_clusters)) stop("没有上皮cluster获准进入恶性判定。")
epithelial_review_index <- match(epithelial$epithelial_cluster_id, epithelial_review$epithelial_cluster)
epithelial$epithelial_subtype <- epithelial_review$epithelial_subtype[epithelial_review_index]
epithelial$epithelial_contamination_status <- epithelial_review$contamination_status[epithelial_review_index]
epithelial$epithelial_review_reason <- epithelial_review$review_reason[epithelial_review_index]
candidate_cells <- colnames(epithelial)[epithelial$epithelial_cluster_id %in% approved_clusters]

resource_summary <- data.frame(
  item = c("candidate_epithelial_cells", "candidate_samples", "candidate_clusters", "requested_infercnv_threads", "requested_copykat_cores"),
  value = c(
    length(candidate_cells), length(unique(epithelial$sample_id[match(candidate_cells, colnames(epithelial))])),
    length(approved_clusters), config$cnv$infercnv_threads, config$cnv$copykat_cores
  ),
  interpretation = c(
    "实际CNV观察细胞数", "逐样本运行单元", "人工批准进入CNV的上皮cluster数",
    "inferCNV并行线程", "CopyKAT并行核数"
  ),
  stringsAsFactors = FALSE
)
f1_write_tsv(resource_summary, file.path(config$paths$malignancy_dir, "F1_CNV_actual_input_resource_summary.tsv"))

cnv_dir <- file.path(config$paths$malignancy_dir, "cnv_by_sample")
dir.create(cnv_dir, recursive = TRUE, showWarnings = FALSE)
existing_scores <- file.exists(file.path(cnv_dir, "infercnv_cell_scores_all.tsv"))
existing_copykat <- file.exists(file.path(cnv_dir, "copykat_cell_calls_all.tsv"))
if (config$cnv$require_explicit_execution_approval && !(existing_scores && existing_copykat) && !args$approve_cnv_execution) {
  stop(
    "已输出F1.6实际输入规模，但尚未运行高资源CNV步骤。批准设备后重新运行并增加：",
    " --approve-cnv-execution"
  )
}

all_cells <- f1_join_assay(all_cells, "RNA")
full_counts <- SeuratObject::LayerData(all_cells, assay = "RNA", layer = "counts")
f1_assert_integer_counts(full_counts, "F1.6 full-cell RNA counts")
gene_order <- f1_read_tsv(config$paths$gene_order)
f1_require_columns(gene_order, c("gene_symbol", "chromosome", "start", "end"), "inferCNV gene order")
gene_order <- gene_order[grepl("^chr([1-9]|1[0-9]|2[0-2])$", gene_order$chromosome), , drop = FALSE]
gene_order_duplicate_count <- sum(duplicated(gene_order$gene_symbol))
gene_order <- gene_order[!duplicated(gene_order$gene_symbol), , drop = FALSE]
matched_gene_order <- gene_order[gene_order$gene_symbol %in% rownames(full_counts), , drop = FALSE]
if (!nrow(matched_gene_order)) stop("inferCNV gene order与RNA gene symbol没有交集。")
mapping_report <- data.frame(
  gene_order_source = config$paths$gene_order,
  source_chr1_22_genes = nrow(gene_order),
  RNA_features = nrow(full_counts),
  matched_genes = nrow(matched_gene_order),
  matched_fraction_of_gene_order = nrow(matched_gene_order) / nrow(gene_order),
  duplicate_gene_symbols_removed = gene_order_duplicate_count,
  chromosomes = paste(unique(matched_gene_order$chromosome), collapse = "|"),
  input_assay_layer = "RNA_counts",
  status = if (nrow(matched_gene_order) >= 1000) "PASS_BASIC_ID_COMPATIBILITY" else "FAIL_LIKELY_ID_MISMATCH",
  stringsAsFactors = FALSE
)
f1_write_tsv(mapping_report, file.path(config$paths$malignancy_dir, "gene_order_mapping_report.tsv"))
if (mapping_report$status == "FAIL_LIKELY_ID_MISMATCH") {
  stop("inferCNV gene-order匹配少于1000个基因，提示gene ID或输入明显不匹配。")
}

choose_reference_cells <- function(sample_id, all_cells, minimum) {
  meta <- all_cells[[]]
  confidence_ok <- meta$annotation_confidence %in% c("high", "medium")
  same_sample <- meta$sample_id == sample_id
  tnk <- rownames(meta)[same_sample & confidence_ok & meta$cell_type_major == "T/NK"]
  bplasma <- rownames(meta)[same_sample & confidence_ok & meta$cell_type_major == "B/Plasma"]
  selected <- tnk
  source <- "same_sample_T_NK"
  if (length(selected) < minimum) {
    selected <- unique(c(selected, bplasma))
    source <- "same_sample_T_NK_plus_B_Plasma"
  }
  if (length(selected) < minimum) {
    pooled <- rownames(meta)[confidence_ok & meta$cell_type_major %in% c("T/NK", "B/Plasma")]
    selected <- unique(c(selected, sort(pooled)))
    source <- "pooled_high_medium_confidence_T_NK_B_Plasma"
  }
  selected <- utils::head(selected, max(minimum, 500L))
  list(cells = selected, source = source)
}

infercnv_cell_scores <- function(infer_object, reference_cells, sample_id, epithelial_cluster) {
  expr <- infer_object@expr.data
  reference_cells <- intersect(reference_cells, colnames(expr))
  if (!length(reference_cells)) stop(sample_id, " inferCNV结果中没有reference cell。")
  ref_center <- Matrix::rowMeans(expr[, reference_cells, drop = FALSE])
  scores <- numeric(ncol(expr))
  chunk_starts <- seq.int(1L, ncol(expr), by = 250L)
  for (start in chunk_starts) {
    index <- start:min(start + 249L, ncol(expr))
    block <- as.matrix(expr[, index, drop = FALSE])
    scores[index] <- colMeans(abs(sweep(block, 1, ref_center, FUN = "-")))
  }
  names(scores) <- colnames(expr)
  threshold <- as.numeric(stats::quantile(scores[reference_cells], 0.95, na.rm = TRUE, type = 7))
  cells <- intersect(names(epithelial_cluster), names(scores))
  data.frame(
    cell_id_final = cells,
    sample_id = sample_id,
    epithelial_cluster = unname(epithelial_cluster[cells]),
    infercnv_cell_burden = unname(scores[cells]),
    reference_background_P95 = threshold,
    above_reference_P95 = unname(scores[cells]) > threshold,
    stringsAsFactors = FALSE
  )
}

extract_copykat_calls <- function(result, candidate_cells, sample_id) {
  prediction <- result$prediction
  if (is.null(prediction)) stop(sample_id, " CopyKAT结果缺少prediction。")
  prediction <- as.data.frame(prediction, stringsAsFactors = FALSE)
  f1_require_columns(prediction, c("cell.names", "copykat.pred"), paste0(sample_id, " CopyKAT prediction"))
  raw_call <- setNames(tolower(as.character(prediction$copykat.pred)), as.character(prediction$cell.names))
  call <- setNames(rep("uncalled", length(candidate_cells)), candidate_cells)
  matched <- intersect(candidate_cells, names(raw_call))
  normalized <- ifelse(
    grepl("aneuploid", raw_call[matched]), "aneuploid",
    ifelse(grepl("diploid", raw_call[matched]), "diploid", "uncalled")
  )
  call[matched] <- normalized
  data.frame(
    cell_id_final = candidate_cells,
    sample_id = sample_id,
    copykat_call = unname(call[candidate_cells]),
    copykat_returned = candidate_cells %in% names(raw_call),
    copykat_note = ifelse(candidate_cells %in% names(raw_call), "returned_by_copykat", "filtered_or_uncalled_by_copykat"),
    stringsAsFactors = FALSE
  )
}

if (existing_scores && existing_copykat) {
  infer_scores_all <- f1_read_tsv(file.path(cnv_dir, "infercnv_cell_scores_all.tsv"))
  copykat_calls_all <- f1_read_tsv(file.path(cnv_dir, "copykat_cell_calls_all.tsv"))
  reference_summary_all <- f1_read_tsv(file.path(cnv_dir, "infercnv_reference_summary.tsv"))
  f1_append_log(config, stage, "复用已保存的逐样本inferCNV与CopyKAT结果")
} else {
  infer_score_rows <- list()
  copykat_rows <- list()
  reference_rows <- list()
  sample_ids <- sort(unique(as.character(epithelial$sample_id[match(candidate_cells, colnames(epithelial))])))

  for (sample_id in sample_ids) {
    message("[", stage, "] CNV分析：", sample_id)
    obs_cells <- candidate_cells[epithelial$sample_id[match(candidate_cells, colnames(epithelial))] == sample_id]
    obs_cluster <- setNames(
      as.character(epithelial$epithelial_cluster_id[match(obs_cells, colnames(epithelial))]),
      obs_cells
    )
    reference <- choose_reference_cells(sample_id, all_cells, config$cnv$minimum_reference_cells)
    reference_cells <- setdiff(intersect(reference$cells, colnames(full_counts)), obs_cells)
    if (length(reference_cells) < config$cnv$minimum_reference_cells) {
      stop(sample_id, "可用高/中置信免疫reference不足", config$cnv$minimum_reference_cells, "个。")
    }
    sample_dir <- file.path(cnv_dir, sample_id)
    infer_dir <- file.path(sample_dir, "infercnv")
    copykat_dir <- file.path(sample_dir, "copykat")
    dir.create(infer_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(copykat_dir, recursive = TRUE, showWarnings = FALSE)

    input_cells <- c(obs_cells, reference_cells)
    infer_counts <- full_counts[matched_gene_order$gene_symbol, input_cells, drop = FALSE]
    gene_order_file <- file.path(sample_dir, paste0(sample_id, "_infercnv_gene_order.tsv"))
    data.table::fwrite(
      matched_gene_order[, c("gene_symbol", "chromosome", "start", "end")],
      gene_order_file,
      sep = "\t",
      quote = FALSE,
      col.names = FALSE
    )
    annotation_file <- file.path(sample_dir, paste0(sample_id, "_infercnv_annotations.tsv"))
    annotation <- data.frame(
      cell = input_cells,
      group = c(paste0("obs__", unname(obs_cluster[obs_cells])), rep("reference", length(reference_cells))),
      stringsAsFactors = FALSE
    )
    data.table::fwrite(annotation, annotation_file, sep = "\t", quote = FALSE, col.names = FALSE)

    infer_object <- infercnv::CreateInfercnvObject(
      raw_counts_matrix = infer_counts,
      annotations_file = annotation_file,
      delim = "\t",
      gene_order_file = gene_order_file,
      ref_group_names = "reference"
    )
    set.seed(config$seed)
    infer_result <- infercnv::run(
      infer_object,
      cutoff = config$cnv$infercnv_cutoff,
      out_dir = infer_dir,
      cluster_by_groups = TRUE,
      cluster_references = TRUE,
      denoise = config$cnv$infercnv_denoise,
      HMM = config$cnv$infercnv_hmm,
      num_threads = config$cnv$infercnv_threads,
      no_plot = FALSE,
      save_rds = TRUE,
      save_final_rds = TRUE,
      resume_mode = TRUE
    )
    saveRDS(infer_result, file.path(sample_dir, paste0(sample_id, "_infercnv_final_object.rds")), compress = FALSE)
    infer_score_rows[[sample_id]] <- infercnv_cell_scores(infer_result, reference_cells, sample_id, obs_cluster)

    copykat_input <- full_counts[, input_cells, drop = FALSE]
    old_wd <- getwd()
    setwd(copykat_dir)
    copykat_result <- tryCatch({
      set.seed(config$seed)
      copykat::copykat(
        rawmat = copykat_input,
        id.type = "S",
        cell.line = "no",
        ngene.chr = config$cnv$copykat_ngene_chr,
        min.gene.per.cell = 200,
        LOW.DR = 0.05,
        UP.DR = 0.10,
        win.size = config$cnv$copykat_win_size,
        norm.cell.names = reference_cells,
        KS.cut = config$cnv$copykat_ks_cut,
        sam.name = sample_id,
        distance = "euclidean",
        output.seg = "FALSE",
        plot.genes = "TRUE",
        genome = "hg20",
        n.cores = config$cnv$copykat_cores
      )
    }, finally = {
      setwd(old_wd)
    })
    saveRDS(copykat_result, file.path(sample_dir, paste0(sample_id, "_copykat_result.rds")), compress = FALSE)
    copykat_rows[[sample_id]] <- extract_copykat_calls(copykat_result, obs_cells, sample_id)
    reference_rows[[sample_id]] <- data.frame(
      sample_id = sample_id,
      observation_cells = length(obs_cells),
      reference_cells = length(reference_cells),
      reference_source = reference$source,
      reference_lineages = paste(sort(unique(all_cells$cell_type_major[match(reference_cells, colnames(all_cells))])), collapse = "|"),
      infercnv_cutoff = config$cnv$infercnv_cutoff,
      infercnv_HMM = config$cnv$infercnv_hmm,
      infercnv_denoise = config$cnv$infercnv_denoise,
      stringsAsFactors = FALSE
    )
    rm(infer_counts, infer_object, infer_result, copykat_input, copykat_result)
    gc(verbose = FALSE)
  }
  infer_scores_all <- do.call(rbind, infer_score_rows)
  copykat_calls_all <- do.call(rbind, copykat_rows)
  reference_summary_all <- do.call(rbind, reference_rows)
  f1_write_tsv(infer_scores_all, file.path(cnv_dir, "infercnv_cell_scores_all.tsv"))
  f1_write_tsv(copykat_calls_all, file.path(cnv_dir, "copykat_cell_calls_all.tsv"))
  f1_write_tsv(reference_summary_all, file.path(cnv_dir, "infercnv_reference_summary.tsv"))
}

f1_require_columns(
  infer_scores_all,
  c("cell_id_final", "sample_id", "epithelial_cluster", "infercnv_cell_burden", "reference_background_P95", "above_reference_P95"),
  "infercnv_cell_scores_all.tsv"
)
f1_require_columns(copykat_calls_all, c("cell_id_final", "sample_id", "copykat_call"), "copykat_cell_calls_all.tsv")
if (!setequal(candidate_cells, infer_scores_all$cell_id_final)) {
  stop("inferCNV cell score没有完整覆盖获准的上皮候选细胞。")
}
if (!setequal(candidate_cells, copykat_calls_all$cell_id_final)) {
  stop("CopyKAT call表没有完整覆盖获准的上皮候选细胞。")
}

infer_cluster_summary <- do.call(rbind, lapply(
  split(infer_scores_all, interaction(infer_scores_all$sample_id, infer_scores_all$epithelial_cluster, drop = TRUE)),
  function(x) data.frame(
    sample_id = x$sample_id[[1]],
    epithelial_cluster = x$epithelial_cluster[[1]],
    cell_count = nrow(x),
    infercnv_burden_median = stats::median(x$infercnv_cell_burden, na.rm = TRUE),
    reference_background_P95 = unique(x$reference_background_P95)[[1]],
    fraction_above_reference_P95 = mean(f1_as_logical(x$above_reference_P95, "above_reference_P95"), na.rm = TRUE),
    broad_segment_support = "pending_researcher_heatmap_review",
    stringsAsFactors = FALSE
  )
))
f1_write_tsv(infer_cluster_summary, file.path(config$paths$malignancy_dir, "infercnv_cluster_summary.tsv"))
f1_write_tsv(copykat_calls_all, file.path(config$paths$malignancy_dir, "copykat_cell_calls.tsv"))

copykat_cluster <- do.call(rbind, lapply(
  split(copykat_calls_all, interaction(copykat_calls_all$sample_id, epithelial$epithelial_cluster_id[match(copykat_calls_all$cell_id_final, colnames(epithelial))], drop = TRUE)),
  function(x) {
    cluster_id <- epithelial$epithelial_cluster_id[match(x$cell_id_final[[1]], colnames(epithelial))]
    data.frame(
      sample_id = x$sample_id[[1]],
      epithelial_cluster = cluster_id,
      copykat_aneuploid = sum(x$copykat_call == "aneuploid"),
      copykat_diploid = sum(x$copykat_call == "diploid"),
      copykat_uncalled = sum(x$copykat_call == "uncalled"),
      stringsAsFactors = FALSE
    )
  }
))
review_template <- merge(infer_cluster_summary, copykat_cluster, by = c("sample_id", "epithelial_cluster"), all = TRUE, sort = TRUE)
review_template$infercnv_broad_segment_support <- ""
review_template$tumor_program_support <- ""
review_template$normal_program_support <- ""
review_template$review_note <- ""
f1_write_tsv(review_template, config$paths$malignancy_review_template)

if (!file.exists(config$paths$malignancy_review_approved)) {
  f1_append_log(config, stage, paste0("CNV结果已生成，等待热图和marker联合审核：", config$paths$malignancy_review_template))
  stop(
    "F1.6已生成inferCNV/CopyKAT摘要和恶性判定模板。请审核热图及marker后建立：",
    config$paths$malignancy_review_approved
  )
}

approved <- f1_read_tsv(config$paths$malignancy_review_approved)
review_fields <- c(
  "sample_id", "epithelial_cluster", "infercnv_broad_segment_support",
  "tumor_program_support", "normal_program_support", "review_note"
)
f1_require_columns(approved, review_fields, "F1_malignancy_cluster_review_approved.tsv")
approved$key <- paste(approved$sample_id, approved$epithelial_cluster, sep = "__")
template_keys <- paste(review_template$sample_id, review_template$epithelial_cluster, sep = "__")
if (anyDuplicated(approved$key) || !setequal(approved$key, template_keys)) {
  stop("恶性审核批准表必须与当前sample_id × epithelial_cluster一对一对应。")
}
allowed_infer <- c("strong", "weak", "absent", "not_evaluable")
if (any(!approved$infercnv_broad_segment_support %in% allowed_infer)) {
  stop("infercnv_broad_segment_support只允许strong、weak、absent或not_evaluable。")
}
approved$tumor_program_support <- f1_as_logical(approved$tumor_program_support, "tumor_program_support")
approved$normal_program_support <- f1_as_logical(approved$normal_program_support, "normal_program_support")

epi_key <- paste(epithelial$sample_id, epithelial$epithelial_cluster_id, sep = "__")
review_index <- match(epi_key, approved$key)
copy_index <- match(colnames(epithelial), copykat_calls_all$cell_id_final)
epithelial$infercnv_broad_segment_support <- approved$infercnv_broad_segment_support[review_index]
epithelial$tumor_program_support <- approved$tumor_program_support[review_index]
epithelial$normal_program_support <- approved$normal_program_support[review_index]
epithelial$copykat_call <- copykat_calls_all$copykat_call[copy_index]
epithelial$infercnv_cell_burden <- infer_scores_all$infercnv_cell_burden[match(colnames(epithelial), infer_scores_all$cell_id_final)]

cluster_approved <- epithelial$epithelial_cluster_id %in% approved_clusters
infer_strong <- epithelial$infercnv_broad_segment_support == "strong"
infer_negative_or_weak <- epithelial$infercnv_broad_segment_support %in% c("weak", "absent")
aneuploid <- epithelial$copykat_call == "aneuploid"
diploid <- epithelial$copykat_call == "diploid"
uncalled <- is.na(epithelial$copykat_call) | epithelial$copykat_call == "uncalled"
tumor_source <- epithelial$group_analysis %in% c("Primary_Tumor", "Peritoneal_Metastasis")
tumor_context <- tumor_source & epithelial$tumor_program_support
normal_context <- epithelial$normal_program_support

label <- rep("epithelial_uncertain", ncol(epithelial))
label[!cluster_approved] <- "exclude_contamination_or_doublet"
label[cluster_approved & infer_negative_or_weak & diploid & normal_context] <- "non_malignant_epithelial"
label[cluster_approved & aneuploid & !infer_strong & tumor_context] <- "malignant_probable_copykat"
label[cluster_approved & infer_strong & (diploid | uncalled) & tumor_context] <- "malignant_probable_infercnv"
label[cluster_approved & infer_strong & aneuploid] <- "malignant_high_confidence"
epithelial$malignancy_label <- label
epithelial$include_in_06a <- label %in% c(
  "malignant_high_confidence", "malignant_probable_infercnv", "malignant_probable_copykat"
)
epithelial$include_in_06b <- label == "malignant_high_confidence"

calling <- data.frame(
  cell_id_final = colnames(epithelial),
  sample_id = epithelial$sample_id,
  patient_id = epithelial$patient_id,
  group_analysis = epithelial$group_analysis,
  epithelial_cluster = epithelial$epithelial_cluster_id,
  infercnv_cell_burden = epithelial$infercnv_cell_burden,
  infercnv_broad_segment_support = epithelial$infercnv_broad_segment_support,
  copykat_call = epithelial$copykat_call,
  tumor_program_support = epithelial$tumor_program_support,
  normal_program_support = epithelial$normal_program_support,
  malignancy_label = epithelial$malignancy_label,
  include_in_06a = epithelial$include_in_06a,
  include_in_06b = epithelial$include_in_06b,
  stringsAsFactors = FALSE
)
f1_write_tsv(calling, file.path(config$paths$malignancy_dir, "malignant_cell_calling_summary.tsv"))
f1_write_tsv(epithelial_review, file.path(config$paths$annotation_dir, "epithelial_contamination_review.tsv"))
f1_save_rds_atomic(epithelial, config$paths$object_05, compress = FALSE)

n_main <- sum(epithelial$include_in_06a)
n_high <- sum(epithelial$include_in_06b)
if (n_main == 0L || n_high == 0L) {
  stop(
    "联合判定后06a或06b没有细胞（06a=", n_main, ", 06b=", n_high,
    "）。05对象和摘要已保存，请检查reference与审核证据，不得放宽规则凑数。"
  )
}
main_object <- subset(epithelial, cells = colnames(epithelial)[epithelial$include_in_06a])
high_object <- subset(epithelial, cells = colnames(epithelial)[epithelial$include_in_06b])
f1_save_rds_atomic(main_object, config$paths$object_06a, compress = FALSE)
f1_save_rds_atomic(high_object, config$paths$object_06b, compress = FALSE)

# 汇总全部公开细胞的关键metadata，方便后续按原始细胞追溯排除原因。
initial <- readRDS(config$paths$object_01)
final_meta <- initial[[]]
full_meta <- all_cells[[]]
for (field in c("seurat_clusters", "cell_type_major", "cell_type_minor", "cell_state", "annotation_confidence")) {
  final_meta[[field]] <- full_meta[[field]][match(rownames(final_meta), rownames(full_meta))]
}
for (field in c("epithelial_cluster_id", "copykat_call", "infercnv_cell_burden", "malignancy_label", "include_in_06a", "include_in_06b")) {
  final_meta[[field]] <- epithelial[[field, drop = TRUE]][match(rownames(final_meta), colnames(epithelial))]
}
final_meta$cell_id_final <- rownames(final_meta)
f1_write_tsv(final_meta, file.path(config$paths$malignancy_dir, "cell_metadata_final.tsv"))

label_counts <- table(epithelial$malignancy_label)
status <- "PASS_WITH_NOTED_LIMITATIONS"
report <- c(
  "# F1 Final Report", "",
  paste0("Generated at: ", f1_now()), "",
  "## Main Result", "",
  paste0("- F1 status: ", status),
  paste0("- QC/doublet-filtered all cells: ", ncol(all_cells)),
  paste0("- Reclustered epithelial cells: ", ncol(epithelial)),
  paste0("- 06a malignant main cells: ", n_main),
  paste0("- 06b high-confidence malignant cells: ", n_high),
  paste0("- Malignancy labels: ", paste0(names(label_counts), "=", as.integer(label_counts), collapse = "; ")), "",
  "## Methods", "",
  "- Fixed per-sample min.cells=3 feature space and five frozen cell thresholds were used.",
  "- scDblFinder was the deletion rule; DoubletFinder was sensitivity-only.",
  "- Per-sample SCTransform v2 used vars.to.regress=NULL; Harmony corrected sample_id only.",
  "- UCell/MLMOD and prognosis were not used in F1 decisions.",
  "- inferCNV and CopyKAT used RNA raw counts and were interpreted jointly with approved biological review.", "",
  "## Limitations", "",
  "- Raw/empty droplets are unavailable; cell calling, emptyDrops, SoupX and CellBender were not evaluable.",
  "- inferCNV and CopyKAT are both RNA-derived and do not constitute independent DNA validation.",
  "- CopyKAT diploid calls cannot exclude near-diploid malignant cells.",
  "- The public matrix and repeated mt_percent<=20 rule limit conclusions about cells above that boundary.", "",
  "## F2 Entry", "",
  "- 06a is the F2 main object; 06b is the high-confidence sensitivity object.",
  "- Enter F2 only after independent script/result review and user approval."
)
writeLines(report, file.path(config$paths$malignancy_dir, "F1_final_report.md"), useBytes = TRUE)
f1_write_parameter_versions(config)
f1_save_session_info(config, "F1_06_malignancy_inference")
f1_append_log(config, stage, sprintf("完成F1：06a=%d，06b=%d；状态=%s", n_main, n_high, status))
message("F1.6完成：", config$paths$object_06a)
