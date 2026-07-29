# F1 静态验证 ---------------------------------------------------------------
#
# 只解析脚本和核对文字契约，不读取真实表达矩阵、不加载缺失R包、不运行分析。

file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), winslash = "/"))
root <- normalizePath(file.path(script_dir, "..", ".."), winslash = "/", mustWork = TRUE)

scripts <- c(
  "F1_config.R", "F1_utils.R", "F1_01_import.R", "F1_02_qc_doublet.R",
  "F1_03_sct_harmony_cluster.R", "F1_04_annotation.R",
  "F1_05_epithelial_recluster.R", "F1_06_malignancy_inference.R", "run_F1.R"
)
for (script in scripts) {
  path <- file.path(script_dir, script)
  if (!file.exists(path)) stop("缺少脚本：", path)
  parse(file = path, encoding = "UTF-8")
  cat("PARSE PASS：", script, "\n")
}

read_all <- function(path) paste(readLines(path, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
combined <- paste(vapply(file.path(script_dir, scripts), read_all, character(1)), collapse = "\n")
required_patterns <- c(
  "vst.flavor = config\\$sct\\$vst_flavor",
  "vars.to.regress = config\\$sct\\$vars_to_regress",
  "group.by.vars = config\\$sct\\$harmony_group",
  "RNA_raw_counts",
  "nfeature_min_inclusive = 500L",
  "nfeature_max_exclusive = 6000L",
  "ncount_min_exclusive = 1000L",
  "mt_max_inclusive = 20",
  "hb_max_exclusive = 5",
  "no_ncount_sensitivity_enabled = TRUE",
  "fixed_qc_pass_no_ncount_sensitivity",
  "no_ncount_sensitivity_extra_vs_main",
  "intersect\\(c\\(\"BCmetric\", \"BCmvn\"\\)",
  "minimum_reliable_lineages = 2L",
  "z = coarse_labels",
  "completed_with_researcher_approved_coarse_labels",
  "infercnv_analysis_mode = \"subclusters\"",
  "infercnv_tumor_subcluster_partition_method = \"leiden\"",
  "infercnv_k_nn = 20L",
  "infercnv_leiden_resolution = \"auto\"",
  "analysis_mode = config\\$cnv\\$infercnv_analysis_mode",
  "tumor_subcluster_partition_method",
  "inspect_subclusters = config\\$cnv\\$infercnv_inspect_subclusters",
  "extract_infercnv_subcluster_membership",
  "infercnv_subcluster_membership.tsv",
  "final_plot_cell_order",
  "contains_final_expr_data_for_publication_redraw",
  "not_evaluable_insufficient_reference",
  "same_patient_normal_gastric",
  "balanced_other_patient_normal_gastric",
  "choose_copykat_known_normals",
  "copykat_external_reference = \"prohibited\"",
  "infercnv_decontX_corrected",
  "normal_context & !epithelial\\$tumor_program_support",
  "tumor_context & !normal_context",
  "aneuploid & !infer_strong",
  "candidate_lineage_marker_detection_pct",
  "same_lineage_normal_reference_cells",
  "same_lineage_normal_reference_samples",
  "intestinal_like_is_metaplasia_or_tumor_ambiguous",
  "--approve-cnv-execution"
)
for (pattern in required_patterns) {
  if (!grepl(pattern, combined, perl = TRUE)) stop("F1脚本缺少冻结契约：", pattern)
}

f12_text <- read_all(file.path(script_dir, "F1_02_qc_doublet.R"))
f14_text <- read_all(file.path(script_dir, "F1_04_annotation.R"))
f16_text <- read_all(file.path(script_dir, "F1_06_malignancy_inference.R"))
if (grepl("celda::decontX", f12_text, fixed = TRUE)) {
  stop("DecontX仍在F1.2运行；应在粗谱系注释批准后的F1.4运行。")
}
if (!grepl("celda::decontX", f14_text, fixed = TRUE)) {
  stop("F1.4缺少注释后逐样本DecontX。")
}
if (!grepl(
  "copykat_input <- full_counts[, sample_cells, drop = FALSE]",
  f16_text,
  fixed = TRUE
)) {
  stop("F1.6的CopyKAT输入未固定为当前样本全部QC后singlet。")
}
if (grepl(
  "copykat_known_normals <- if (infercnv_evaluable) reference_cells",
  f16_text,
  fixed = TRUE
)) {
  stop("CopyKAT仍复用inferCNV reference，可能把外部样本带入CopyKAT。")
}
if (!grepl(
  "target_sample_id != copykat_manifest$reference_sample_id",
  f16_text,
  fixed = TRUE
)) {
  stop("F1.6缺少CopyKAT已知正常细胞必须来自当前样本的运行时检查。")
}

plan_paths <- c(
  file.path(root, "AGENTS.md"),
  file.path(root, "conventions.txt"),
  file.path(root, "project_overview.txt"),
  file.path(root, "胃癌MLMOD亚群主线研究方案.txt"),
  file.path(root, "data", "metadata", "pipeline_parameters.yaml"),
  file.path(root, "reports", "environment_setup", "F1_execution_plan_for_review.md")
)
plan_text <- paste(vapply(plan_paths, read_all, character(1)), collapse = "\n")
if (grepl("F1主归一化保持LogNormalize|主流程采用LogNormalize \\+ Harmony|SCTransform不例行运行", plan_text)) {
  stop("方案中仍残留LogNormalize主线/SCTransform仅敏感性的旧表述。")
}
if (!grepl("SCTransform v2", plan_text, fixed = TRUE) || !grepl("RNA raw counts", plan_text, fixed = TRUE)) {
  stop("方案未完整登记SCTransform主线和RNA raw counts用途。")
}

main_plan <- read_all(file.path(root, "胃癌MLMOD亚群主线研究方案.txt"))
if (length(gregexpr("\\[TABLE", main_plan, perl = TRUE)[[1]]) !=
    length(gregexpr("\\[/TABLE\\]", main_plan, perl = TRUE)[[1]])) {
  stop("主线方案的[TABLE]与[/TABLE]数量不一致。")
}

cat("STATIC CONTRACT PASS：F1脚本可解析；QC→粗注释→DecontX顺序、SCTransform/raw-count路由、inferCNV subclusters与CopyKAT当前样本边界一致。\n")
