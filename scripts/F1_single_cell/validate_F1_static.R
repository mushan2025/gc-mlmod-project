# F1 静态验证 ---------------------------------------------------------------
#
# 只解析脚本和核对文字契约，不读取真实表达矩阵、不加载缺失R包、不运行分析。

file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), winslash = "/"))
root <- normalizePath(file.path(script_dir, "..", ".."), winslash = "/", mustWork = TRUE)

scripts <- c(
  "F1_config.R", "F1_utils.R", "F1_01_import.R", "F1_02_qc_doublet.R",
  "F1_03_sct_harmony_cluster.R", "F1_04_prepare_annotation_review.R",
  "F1_04_annotation.R", "F1_04_validate_decontx_outputs.R",
  "F1_04_ambient_impact_review.R",
  "F1_05_epithelial_recluster.R",
  "F1_05_validate_and_ambient_review.R",
  "F1_05_prepare_annotation_review.R",
  "F1_06_malignancy_inference.R", "run_F1.R"
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
  "min_cells = config\\$sct\\$minimum_cells_per_gene",
  "features = pca_features",
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
  "reuse.pANN = NULL",
  "vapply\\(config\\$sct\\$resolutions",
  "minimum_reliable_lineages = 2L",
  "f1_cluster_composition",
  "f1_apply_f1_annotation_draft",
  "z = model_clusters",
  "completed_with_researcher_approved_seurat_clusters",
  "researcher_approved_seurat_cluster_partition",
  "minimum_detected_features_for_sample_sct <- 500L",
  "minimum_cells_per_gene_for_sct <- config\\$sct\\$minimum_cells_per_gene",
  "epithelial_sample_input_selection.tsv",
  "epithelial_excluded_cells.tsv",
  "epithelial_ambient_annotation_impact_summary.tsv",
  "infercnv_analysis_mode = \"subclusters\"",
  "infercnv_tumor_subcluster_partition_method = \"leiden\"",
  "infercnv_k_nn = 20L",
  "infercnv_leiden_resolution = \"auto\"",
  "infercnv_observation_group = \"observations\"",
  "infercnv_hmm = TRUE",
  "infercnv_hmm_type = \"i6\"",
  "infercnv_hmm_report_by = \"subcluster\"",
  "infercnv_bayes_max_p_normal = 0.5",
  "infercnv_reassign_cnvs = TRUE",
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
  "A_self_estimated",
  "B_same_sample_immune",
  "C_normal_gastric_holdout_fold",
  "known_normal_cells = character\\(\\)",
  "copykat_args\\$norm.cell.names <- known_normal_cells",
  "copykat_holdout_seed = 42L",
  "not_evaluable_baseline_suspect",
  "copykat_external_reference = \"prohibited\"",
  "infercnv_decontX_corrected",
  "normal_context & !epithelial\\$tumor_program_support",
  "tumor_context & !normal_context",
  "aneuploid & infer_not_strong",
  "CopyKAT单独支持的细胞不得进入06a",
  "06a每个细胞必须至少具有weak inferCNV热图大片段支持",
  "heatmap_broad_segment_description",
  "normal_gastric_epithelial_comparison",
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
if (grepl("norm.cell.names = copykat_known_normals", f16_text, fixed = TRUE)) {
  stop("CopyKAT仍把单一known-normal规则同时当作主分析；A臂必须省略norm.cell.names。")
}
if (!grepl(
  "group = c(\n      rep(config$cnv$infercnv_observation_group, length(obs_cells))",
  f16_text,
  fixed = TRUE
)) {
  stop("inferCNV观察细胞未固定为每样本单一observations组。")
}
if (grepl("paste0(\"obs__\", unname(obs_cluster[obs_cells]))", f16_text, fixed = TRUE)) {
  stop("inferCNV仍按原上皮cluster拆分observation group。")
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
if (!grepl("每个样本的全部候选上皮细胞只写入一个`observations`组", plan_text, fixed = TRUE)) {
  stop("方案未登记inferCNV每样本单一观察组。")
}
if (!grepl("A臂为主分析", plan_text, fixed = TRUE) ||
    !grepl("C臂只用于Normal_Gastric样本", plan_text, fixed = TRUE)) {
  stop("方案未完整登记CopyKAT A/B/C三臂。")
}
main_plan <- read_all(file.path(root, "胃癌MLMOD亚群主线研究方案.txt"))
if (grepl(
  "06a_malignant_epithelial_main.rds\\t固定F2主对象，仅包含.*malignant_probable_copykat",
  main_plan,
  perl = TRUE
)) {
  stop("主线方案仍允许CopyKAT单独支持细胞进入06a。")
}

if (length(gregexpr("\\[TABLE", main_plan, perl = TRUE)[[1]]) !=
    length(gregexpr("\\[/TABLE\\]", main_plan, perl = TRUE)[[1]])) {
  stop("主线方案的[TABLE]与[/TABLE]数量不一致。")
}

cat("STATIC CONTRACT PASS：F1脚本可解析；inferCNV单一观察组+i6 HMM、热图主判、CopyKAT A/B/C三臂及06a inferCNV准入边界一致。\n")
