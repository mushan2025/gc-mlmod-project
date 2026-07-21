# F1 静态验证 ---------------------------------------------------------------
#
# 只解析脚本和核对文字契约，不读取真实表达矩阵、不加载缺失R包、不运行分析。

file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), winslash = "/"))
root <- normalizePath(file.path(script_dir, "..", ".."), winslash = "/", mustWork = TRUE)

scripts <- c(
  "F1_config.R", "F1_utils.R", "F1_01_import.R", "F1_02_qc_doublet_ambient.R",
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
  "--approve-cnv-execution"
)
for (pattern in required_patterns) {
  if (!grepl(pattern, combined, perl = TRUE)) stop("F1脚本缺少冻结契约：", pattern)
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

cat("STATIC CONTRACT PASS：F1脚本可解析，SCTransform/raw-count路由与方案一致。\n")
