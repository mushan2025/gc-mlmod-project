# F1 总入口 ------------------------------------------------------------------
#
# 默认只显示执行计划和缺失条件；只有显式增加 --execute 才按F1.1到F1.6运行。
# 可用 --from=F1.3 --to=F1.5 从已完成的中间对象恢复。

file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), winslash = "/"))
source(file.path(script_dir, "F1_config.R"), encoding = "UTF-8")
source(file.path(script_dir, "F1_utils.R"), encoding = "UTF-8")

args <- f1_parse_args()
config <- f1_build_config(args$project_root)
stage_order <- c("F1.1", "F1.2", "F1.3", "F1.4", "F1.5", "F1.6")
script_map <- c(
  F1.1 = "F1_01_import.R",
  F1.2 = "F1_02_qc_doublet.R",
  F1.3 = "F1_03_sct_harmony_cluster.R",
  F1.4 = "F1_04_annotation.R",
  F1.5 = "F1_05_epithelial_recluster.R",
  F1.6 = "F1_06_malignancy_inference.R"
)
description <- c(
  F1.1 = "导入F0批准的gene-by-cell raw count矩阵",
  F1.2 = "固定QC、scDblFinder与DoubletFinder",
  F1.3 = "逐样本SCTransform v2、PCA、Harmony和Leiden",
  F1.4 = "主要谱系注释后逐样本DecontX评估",
  F1.5 = "上皮提取、SCTransform v2二次聚类",
  F1.6 = "raw-count inferCNV/CopyKAT与联合恶性判定"
)

if (!args$from %in% stage_order || !args$to %in% stage_order) {
  stop("--from和--to必须位于F1.1至F1.6。")
}
from_index <- match(args$from, stage_order)
to_index <- match(args$to, stage_order)
if (from_index > to_index) stop("--from不能晚于--to。")
selected <- stage_order[from_index:to_index]

cat("\nF1 SCTransform正式流程\n")
cat("项目根目录：", config$project_root, "\n", sep = "")
cat("本次范围：", paste(selected, collapse = " -> "), "\n\n", sep = "")
for (stage in selected) {
  packages <- config$packages[[stage]]
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  cat(
    sprintf(
      "%s  %-42s  依赖：%s\n",
      stage,
      description[[stage]],
      if (length(missing)) paste0("缺少 ", paste(missing, collapse = ", ")) else "就绪"
    )
  )
}
f0 <- f1_check_f0_ready(config, stop_on_failure = FALSE)
cat("\nF0准入：", if (f0$ready) "就绪" else "尚未就绪", "\n", sep = "")
if (length(f0$problems)) for (problem in f0$problems) cat("  - ", problem, "\n", sep = "")

if (!args$execute) {
  cat("\n当前为只读计划模式：没有读取表达矩阵，也没有写分析结果。\n")
  cat("正式执行需显式增加 --execute；F1.6首次CNV计算还需 --approve-cnv-execution。\n")
  quit(save = "no", status = 0)
}

if (!f0$ready) stop("F0正式输出尚未允许进入F1。")
required_packages <- unique(unlist(config$packages[selected], use.names = FALSE))
f1_require_packages(required_packages, paste(selected, collapse = "-"))

rscript_candidates <- c(
  file.path(R.home("bin"), "Rscript.exe"),
  file.path(R.home("bin"), "x64", "Rscript.exe"),
  Sys.which("Rscript")
)
rscript_candidates <- unique(rscript_candidates[nzchar(rscript_candidates)])
rscript_path <- rscript_candidates[file.exists(rscript_candidates)][[1]]

for (stage in selected) {
  script <- file.path(script_dir, script_map[[stage]])
  command_args <- c(
    shQuote(script),
    "--execute",
    shQuote(paste0("--project-root=", config$project_root))
  )
  if (stage == "F1.6" && args$approve_cnv_execution) {
    command_args <- c(command_args, "--approve-cnv-execution")
  }
  cat("\n开始", stage, "：", description[[stage]], "\n")
  status <- system2(rscript_path, args = command_args)
  if (!identical(status, 0L)) {
    stop(stage, "未完成。若提示等待人工注释/恶性审核，这是预期科学审核点；完成批准表后从本节恢复。")
  }
}
cat("\n所选F1步骤已完成。\n")
