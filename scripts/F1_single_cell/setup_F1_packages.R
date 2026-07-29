# Ubuntu服务器F1 R包环境 ------------------------------------------------------
#
# 用法：
#   /opt/R/4.4.3/bin/Rscript setup_F1_packages.R --project-root=/path/to/project
#
# 只安装F1真正使用的包，并尽量固定到本地已审核环境的直接依赖版本。
# Linux系统依赖由setup_F1_ubuntu22.sh负责；本脚本不修改分析参数。

args <- commandArgs(trailingOnly = TRUE)
root_arg <- grep("^--project-root=", args, value = TRUE)
if (length(root_arg) != 1L) {
  stop("必须提供唯一的--project-root=/path/to/project。")
}
project_root <- normalizePath(
  sub("^--project-root=", "", root_arg),
  winslash = "/",
  mustWork = TRUE
)

if (!identical(as.character(getRversion()), "4.4.3")) {
  stop("F1服务器环境要求R 4.4.3；当前为", getRversion(), "。")
}
if (!identical(.Platform$OS.type, "unix")) {
  stop("本脚本只用于Linux服务器。")
}

ncpus <- suppressWarnings(as.integer(Sys.getenv("F1_INSTALL_NCPUS", "12")))
if (is.na(ncpus) || ncpus < 1L || ncpus > 24L) {
  stop("F1_INSTALL_NCPUS必须为1到24之间的整数。")
}
cran_repo <- Sys.getenv(
  "F1_CRAN_REPO",
  "https://mirrors.tuna.tsinghua.edu.cn/CRAN"
)
bioc_mirror <- Sys.getenv(
  "F1_BIOC_MIRROR",
  "https://bioconductor.org"
)
options(
  repos = c(CRAN = cran_repo),
  BioC_mirror = bioc_mirror,
  Ncpus = ncpus,
  timeout = 3600
)
Sys.setenv(MAKEFLAGS = paste0("-j", max(1L, floor(48L / ncpus))))

install_if_missing <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) {
    install.packages(
      package,
      dependencies = c("Depends", "Imports", "LinkingTo"),
      Ncpus = ncpus
    )
  }
}

install_if_missing("BiocManager")
install_if_missing("remotes")
BiocManager::install(version = "3.20", ask = FALSE, update = FALSE)

# glmGamPoi 1.18.0仍按C++11构建；15.4系列RcppArmadillo已要求C++14。
# 固定到本地已验证的15.2.4.1，避免在Linux源码编译时出现标准不兼容。
rcpp_armadillo_version <- "15.2.4.1"
installed_rcpp_armadillo <- if (
  requireNamespace("RcppArmadillo", quietly = TRUE)
) {
  as.character(utils::packageVersion("RcppArmadillo"))
} else {
  ""
}
if (!identical(installed_rcpp_armadillo, rcpp_armadillo_version)) {
  remotes::install_version(
    "RcppArmadillo",
    version = rcpp_armadillo_version,
    dependencies = c("Depends", "Imports", "LinkingTo"),
    upgrade = "never",
    Ncpus = ncpus
  )
}

bioc_packages <- c(
  "glmGamPoi", "SingleCellExperiment", "SummarizedExperiment",
  "BiocParallel", "scDblFinder", "celda", "infercnv"
)
missing_bioc <- bioc_packages[
  !vapply(bioc_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_bioc)) {
  BiocManager::install(
    missing_bioc,
    ask = FALSE,
    update = FALSE,
    Ncpus = ncpus
  )
}

cran_versions <- c(
  Seurat = "5.5.0",
  SeuratObject = "5.4.0",
  sctransform = "0.4.3",
  harmony = "1.2.4",
  Matrix = "1.7-5",
  data.table = "1.18.2.1",
  R.utils = "2.13.0",
  ggplot2 = "4.0.3",
  patchwork = "1.3.2",
  leidenbase = "0.1.36",
  future = "1.70.0"
)
for (package in names(cran_versions)) {
  wanted <- cran_versions[[package]]
  installed <- if (requireNamespace(package, quietly = TRUE)) {
    as.character(utils::packageVersion(package))
  } else {
    ""
  }
  # R把1.7-5显示为1.7.5；比较时统一标点。
  normalized_match <- gsub("-", ".", installed, fixed = TRUE) ==
    gsub("-", ".", wanted, fixed = TRUE)
  if (!normalized_match) {
    remotes::install_version(
      package,
      version = wanted,
      dependencies = c("Depends", "Imports", "LinkingTo"),
      upgrade = "never",
      Ncpus = ncpus
    )
  }
}

github_specs <- c(
  DoubletFinder =
    "chris-mcginnis-ucsf/DoubletFinder@1B244D8F0D54B4B1CB4365639931BBB16F01E1CD",
  copykat =
    "navinlabcode/copykat@12B7C7E15D42596296E46819C64ACA347BFDE2E5"
)
github_versions <- c(DoubletFinder = "2.0.6", copykat = "1.1.0")
for (package in names(github_specs)) {
  installed <- if (requireNamespace(package, quietly = TRUE)) {
    as.character(utils::packageVersion(package))
  } else {
    ""
  }
  if (!identical(installed, github_versions[[package]])) {
    remotes::install_github(
      github_specs[[package]],
      dependencies = c("Depends", "Imports", "LinkingTo"),
      upgrade = "never",
      Ncpus = ncpus
    )
  }
}

registry_path <- file.path(project_root, "environment", "F1", "required_packages.tsv")
registry <- read.delim(
  registry_path,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)
packages <- setdiff(registry$package, "R")
versions <- vapply(packages, function(package) {
  if (!requireNamespace(package, quietly = TRUE)) {
    return("NOT_AVAILABLE")
  }
  as.character(utils::packageVersion(package))
}, character(1))
if (any(versions == "NOT_AVAILABLE")) {
  stop(
    "仍缺少F1必需包：",
    paste(packages[versions == "NOT_AVAILABLE"], collapse = ", ")
  )
}

log_dir <- file.path(project_root, "logs", "F1_single_cell")
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
version_table <- data.frame(
  package = c("R", packages),
  version = c(as.character(getRversion()), versions),
  platform = R.version$platform,
  checked_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  stringsAsFactors = FALSE
)
write.table(
  version_table,
  file.path(log_dir, "server_environment_versions.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
capture.output(
  sessionInfo(),
  file = file.path(log_dir, "server_environment_sessionInfo.txt")
)
cat("F1服务器R环境检查通过，共", length(packages), "个直接包。\n", sep = "")
