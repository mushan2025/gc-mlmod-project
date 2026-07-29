# F1 统一参数表 --------------------------------------------------------------
#
# 目的：让六个脚本使用同一套QC边界、SCTransform设置、随机种子和文件路径。
# 这里登记的是执行前已经确认的规则；正式分析时不要根据MLMOD分数或结果好坏修改。

f1_env_integer <- function(name, default, minimum = 1L, maximum = 48L) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) return(as.integer(default))
  parsed <- suppressWarnings(as.integer(value))
  if (is.na(parsed) || parsed < minimum || parsed > maximum) {
    stop(
      name, "必须是", minimum, "到", maximum,
      "之间的整数；当前值为：", value
    )
  }
  parsed
}

f1_env_numeric <- function(name, default, minimum, maximum) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) return(as.numeric(default))
  parsed <- suppressWarnings(as.numeric(value))
  if (!is.finite(parsed) || parsed < minimum || parsed > maximum) {
    stop(
      name, "必须是", minimum, "到", maximum,
      "之间的数值；当前值为：", value
    )
  }
  parsed
}

f1_build_config <- function(project_root) {
  root <- normalizePath(project_root, winslash = "/", mustWork = TRUE)

  list(
    project_root = root,
    seed = 42L,
    paths = list(
      processed_manifest = file.path(root, "data", "metadata", "processed_input_manifest.tsv"),
      sample_info = file.path(root, "data", "metadata", "sample_info.tsv"),
      data_audit = file.path(root, "data", "metadata", "data_audit.tsv"),
      f0_gate = file.path(root, "results", "F0_audit", "F0_gate_checklist.tsv"),
      f0_report = file.path(root, "results", "F0_audit", "F0_execution_report.md"),
      marker_panel = file.path(root, "data", "metadata", "cell_type_marker_panel.tsv"),
      gene_order = file.path(
        root, "data", "public_downloads", "inferCNV",
        "infercnv_gene_order_hg38_gencode_v44_gene_symbols_chr1_22.tsv"
      ),
      package_registry = file.path(root, "environment", "F1", "required_packages.tsv"),
      object_dir = file.path(root, "objects", "F1_single_cell"),
      qc_dir = file.path(root, "results", "F1_qc"),
      annotation_dir = file.path(root, "results", "F1_annotation"),
      malignancy_dir = file.path(root, "results", "F1_malignancy"),
      log_dir = file.path(root, "logs", "F1_single_cell"),
      object_01 = file.path(root, "objects", "F1_single_cell", "01_all_cells_raw_or_initial.rds"),
      object_02 = file.path(root, "objects", "F1_single_cell", "02_all_cells_qc_filtered.rds"),
      object_03a = file.path(root, "objects", "F1_single_cell", "03a_all_cells_sct_harmony_clustered.rds"),
      object_03 = file.path(root, "objects", "F1_single_cell", "03_all_cells_integrated_annotated.rds"),
      object_04 = file.path(root, "objects", "F1_single_cell", "04_epithelial_reclustered.rds"),
      object_05 = file.path(root, "objects", "F1_single_cell", "05_malignant_epithelial.rds"),
      object_06a = file.path(root, "objects", "F1_single_cell", "06a_malignant_epithelial_main.rds"),
      object_06b = file.path(root, "objects", "F1_single_cell", "06b_malignant_epithelial_high_confidence_only.rds"),
      decontx_corrected_dir = file.path(root, "objects", "F1_single_cell", "decontX_corrected_by_sample"),
      ambient_cell_estimates = file.path(root, "results", "F1_qc", "ambient_rna_cell_estimates.tsv"),
      ambient_summary = file.path(root, "results", "F1_qc", "ambient_rna_summary_by_sample.tsv"),
      annotation_template = file.path(root, "results", "F1_annotation", "F1_cluster_annotation_template.tsv"),
      annotation_approved = file.path(root, "data", "metadata", "F1_cluster_annotation_approved.tsv"),
      epithelial_review_template = file.path(root, "results", "F1_annotation", "F1_epithelial_cluster_review_template.tsv"),
      epithelial_review_approved = file.path(root, "data", "metadata", "F1_epithelial_cluster_review_approved.tsv"),
      malignancy_review_template = file.path(root, "results", "F1_malignancy", "F1_malignancy_cluster_review_template.tsv"),
      malignancy_review_approved = file.path(root, "data", "metadata", "F1_malignancy_cluster_review_approved.tsv"),
      analysis_log = file.path(root, "logs", "F1_single_cell", "analysis_log.md")
    ),
    qc = list(
      min_cells_per_feature = 3L,
      nfeature_min_inclusive = 500L,
      nfeature_max_exclusive = 6000L,
      ncount_min_exclusive = 1000L,
      mt_max_inclusive = 20,
      hb_max_exclusive = 5,
      no_ncount_sensitivity_enabled = TRUE,
      globin_panel = c(
        "HBA1", "HBA2", "HBB", "HBD", "HBE1",
        "HBG1", "HBG2", "HBM", "HBQ1", "HBZ"
      )
    ),
    doublet = list(
      scdblfinder_dbr_per_1k = 0.008,
      scdblfinder_workers = f1_env_integer(
        "F1_SCDBLFINDER_WORKERS", 1L, minimum = 1L, maximum = 16L
      ),
      doubletfinder_pN = 0.25,
      doubletfinder_pcs = 1:20,
      minimum_cells_for_doubletfinder = 200L
    ),
    ambient = list(
      minimum_reliable_lineages = 2L,
      non_lineage_labels = c("Uncertain", "Mixed_or_doublet_suspect")
    ),
    sct = list(
      vst_flavor = "v2",
      variable_features_n = 3000L,
      vars_to_regress = NULL,
      conserve_memory = TRUE,
      pca_npcs = 50L,
      main_dims = 1:30,
      resolutions = c(0.2, 0.4, 0.6, 0.8, 1.0),
      default_resolution = 0.6,
      leiden_algorithm = 4L,
      harmony_group = "sample_id"
    ),
    execution = list(
      # 默认保持单进程；服务器启动脚本显式提高并行度，并把实际值写入参数表。
      future_workers = f1_env_integer(
        "F1_FUTURE_WORKERS", 1L, minimum = 1L, maximum = 12L
      ),
      future_globals_max_gb = f1_env_numeric(
        "F1_FUTURE_GLOBALS_MAX_GB", 32, minimum = 4, maximum = 80
      )
    ),
    markers = list(
      min_pct = 0.25,
      logfc_threshold = 0.25,
      only_positive = TRUE,
      top_n_per_cluster = 30L
    ),
    cnv = list(
      minimum_reference_cells = 50L,
      maximum_reference_cells = 500L,
      infercnv_cutoff = 0.1,
      infercnv_threads = f1_env_integer(
        "F1_INFERCNV_THREADS", 4L, minimum = 1L, maximum = 32L
      ),
      infercnv_hmm = FALSE,
      infercnv_denoise = TRUE,
      infercnv_analysis_mode = "subclusters",
      infercnv_tumor_subcluster_partition_method = "leiden",
      infercnv_k_nn = 20L,
      infercnv_leiden_resolution = "auto",
      infercnv_leiden_method = "PCA",
      infercnv_leiden_function = "CPM",
      infercnv_inspect_subclusters = TRUE,
      infercnv_output_format = "pdf",
      infercnv_resume_mode = FALSE,
      copykat_ngene_chr = 5L,
      copykat_win_size = 25L,
      copykat_ks_cut = 0.1,
      copykat_cores = f1_env_integer(
        "F1_COPYKAT_CORES", 4L, minimum = 1L, maximum = 32L
      ),
      require_explicit_execution_approval = TRUE
    ),
    packages = list(
      F1.1 = c("Seurat", "SeuratObject", "Matrix", "data.table", "R.utils"),
      F1.2 = c(
        "Seurat", "SeuratObject", "Matrix", "ggplot2", "patchwork",
        "SingleCellExperiment", "SummarizedExperiment", "BiocParallel",
        "scDblFinder", "DoubletFinder"
      ),
      F1.3 = c(
        "Seurat", "SeuratObject", "sctransform", "glmGamPoi", "harmony",
        "leidenbase", "future", "ggplot2", "patchwork", "data.table"
      ),
      F1.4 = c(
        "Seurat", "SeuratObject", "Matrix", "SingleCellExperiment",
        "SummarizedExperiment", "celda", "presto", "ggplot2", "patchwork",
        "data.table"
      ),
      F1.5 = c(
        "Seurat", "SeuratObject", "sctransform", "glmGamPoi", "harmony",
        "leidenbase", "future",
        "ggplot2", "patchwork", "data.table"
      ),
      F1.6 = c(
        "Seurat", "SeuratObject", "Matrix", "data.table", "infercnv", "copykat"
      )
    )
  )
}
