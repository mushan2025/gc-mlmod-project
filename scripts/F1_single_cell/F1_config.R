# F1 统一参数表 --------------------------------------------------------------
#
# 目的：让六个脚本使用同一套QC边界、SCTransform设置、随机种子和文件路径。
# 这里登记的是执行前已经确认的规则；正式分析时不要根据MLMOD分数或结果好坏修改。

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
      globin_panel = c(
        "HBA1", "HBA2", "HBB", "HBD", "HBE1",
        "HBG1", "HBG2", "HBM", "HBQ1", "HBZ"
      )
    ),
    doublet = list(
      scdblfinder_dbr_per_1k = 0.008,
      doubletfinder_pN = 0.25,
      doubletfinder_pcs = 1:20,
      minimum_cells_for_doubletfinder = 200L
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
    markers = list(
      min_pct = 0.25,
      logfc_threshold = 0.25,
      only_positive = TRUE,
      top_n_per_cluster = 30L
    ),
    cnv = list(
      minimum_reference_cells = 50L,
      infercnv_cutoff = 0.1,
      infercnv_threads = 4L,
      infercnv_hmm = FALSE,
      infercnv_denoise = TRUE,
      copykat_ngene_chr = 5L,
      copykat_win_size = 25L,
      copykat_ks_cut = 0.1,
      copykat_cores = 4L,
      require_explicit_execution_approval = TRUE
    ),
    packages = list(
      F1.1 = c("Seurat", "SeuratObject", "Matrix", "data.table", "R.utils"),
      F1.2 = c(
        "Seurat", "SeuratObject", "Matrix", "ggplot2", "patchwork",
        "SingleCellExperiment", "SummarizedExperiment", "BiocParallel",
        "scDblFinder", "celda", "DoubletFinder"
      ),
      F1.3 = c(
        "Seurat", "SeuratObject", "sctransform", "glmGamPoi", "harmony",
        "leidenbase", "ggplot2", "patchwork"
      ),
      F1.4 = c("Seurat", "SeuratObject", "ggplot2", "patchwork", "data.table"),
      F1.5 = c(
        "Seurat", "SeuratObject", "sctransform", "glmGamPoi", "harmony", "leidenbase",
        "ggplot2", "patchwork", "data.table"
      ),
      F1.6 = c(
        "Seurat", "SeuratObject", "Matrix", "data.table", "infercnv", "copykat"
      )
    )
  )
}
