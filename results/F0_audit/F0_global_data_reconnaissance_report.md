# F0 Global Data Reconnaissance Report

Run ID: 20260721_230132
Generated at: 2026-07-21T23:14:24+08:00

## Current Usable Data

- GSE183904: 40 Cell Ranger-derived public called/retained-cell raw gene-count CSV.gz matrices; all 26571 archived rows are preserved, and one per-sample min.cells=3 working feature space is used for fixed-QC recalculation; group counts: {'Normal_Gastric': 10, 'Primary_Tumor': 26, 'Peritoneal_Metastasis': 3, 'Normal_Peritoneum': 1}.
- GSE183904 is the only dataset eligible to start F1 after F0 gate approval.
- GSE239676 passed a structure-only preaudit (8,630 features, 222,240 cells, 20 patients, PC/LM labels) but remains isolated until F2.4 approval.
- Other bulk, multiomics and external resources remain candidates until section-specific audits.

## GSE183904 Processing Provenance

- Wet-lab acquisition, dissociation, 10x 5-prime library preparation, HiSeq4000 sequencing, Cell Ranger v3.0/hg38 count generation and author Seurat thresholds are source-audited.
- Fixed-QC recalculation: F0_fixed_QC_recalculation_pass_all_samples; source-reported nFeature/percent.mt rules pass=158608; final source-plus-project fixed-QC pass=153320; fail=5321; sample1 pilot=pass.
- Working feature-space recalculation: F0_working_feature_space_recalculation_pass_all_samples; retained features per sample=12982-21683; archived sample-by-gene rows below 3 detected cells=320543.
- Author doublet method: DoubletFinder was used on gene-expression data and identified doublets were removed from the author's analysis data set.
- Public versus paper final tissue-cell count: 158641 versus 152423; difference=6218.
- Cell Ranger cell-calling details, exact DoubletFinder parameters/barcodes, ambient correction status and the public export timing relative to doublet removal remain unresolved; true knee/cell calling, emptyDrops, SoupX and CellBender are not_evaluable_input_limited.
- Author SCTransform, integration, clustering and annotation describe the author's downstream object and are not embedded in the public integer-count CSV values.

## Boundaries

- F0 does not produce biological conclusions.
- GSE183904 lacks FASTQ/raw/empty droplets in current public inputs.
- Normal_Peritoneum is reference/display only; PM sample-level inference is directional because n=3.
- The fixed-QC counts describe an independent recalculation on currently public cells; they do not prove when the author filtered cells or reproduce excluded barcodes.
- Any formal/precheck mismatch, non-evaluable fixed rule or sample1 regression failure blocks F1.

## F0 Gate Summary

- F0_scRNA_F1_gate: PASS_WITH_NOTED_LIMITATIONS
- F0_project_data_inventory_status: partial_with_pending_local_inputs
- GSE183904 contains 40 samples; 40 are allowed into F1 object construction if approved, 0 are excluded or pending, and 39 are allowed into the main group comparison.
- Main comparison groups are Normal_Gastric, Primary_Tumor and Peritoneal_Metastasis; Normal_Peritoneum is reference-only.
- Blocking checklist failures: 0
