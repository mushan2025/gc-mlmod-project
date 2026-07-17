# F0 Formal Execution Plan For Review

Status: scripts and plan prepared for Claude Code review and user approval; full F0 has not been executed.

## 1. Objective And Boundary

F0 establishes the data and reproducibility baseline required before F1. It will:

- verify the GSE183904 archive, 40 sample matrices, GEO sample mapping and checksums;
- perform a full-stream structural and numeric audit of every matrix;
- recalculate the approved fixed-QC rule on every currently public cell, only to verify that the rule is computable and to estimate its effect;
- document the known and unknown processing history from tissue acquisition to public export;
- produce the project data inventory, readiness map, method-prior decisions and F0 gate report.

F0 will not remove cells, create a filtered Seurat object, run doublet/ambient algorithms, start F1, or make biological claims. Actual cell exclusion begins only in F1 after the F0 gate is independently reviewed and approved.

## 2. Technical Judgment

### Premise verification

- The local archive is expected to contain 40 gene-by-cell `csv.gz` matrices with 26,571 rows and 158,641 public cells in total.
- Current public inputs do not contain FASTQ files, raw/empty droplets, Cell Ranger cell-calling records, author-excluded barcode lists or complete DoubletFinder parameters.
- The approved QC definition is already frozen in the main plan and therefore must be implemented deterministically rather than selected from F0 outcomes.

### Quantitative evidence

The full-file sample1 pilot gave:

- 2,685 public cells and 26,571 archived feature rows;
- 19,294 features after per-sample `min.cells=3`;
- 2,684 cells passing the source-reported `nFeature` and mitochondrial rules;
- 53 additional cells removed by `nCount>1000` and 0 additional cells removed by `HB_percent<5` in sequential accounting;
- 2,631 cells passing the final fixed rule;
- exact QC globin intersection `HBA1,HBA2,HBB,HBD`;
- false broad-prefix genes `HBEGF,HBP1,HBS1L` excluded from `HB_percent`.

The read-only validation script also tests every inequality boundary on synthetic cells and exact globin matching on a synthetic matrix.
It additionally sends a synthetic 40-sample audit contract through Step3 processing-history construction and the ten-item Step4 gate, so renamed fields cannot silently break downstream decisions.

### Conclusion

F0 is suitable for Windows-native streaming execution on the default laptop. R, WSL, GPU and a temporary server are not required. Full execution remains blocked until Claude Code review and explicit user approval.
Immediately before an approved formal run, available RAM and D-drive space must be measured again; the current resource row is evidence for planning, not permission to ignore changed machine state.

## 3. Inputs

The complete checklist is `reports/environment_setup/F0_input_file_checklist.tsv`. Critical inputs are:

- `data/public_downloads/GSE183904_RAW.tar`
- `data/public_downloads/GEO_metadata/GSE183904_series_matrix.txt.gz`
- `docs/source_verification/GSE183904_processing_history_source_audit.tsv`
- existing manifests under `data/metadata/`
- preregistered structure checks under `results/F0_audit/`
- environment baseline files under `environment/`
- read-only `data/metadata/cell_type_marker_panel.tsv`

No network download is part of this F0 execution.

## 4. Scripts And Stage Order

The wrapper runs four stages as one F0 approval unit:

1. `F0_step1_structure_and_extract.py`
   - verify required paths and archive SHA256;
   - verify exactly 40 unique `csv.gz` members;
   - extract compressed matrices to `data/processed_input/GSE183904/`;
   - write the processed-input manifest with uppercase SHA256 values.

2. `F0_step2_sample_info_and_audit.py`
   - derive sample and group metadata from independent filename and GEO-title evidence;
   - full-stream every value in all 40 matrices without loading a merged matrix;
   - check orientation, dimensions, integer values, missing/invalid/negative values, duplicate genes/barcodes, gene order, sparsity and normalization-artifact sentinels;
   - preserve raw full-matrix `nCount` summaries only for input-shape auditing;
   - construct one per-sample `min.cells=3` QC working feature space;
   - calculate the only QC metric set in that working space;
   - calculate fixed-rule pass/fail counts and run the frozen sample1 regression.

3. `F0_step3_inventory_and_markers.py`
   - build dataset, file, metadata and external-resource inventories;
   - record uppercase SHA256 for every formal F0 script and the read-only validator in `F0_file_manifest.tsv`;
   - audit marker-panel structure without modifying the panel;
   - combine primary-source processing history with F0 observations;
   - keep author reports, F0 recalculation, cross-source inference and unresolved history distinct.

4. `F0_step4_decisions_and_gate.py`
   - write readiness, method-prior, decision-evidence and exclusion tables;
   - evaluate a ten-item F0 gate checklist;
   - write the global reconnaissance and execution reports;
   - stop without starting F1.

All stages default to a non-writing dry run unless `--execute` is supplied. The wrapper stops at the first exception or blocking condition.

## 5. Fixed-QC Recalculation Contract

### Immutable raw input

All 26,571 archived rows remain unchanged. `raw_full_nCount_*` is calculated from those rows only to assess count-matrix shape and possible prior normalization; it never determines cell removal.

### One QC working feature space

For each sample independently, retain features detected in at least 3 cells. In this one working space calculate:

- `nCount_RNA`
- `nFeature_RNA`
- `mt_percent`
- `HB_percent`

There is no second cell-QC metric space.

### Fixed cell rule

A cell passes only when all conditions hold:

```text
500 <= nFeature_RNA < 6000
nCount_RNA > 1000
mt_percent <= 20
HB_percent < 5
```

The `nFeature` and mitochondrial conditions are source-reported. The `nCount` and globin conditions are project additions. Per-rule failure counts and the unique final union are both recorded so overlapping failures are not double-counted.

### Frozen globin definition

`HB_percent` uses exact, case-normalized matching to:

```text
HBA1,HBA2,HBB,HBD,HBE1,HBG1,HBG2,HBM,HBQ1,HBZ
```

For every sample the audit records expected genes, genes present in the archived matrix, genes retained in the `min.cells=3` working space, missing genes and excluded false `^HB` matches. Broad prefix matching is forbidden.

### Downstream boundary

The per-sample `min.cells=3` rule is not a permanent gene filter for later differential expression, pseudobulk or scoring. Those methods must use their own approved expression-coverage rules.

## 6. Processing-History Interpretation

The audit may support all of the following simultaneously:

- the source reports Cell Ranger v3.0/hg38 and Seurat QC thresholds;
- the public matrices are nonnegative-integer called/retained-cell count matrices;
- the 158,641 public cells differ from the 152,423 cells in the paper's final 40-tissue-sample object;
- the exact timing and barcode effect of author QC, DoubletFinder and public export remain unresolved.

F0 must not attribute the 6,218-cell difference to one specific step. Recalculated pass counts describe what the approved rule does to the currently public matrices; they do not reconstruct cells that were removed before public export.

In `F0_author_processing_audit.tsv`, `author_reported_status` is reserved for what a reviewed source says. F0 recomputation and cross-source reconciliation outcomes are stored separately in `record_status`; generated F0 rows use an explicit `not_applicable_*` value in the author-status field.

## 7. Commands

Read-only input/output-contract dry run:

```powershell
& 'C:\Users\14799\AppData\Local\Programs\Python\Python310\python.exe' scripts/F0_setup/run_F0_full_audit.py --project-root .
```

Read-only code validation, including the real sample1 regression:

```powershell
& 'C:\Users\14799\AppData\Local\Programs\Python\Python310\python.exe' scripts/F0_setup/validate_F0_readonly.py --project-root .
```

Formal execution, permitted only after Claude Code review and user approval:

```powershell
& 'C:\Users\14799\AppData\Local\Programs\Python\Python310\python.exe' scripts/F0_setup/run_F0_full_audit.py --project-root . --execute
```

## 8. Blocking Conditions

F0 pauses if any of the following occurs:

- archive missing, unreadable, checksum mismatch, duplicate member basename, or member count other than 40;
- archive filename, GSM accession, `sample_id` and GEO-title mapping disagree;
- any group remains `Unclear`;
- matrix orientation, dimensions, gene order or preregistered structural fields disagree without an approved explanation;
- missing, nonnumeric, noninteger or negative values; duplicate genes/barcodes; or an unevaluable numeric distribution;
- a normalization-artifact sentinel triggers;
- the `min.cells=3` working feature space or any fixed-QC metric is not evaluable;
- the frozen globin definition changes or no panel gene can be recognized;
- sample1 does not reproduce 19,294 working features, 2,684 source-rule pass cells, 53 additional `nCount` failures, 0 additional HB failures and 2,631 final pass cells;
- required processing-history stages or downstream constraints are missing;
- any formal F0 script is missing from `F0_file_manifest.tsv` or lacks a valid uppercase SHA256;
- any required F0 contract table or report is not generated.

Known unavailable upstream information is nonblocking only when explicitly labeled and mapped to a conservative F1 action. It yields `PASS_WITH_NOTED_LIMITATIONS`, not a silent clean pass.

## 9. F1 Constraints Carried Forward

If F0 is approved:

- F1 preserves the full raw matrices and independently applies the same one-space fixed QC rule;
- scDblFinder is the primary per-sample doublet call; DoubletFinder is sensitivity only;
- DecontX is run on filtered raw integer counts after fixed QC and primary doublet removal;
- ambient scores do not delete cells, and raw counts remain the main matrix;
- true barcode knee/cell calling, emptyDrops, SoupX and CellBender remain `not_evaluable_input_limited` without raw droplets;
- Normal_Peritoneum remains reference-only, and PM sample-level inference remains directional because `n=3`.

## 10. Expected Outputs

The formal list is `reports/environment_setup/F0_expected_output_manifest.tsv`. Core outputs include:

- `data/metadata/sample_info.tsv`
- `data/metadata/data_audit.tsv`
- `data/metadata/processed_input_manifest.tsv`
- `data/metadata/F0_author_processing_audit.tsv`
- `data/metadata/F0_data_readiness_by_F_section.tsv`
- `data/metadata/F0_method_prior_decision.tsv`
- `data/metadata/decision_evidence_log.tsv`
- `results/F0_audit/F0_gate_checklist.tsv`
- `results/F0_audit/F0_global_data_reconnaissance_report.md`
- `results/F0_audit/F0_execution_report.md`
- `logs/F0_setup/analysis_log.md`

## 11. Claude Code Review Requests

Please verify that:

- F0 recalculates but does not execute F1 cell deletion;
- the exact inequality boundaries match the main plan;
- all QC metrics are calculated only after per-sample `min.cells=3`;
- raw rows remain intact and are not confused with the QC working space;
- globin matching is exact and auditable;
- sample1 regression is an enforced gate rather than narrative documentation;
- author-reported processing and F0-observed effects are not conflated;
- the ten gate items, output contracts, checksums and pause conditions are closed;
- the laptop resource estimate remains credible;
- no script starts F1 or uses MLMOD/outcome information.
