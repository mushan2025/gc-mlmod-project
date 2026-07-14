# F0 Execution Plan For Review

Status: prepared for Claude Code review and user approval; not executed.

## Objective

F0 will produce the project-level data reconnaissance and reproducibility baseline required before F1. It does not generate biological conclusions and does not start F1.

## Execution Command

Dry-run review command:

```powershell
& 'C:\Users\14799\AppData\Local\Programs\Python\Python310\python.exe' scripts/F0_setup/run_F0_full_audit.py --project-root .
```

Formal execution command after Claude Code review and user approval:

```powershell
& 'C:\Users\14799\AppData\Local\Programs\Python\Python310\python.exe' scripts/F0_setup/run_F0_full_audit.py --project-root . --execute
```

The wrapper runs four independently reviewable stages in order:

```powershell
& 'C:\Users\14799\AppData\Local\Programs\Python\Python310\python.exe' scripts/F0_setup/F0_step1_structure_and_extract.py --project-root . --execute
& 'C:\Users\14799\AppData\Local\Programs\Python\Python310\python.exe' scripts/F0_setup/F0_step2_sample_info_and_audit.py --project-root . --execute
& 'C:\Users\14799\AppData\Local\Programs\Python\Python310\python.exe' scripts/F0_setup/F0_step3_inventory_and_markers.py --project-root . --execute
& 'C:\Users\14799\AppData\Local\Programs\Python\Python310\python.exe' scripts/F0_setup/F0_step4_decisions_and_gate.py --project-root . --execute
```

Each stage defaults to dry-run without `--execute`, checks its predecessor
outputs, and stops the wrapper on a blocking or pause condition. The four
stages remain one F0 approval unit; they do not create intermediate user gates.

## Selected Environment

Default environment: `laptop_thinkbook16p`, Windows native.

Rationale:
- F0 uses streaming archive, gzip, metadata and manifest checks.
- No GPU is useful for F0.
- No WSL or temporary Linux server is required for the planned F0 script.
- The existing F0 resource assessment estimates peak RAM below 1.5 GB and temporary disk below 3 GB, well within the current local free disk/RAM baseline.

## Inputs

The formal input checklist is in `reports/environment_setup/F0_input_file_checklist.tsv`.

Key inputs:
- `data/public_downloads/GSE183904_RAW.tar`
- `data/public_downloads/GEO_metadata/GSE183904_series_matrix.txt.gz`
- existing manifests under `data/metadata/`
- existing precheck tables under `results/F0_audit/`
- environment baseline files under `environment/`
- completed structure-only GSE239676 preaudit rows in
  `results/F0_audit/predownloaded_resource_structure_audit.tsv`

## Script Behavior

The staged scripts will:
- create/confirm the planned project directories;
- verify required inputs are present before writing outputs;
- compute SHA256 for `GSE183904_RAW.tar`;
- read the tar member list and extract exactly the 40 `*.csv.gz` sample matrices into `data/processed_input/GSE183904/`;
- keep `csv.gz` files compressed and not create permanent plain CSV files;
- parse GSE183904 GEO sample metadata and title-to-patient mapping;
- generate `sample_info.tsv`;
- stream-audit each extracted matrix for orientation, row/column consistency,
  exact missing/noninteger/invalid/negative-value classes, MT/HB gene detection
  and gene-order SHA256;
- compare the formal audit to `gse183904_csv_structure_precheck.tsv`;
- normalize SHA256 comparisons case-insensitively and write every generated
  SHA256 field in uppercase;
- include the GSE239676 structure preaudit without opening biological
  expression results for tuning (8,630 features, 222,240 cells, 20 patients,
  and peritoneal/liver-metastasis labels represented);
- audit marker-panel required fields, duplicate cell types, positive-marker
  format and evidence-ID integrity without modifying the read-only panel;
- create project inventory, file manifest, metadata inventory, author processing audit, readiness-by-F table, external resource inventory, method-prior table, decision evidence log, excluded samples table, F0 gate checklist and F0 reports.

## Pause Conditions

F0 must pause and not advance to F1 if any of the following occur:
- `GSE183904_RAW.tar` is missing, unreadable, SHA256-mismatched, or does not contain 40 `csv.gz` sample matrices;
- GSE183904 GEO metadata cannot be matched to GSM/sample files;
- any sample has `group_analysis = Unclear`;
- formal stream audit conflicts with the precheck table for key fields without an accepted explanation;
- any include-in-F1 sample is not nonnegative integer count-like or has `normalization_artifact_flag = true`;
- required gate outputs are missing.

Marker-panel content issues are nonblocking and yield
`PASS_WITH_NOTED_ISSUES`; they are written only to the conditional
`results/F0_audit/marker_panel_issue_report.tsv`.

`data_audit.tsv` replaces the former row-level numeric approximations with
`integer_check_method`, `total_values_checked`, `missing_value_count`,
`noninteger_float_count`, `invalid_nonnumeric_count`,
`negative_integer_count`, and `numeric_anomaly_count`. Legacy precheck fields
`missing_value_rows` and `invalid_numeric_value_rows` are mapped explicitly to
zero/nonzero compatibility checks and recorded in
`legacy_numeric_precheck_status`.

## Expected Outputs

The formal output list is in `reports/environment_setup/F0_expected_output_manifest.tsv`.

F0 completion requires at minimum:
- `data/metadata/sample_info.tsv`
- `data/metadata/data_audit.tsv`
- `data/metadata/processed_input_manifest.tsv`
- `data/metadata/F0_dataset_inventory.tsv`
- `data/metadata/F0_file_manifest.tsv`
- `data/metadata/F0_metadata_field_inventory.tsv`
- `data/metadata/F0_author_processing_audit.tsv`
- `data/metadata/F0_data_readiness_by_F_section.tsv`
- `data/metadata/F0_external_resource_inventory.tsv`
- `data/metadata/F0_method_prior_decision.tsv`
- `data/metadata/decision_evidence_log.tsv`
- `data/metadata/excluded_samples.tsv`
- `results/F0_audit/F0_gate_checklist.tsv`
- `results/F0_audit/F0_global_data_reconnaissance_report.md`
- `results/F0_audit/F0_execution_report.md`
- `logs/F0_setup/analysis_log.md`

## Review Requests For Claude Code

Please review whether:
- the script faithfully implements the F0 plan without starting F1;
- outputs and fields match the plan and gate requirements;
- the script avoids overwriting researcher-confirmed input files;
- resource estimates remain reasonable for the default laptop;
- sample grouping, Normal_Peritoneum handling and PM n=3 limitations are recorded correctly;
- `data_audit.tsv` logic is strong enough to support F1 intake decisions;
- missing R packages are not treated as F0 blockers but remain F1 blockers where relevant.
