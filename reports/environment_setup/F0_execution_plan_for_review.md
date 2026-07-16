# F0 Execution Plan For Review

Status: prepared for Claude Code review and user approval; not executed.

Round-2 correction record: commit `59bc74b` implemented exact anomaly
classification but silently omitted the broader all-cell nCount, distribution
and duplicate metrics that had been promised during review. This revision
restores that omitted scope explicitly; the prior narrow summary must not be
treated as evidence that the complete plan had already been implemented.

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
- Step 2 uses the locked project Python 3.10.11 and NumPy 2.2.6 for
  row-wise integer parsing and per-cell accumulation; R is not required.
- The resource assessment has been updated from the strict single-sample
  pilot described below; the expanded audit remains locally feasible.

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
- cross-check, for every GSM accession, the `sample_id` parsed from the step-1
  archive filename against the independent `sampleN` token in the GEO title;
  no positional fallback is allowed, and a missing or unequal value pauses F0;
- generate `sample_info.tsv`;
- stream-audit each extracted matrix for orientation, row/column consistency,
  exact missing/noninteger/invalid/negative-value classes, global min/max,
  integer and zero rates, duplicate genes/barcodes, MT/HB gene detection and
  gene-order SHA256;
- accumulate `nCount` for every cell over all 26,571 gene rows and report its
  min/Q1/median/Q3/max, distinct count and preregistered concentration metrics;
- calculate every per-gene mean and record a numerical sparsity/right-tail
  summary rather than a qualitative assertion;
- compare the formal audit to `gse183904_csv_structure_precheck.tsv`;
- normalize SHA256 comparisons case-insensitively and write every generated
  SHA256 field in uppercase;
- include the GSE239676 structure preaudit without opening biological
  expression results for tuning (8,630 features, 222,240 cells, 20 patients,
  and peritoneal/liver-metastasis labels represented); explicitly carry forward
  that its restricted feature space requires F2.4 signature-coverage and
  fixed-UCell-maxRank assessment and does not support direct cross-cohort
  absolute-score comparison;
- audit marker-panel required fields, duplicate cell types, positive-marker
  format and evidence-ID integrity without modifying the read-only panel;
- create project inventory, file manifest, metadata inventory, author processing audit, readiness-by-F table, external resource inventory, method-prior table, decision evidence log, excluded samples table, F0 gate checklist and F0 reports.

## Pause Conditions

F0 must pause and not advance to F1 if any of the following occur:
- `GSE183904_RAW.tar` is missing, unreadable, SHA256-mismatched, or does not contain 40 `csv.gz` sample matrices;
- GSE183904 GEO metadata cannot be matched to GSM/sample files;
- any sample has `group_analysis = Unclear`;
- formal stream audit conflicts with the precheck table for key fields without an accepted explanation;
- a filename-derived sample ID is missing from or differs from its GEO title-derived sample ID;
- any matrix is not nonnegative integer count-like, contains duplicate gene
  names or within-file cell barcodes, has a non-evaluable distribution, or has
  `normalization_artifact_flag = true`;
- the sparse/right-skew distribution sentinel requires review;
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

The strict implementation additionally records `integer_value_rate`,
`min_value`, `max_value`, `has_negative_value`, duplicate counts, all-cell
`per_cell_nCount_*` fields, a numerical `per_gene_mean_distribution_note`, and
the evidence underlying `normalization_artifact_flag`.

## Preregistered Matrix-Type Semantics

Two evidence layers are kept separate:

- `observed_numeric_type = nonnegative_integer_count_like` is the format-level
  observation and is compared with the legacy precheck column named
  `suspected_matrix_type`.
- `suspected_matrix_type = author_filtered_raw_gene_count_matrix` is written
  only after the format/distribution checks, manifest/GEO mapping, legacy
  precheck and preregistered GEO processing boundary all pass.
- `format_decision_scope = file_format_only`; the final
  `decision_scope = file_format_and_public_processing_boundary`.

This resolves the previous plan contradiction without forcing one field to
carry both storage-format and public-processing semantics.

## Preregistered Normalization-Artifact Rules

These rules were fixed before inspecting the newly computed nCount results.
Any trigger sets `normalization_artifact_flag = true` and pauses Step 2:

1. All cells have exactly the same nCount (`n_distinct = 1`).
2. For a file with at least 200 cells, at least 80% of cells fall within the
   ten most frequent nCount values that are exact multiples of 100.
3. For a file with at least 200 cells, at least 90% of cells are within +/-1%
   of 1,000, 10,000, 100,000 or 1,000,000 and relative IQR is at most 0.02.
4. For a file with at least 200 cells, relative IQR is at most 0.01 and the
   full relative range is at most 0.05, covering a near-constant non-round
   library-size target.

The fixed target rationale is traceable to the official
[Seurat `NormalizeData` reference](https://satijalab.org/seurat/reference/normalizedata),
which documents a default LogNormalize scale factor of 10,000 and an RC/CPM
scale factor of 1,000,000, and the official
[Scanpy `normalize_total` reference](https://scanpy.readthedocs.io/en/stable/generated/scanpy.pp.normalize_total.html),
which documents fixed per-cell total-count normalization. The concentration
percentages and dispersion cutoffs are conservative project sentinels, not
universal biological laws; a trigger requires review and does not identify a
specific normalization package.

`ncount_range_status` separately flags min < 500, median > 100,000 or max >
1,000,000 for review. A range warning alone does not prove normalization and
does not delete cells. If numeric or column anomalies make all-cell nCount
incomplete, the artifact status is `not_evaluable` and Step 2 pauses.

For the descriptive per-gene distribution check, the script reports matrix
zero rate, fractions of gene means <=0.01 and <=0.1, quartiles, mean, maximum
and max/median. `consistent_with_sparse_right_skew` requires zero rate >=0.50,
gene-mean max > Q3 and gene-mean mean > median; otherwise the status is
`review_required` and Step 2 pauses. This sentinel is not by itself proof that
the matrix contains raw counts.

## Post-Preregistration Strict Pilot

After the preceding thresholds were fixed in the plan and code, project Python
3.10.11 plus NumPy 2.2.6 strictly audited
`GSM5573466_sample1.csv.gz` without sampling:

- `total_values_checked = 71,343,135`; elapsed scan time was 4.39 seconds.
- OS polling observed a 34,803,712-byte peak working set.
- Matrix min/max were 0/22,810; all four anomaly counts and both duplicate
  counts were zero; integer rate was 1.
- Per-cell nCount min/Q1/median/Q3/max were
  636/2,585/4,093/6,849/56,631.
- `ncount_distinct_count = 2,323`, relative IQR = 1.04178, relative range =
  13.6807, dominant-round fraction = 0, and no artifact rule triggered.
- Zero-value rate was 0.934476 and the numerical per-gene distribution status
  was `consistent_with_sparse_right_skew`.
- MT/HB counts were 13/10, uppercase gene-order SHA256 matched the legacy
  precheck, and both matrix-type evidence layers reached their preregistered
  expected values.

The 40 matrices contain 4,215,250,011 values, giving a value-count linear
extrapolation of about 4.3 minutes. The registered 5-10 minute estimate adds
file and orchestration overhead. No threshold was changed after this pilot and
the sampled fallback is not requested.

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
