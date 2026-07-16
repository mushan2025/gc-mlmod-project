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
- the preregistered source-evidence table
  `docs/source_verification/GSE183904_processing_history_source_audit.tsv`
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
- recompute per-cell `nCount`, `nFeature` and mitochondrial-transcript
  percentage in two explicit spaces: `public_full_feature_space`, which
  describes the archived matrix and feeds independent F1 QC, and
  `author_like_min_cells3_feature_space`, which first retains only features
  detected in at least three cells and is used solely to test the reported
  `500 <= nFeature < 6000` and `percent.mt <= 20` provenance boundary;
- count positive cells for every feature in every sample and compare public
  rows directly with the paper's per-sample `min.cells = 3` analysis rule;
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
- build an end-to-end GSE183904 processing lineage from tissue acquisition
  through public-matrix export and the author's later analysis, keeping four
  evidence classes separate: author report, F0 observation, cross-source
  inference and unresolved history;
- reconcile the 158,641 public cells with the paper's 152,423-cell final
  40-tissue-sample object without attributing the 6,218-cell difference to a
  specific step when barcode-level evidence is unavailable;
- treat author SCTransform, anchor integration, clustering and annotation as
  post-export analysis context rather than transformations embedded in the
  public raw-count CSV values;
- convert each unresolved matrix-relevant history item into an explicit F1
  constraint, including the six-layer per-sample QC framework, independent
  per-capture doublet assessment, retained-cell DecontX diagnostics, and
  `not_evaluable_input_limited` status for true knee/cell calling, emptyDrops,
  SoupX and CellBender;
- create project inventory, file manifest, metadata inventory, author processing audit, readiness-by-F table, external resource inventory, method-prior table, decision evidence log, excluded samples table, F0 gate checklist and F0 reports.

## Pause Conditions

F0 must pause and not advance to F1 if any of the following occur:
- `GSE183904_RAW.tar` is missing, unreadable, SHA256-mismatched, or does not contain 40 `csv.gz` sample matrices;
- GSE183904 GEO metadata cannot be matched to GSM/sample files;
- any sample has `group_analysis = Unclear`;
- formal stream audit conflicts with the precheck table for key fields without an accepted explanation;
- either cell-QC feature space is non-evaluable, or a measured threshold
  mismatch cannot be localized and reported; an evaluable mismatch is a
  nonblocking provenance limitation and does not authorize silent deletion;
- the per-sample detected-cell count cannot be evaluated for every feature row;
- required processing-history stages are missing, core provenance fields are
  empty, or an unresolved matrix-relevant item lacks a conservative downstream
  action;
- a filename-derived sample ID is missing from or differs from its GEO title-derived sample ID;
- any matrix is not nonnegative integer count-like, contains duplicate gene
  names or within-file cell barcodes, has a non-evaluable distribution, or has
  `normalization_artifact_flag = true`;
- the sparse/right-skew distribution sentinel requires review;
- required gate outputs are missing.

Marker-panel content issues are nonblocking and yield
`PASS_WITH_NOTED_ISSUES`; they are written only to the conditional
`results/F0_audit/marker_panel_issue_report.tsv`.

Genuine residual provenance unknowns are also nonblocking only when they are
explicitly evidenced and mapped to an F1 constraint. Their presence makes the
overall result `PASS_WITH_NOTED_LIMITATIONS`, not a clean `PASS`. An unknown
must never be converted into an assertion that the author did or did not apply
the step.

Public feature rows below the author's per-sample three-cell rule are a
measured export-boundary mismatch, not a structural count failure. When fully
evaluable, this yields `PASS_WITH_NOTED_ISSUES` and requires F1 to preserve all
archived rows, use the author-like space only for provenance, and define each
downstream feature-eligibility rule separately. It does not authorize deleting
cells, permanently deleting genes, or modifying archived counts.

A cell-threshold mismatch that appears only after the author-like min.cells=3
feature restriction is likewise a measured, space-dependent processing-boundary
result. It is reported in both spaces and yields `PASS_WITH_NOTED_ISSUES` when
fully evaluable; it does not invalidate the integer counts or force exact
reproduction of the author's object.

`data_audit.tsv` replaces the former row-level numeric approximations with
`integer_check_method`, `total_values_checked`, `missing_value_count`,
`noninteger_float_count`, `invalid_nonnumeric_count`,
`negative_integer_count`, and `numeric_anomaly_count`. Legacy precheck fields
`missing_value_rows` and `invalid_numeric_value_rows` are mapped explicitly to
zero/nonzero compatibility checks and recorded in
`legacy_numeric_precheck_status`.

The strict implementation additionally records `integer_value_rate`,
`min_value`, `max_value`, `has_negative_value`, duplicate counts, all-cell
`per_cell_nCount_*`, `public_full_feature_*` and `author_like_*` fields, the
space-specific counts mismatching each author-reported cell-QC threshold,
the number of cells whose nFeature/nCount/percent.mt changes between spaces and
the maximum observed change,
per-sample counts of features detected in 0, 1 or 2 cells, the resulting
comparison with the author-reported `min.cells = 3` feature rule, a numerical
`per_gene_mean_distribution_note`, and the evidence underlying
`normalization_artifact_flag`.

## Preregistered Matrix-Type Semantics

Two evidence layers are kept separate:

- `observed_numeric_type = nonnegative_integer_count_like` is the format-level
  observation and is compared with the legacy precheck column named
  `suspected_matrix_type`.
- `author_cell_qc_reproduction_status_public_space` and
  `author_cell_qc_reproduction_status_author_like_space` separately report
  whether retained public cells satisfy `500 <= nFeature < 6000` and
  `percent.mt <= 20`; the overall status is `pass`, `measured_mismatch` or
  `not_evaluable`.
- `author_feature_filter_reproduction_status` independently reports whether
  every public sample-by-feature row satisfies the paper's per-sample
  three-cell feature rule. A measured mismatch is a provenance limitation,
  not a count-format failure.
- `suspected_matrix_type = public_called_cell_raw_gene_count_matrix` records
  only the public input shape. Cell-threshold and feature-filter boundaries are
  never compressed into this field.
- `audit_decision = enter_full_F1_independent_reQC` means that the matrix can
  enter full downstream analysis after an independent project QC; it does not
  assert exact author-object reproduction.
- `format_decision_scope = file_format_only`; the final
  `decision_scope = file_format_and_public_processing_boundary`.

This resolves the previous plan contradiction without forcing one field to
carry storage format, cell filtering and feature filtering simultaneously.

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

- `total_values_checked = 71,343,135`; the expanded repeat scan took 4.73 seconds.
- OS polling observed a 35,151,872-byte peak working set.
- Matrix min/max were 0/22,810; all four anomaly counts and both duplicate
  counts were zero; integer rate was 1.
- Per-cell nCount min/Q1/median/Q3/max were
  636/2,585/4,093/6,849/56,631.
- Per-cell nFeature min/Q1/median/Q3/max were
  500/1,070/1,481/2,033/5,972; no cell was below 500 or at least 6,000.
- Per-cell percent.mt min/Q1/median/Q3/max were
  0/5.6228/8.2662/11.8887/19.9971; no cell exceeded 20%, so all 2,685
  cells passed the recomputed author-reported cell-QC baseline in
  `public_full_feature_space`.
- The public file contained 5,770 genes detected in 0 cells, 685 in 1 cell and
  822 in 2 cells: 7,277 feature rows below the paper's per-sample three-cell
  rule. After restricting to the resulting 19,294 author-like features,
  nFeature min/Q1/median/Q3/max became 500/1,070/1,480/2,033/5,969; 1,197 cells
  changed nFeature (maximum decrease 41), and percent.mt max became 20.008763,
  producing one near-boundary mismatch. Thus the archived public space broadly
  supports the reported cell thresholds for this sample, but the exact result
  depends on feature-space/order and cannot be represented by one compliance
  field.
- `ncount_distinct_count = 2,323`, relative IQR = 1.04178, relative range =
  13.6807, dominant-round fraction = 0, and no artifact rule triggered.
- Zero-value rate was 0.934476 and the numerical per-gene distribution status
  was `consistent_with_sparse_right_skew`.
- MT/HB counts were 13/10 and uppercase gene-order SHA256 matched the legacy
  precheck. The raw-count format layer passed; the feature-filter comparison
  exposed why the former broad `author_filtered_raw_gene_count_matrix` label
  must be replaced by a format state plus two processing-boundary states.

The 40 matrices contain 4,215,250,011 values, giving a value-count linear
extrapolation of about 4.7 minutes. The registered 5-10 minute estimate adds
file and orchestration overhead. No biological or anomaly threshold was tuned
to the pilot; the added three-cell comparison directly implements the author's
reported rule. The sampled fallback is not requested.

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
- both cell-QC feature spaces are measured independently and an evaluable
  mismatch remains nonblocking but visible;
- missing R packages are not treated as F0 blockers but remain F1 blockers where relevant.
- processing-history claims are correctly attributed to primary sources, F0
  observations or bounded cross-source inference;
- every remaining Cell Ranger, doublet, ambient-RNA and public-export unknown
  has an adequate F1 constraint without claiming exact reproduction of the
  author's final object.
