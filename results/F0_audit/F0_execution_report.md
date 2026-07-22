# F0 Execution Report

Run ID: 20260721_230132
Generated at: 2026-07-21T23:14:24+08:00

## Inputs Checked

- GSE183904 archive, GEO metadata, original-paper processing-history source audit
- Existing manifests/prechecks
- GSE239676 structure-preaudit records
- Read-only marker panel plus integrity-only evidence references

## Runtime

- Python 3.10.11; NumPy 2.2.6; Windows-native execution path.

## Execution Note

- The initial run paused because sample39 and sample40 retain R-converted 10x GEM suffixes in their column names (BARCODE.1_sample). The barcode validator was extended to recognize this valid form, the cell-by-gene rejection test remained intact, and the complete rerun passed without changing counts or QC rules.

## Main Observations

- sample_info rows: 40
- data_audit rows: 40
- enter_full_F1_independent_reQC samples: 40
- group counts: {'Normal_Gastric': 10, 'Primary_Tumor': 26, 'Peritoneal_Metastasis': 3, 'Normal_Peritoneum': 1}
- object-eligible / excluded-or-pending / main-group samples: 40 / 0 / 39
- author processing-history rows for GSE183904: 20
- fixed-QC recalculation: F0_fixed_QC_recalculation_pass_all_samples; source-rule pass=158608; final pass/fail=153320/5321; sample1 pilot=pass
- min.cells=3 working feature space: F0_working_feature_space_recalculation_pass_all_samples; retained feature range=12982-21683
- public/paper tissue-cell reconciliation: 158641 / 152423
- unresolved upstream details are carried forward as explicit F1 constraints

## Gate Decision

F0_scRNA_F1_gate: PASS_WITH_NOTED_LIMITATIONS
F0_project_data_inventory_status: partial_with_pending_local_inputs

F1 may start only after Claude Code reviews the executed outputs and the user approves the gate.
