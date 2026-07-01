# Download Resources

This folder contains resumable download scripts for large public resources.

Rules:

1. Run these scripts as standalone tasks, not mixed with plan editing.
2. Keep logs under `logs/download_resources/`.
3. Keep download manifests under `reports/download_resources/`.
4. Do not mark a large resource as ready for analysis until size and checksum validation pass.
5. For GEO files without upstream MD5, record local SHA256 and gzip readability.
6. For ENA FASTQ files, validate both expected byte size and ENA MD5.

Examples:

```powershell
powershell -ExecutionPolicy Bypass -File scripts\download_resources\download_gse239676.ps1 -PlanOnly
powershell -ExecutionPolicy Bypass -File scripts\download_resources\download_gse239676.ps1
powershell -ExecutionPolicy Bypass -File scripts\download_resources\download_srp444325_fastq.ps1 -PlanOnly
powershell -ExecutionPolicy Bypass -File scripts\download_resources\download_srp444325_fastq.ps1
```

`download_gse239676.ps1` downloads GEO supplementary files:

- `GSE239676_barcodes.tsv.gz`
- `GSE239676_features.tsv.gz`
- `GSE239676_meta.tsv.gz`
- `GSE239676_count_matrix.mtx.gz`

GSE239676 read method:

- The GEO filenames all end with `.gz`, but the files must be read according to their actual file headers and the manifest `gzip_status`, not by suffix alone.
- `GSE239676_count_matrix.mtx.gz` is a real gzip-compressed MatrixMarket file. Read it with gzip-aware MatrixMarket readers, for example `gzfile()` plus `Matrix::readMM()` in R.
- `GSE239676_barcodes.tsv.gz`, `GSE239676_features.tsv.gz`, and `GSE239676_meta.tsv.gz` currently have no gzip header and are readable UTF-8/plain-text TSV-style files despite the `.gz` suffix. Read these with ordinary text/TSV readers, not `gzfile()`.
- Do not rename or decompress these files during resource preparation. F0/F2 audit code should use `reports/download_resources/GSE239676_download_manifest.tsv` to choose the reader and record the observed `gzip_status`.

`download_srp444325_fastq.ps1` reads:

- `data/public_downloads/SRP444325/SRP444325_ENA_read_run.tsv`

and downloads the main paired-end R1/R2 FASTQ files using resume mode.

SRP444325 note:

- The ENA table contains one extra 115-byte `SRR24947498.fastq.gz` entry in addition to the expected paired-end `SRR24947498_1.fastq.gz` and `SRR24947498_2.fastq.gz`.
- The extra non-R1/R2 file is recorded in the plan but excluded from the main download by default.
- Use `-IncludeNonPairedArtifacts` only if an audit later decides this placeholder file must be fetched for provenance.

Useful options:

```powershell
# Write URL/size/checksum plan without downloading.
powershell -ExecutionPolicy Bypass -File scripts\download_resources\download_srp444325_fastq.ps1 -PlanOnly

# Smoke test: process only the first two FASTQ files.
powershell -ExecutionPolicy Bypass -File scripts\download_resources\download_srp444325_fastq.ps1 -MaxFiles 2
```

Completion criteria:

- GSE239676 is usable only after `reports/download_resources/GSE239676_download_manifest.tsv` shows `size_ok` for all four supplementary files and a readable content status for each file. Acceptable `gzip_status` values are `gzip_read_ok` for true gzip files and `plain_text_read_ok_gz_suffix` for plain-text files that retain a `.gz` suffix on GEO. Any `gzip_read_failed`, `gzip_and_plain_text_read_failed`, `missing`, or `skipped` status must block use until resolved or explicitly approved.
- SRP444325 is usable for F2.1 raw reprocessing only after `reports/download_resources/SRP444325_fastq_download_manifest.tsv` shows `size_ok` and `md5_ok` for all 30 main paired FASTQ files.
