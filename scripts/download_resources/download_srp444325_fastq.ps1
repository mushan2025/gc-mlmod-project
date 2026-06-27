param(
    [string]$ProjectRoot = (Resolve-Path (Join-Path $PSScriptRoot "..\..")).Path,
    [string]$EnaRunTable = "",
    [switch]$Force,
    [switch]$PlanOnly,
    [int]$MaxFiles = 0,
    [switch]$IncludeNonPairedArtifacts
)

$ErrorActionPreference = "Stop"

$dataset = "SRP444325"
$outDir = Join-Path $ProjectRoot "data\public_downloads\SRP444325\fastq"
$metaDir = Join-Path $ProjectRoot "data\public_downloads\SRP444325"
$logDir = Join-Path $ProjectRoot "logs\download_resources\SRP444325"
$reportDir = Join-Path $ProjectRoot "reports\download_resources"
$manifest = if ($PlanOnly) {
    Join-Path $reportDir "SRP444325_fastq_download_plan.tsv"
} else {
    Join-Path $reportDir "SRP444325_fastq_download_manifest.tsv"
}

New-Item -ItemType Directory -Force -Path $outDir, $metaDir, $logDir, $reportDir | Out-Null

if ([string]::IsNullOrWhiteSpace($EnaRunTable)) {
    $EnaRunTable = Join-Path $metaDir "SRP444325_ENA_read_run.tsv"
}

if (-not (Test-Path $EnaRunTable)) {
    throw "ENA run table not found: $EnaRunTable. Download it first from ENA Portal API."
}

function Write-Log {
    param([string]$Message)
    $stamp = Get-Date -Format "yyyy-MM-dd HH:mm:ss"
    $line = "[$stamp] $Message"
    Write-Host $line
    Add-Content -Path (Join-Path $logDir "download.log") -Value $line -Encoding UTF8
}

function Normalize-FtpUrl {
    param([string]$FtpPath)
    if ($FtpPath -match '^https?://') { return $FtpPath }
    if ($FtpPath -match '^ftp://') { return $FtpPath -replace '^ftp://', 'https://' }
    return "https://$FtpPath"
}

$runs = Import-Csv -Path $EnaRunTable -Delimiter "`t"
$jobs = New-Object System.Collections.Generic.List[object]

foreach ($run in $runs) {
    $urls = @($run.fastq_ftp -split ';' | Where-Object { -not [string]::IsNullOrWhiteSpace($_) })
    $md5s = @($run.fastq_md5 -split ';' | Where-Object { -not [string]::IsNullOrWhiteSpace($_) })
    $bytes = @($run.fastq_bytes -split ';' | Where-Object { $_ -match '^\d+$' })

    if (($urls.Count -ne $md5s.Count) -or ($urls.Count -ne $bytes.Count)) {
        Write-Log "WARNING: URL/MD5/byte count mismatch for $($run.run_accession): urls=$($urls.Count), md5=$($md5s.Count), bytes=$($bytes.Count)"
    }

    for ($i = 0; $i -lt $urls.Count; $i++) {
        $url = Normalize-FtpUrl $urls[$i]
        $name = Split-Path $url -Leaf
        $includeInMainDownload = $true
        $exclusionReason = ""
        if (($run.library_layout -eq "PAIRED") -and ($name -notmatch '_[12]\.fastq\.gz$')) {
            $includeInMainDownload = $false
            $exclusionReason = "paired_end_non_R1_R2_ena_extra_file_or_placeholder"
        }
        $jobs.Add([pscustomobject]@{
            run_accession = $run.run_accession
            sample_accession = $run.sample_accession
            experiment_accession = $run.experiment_accession
            file_name = $name
            source_url = $url
            expected_md5 = if ($i -lt $md5s.Count) { $md5s[$i] } else { "" }
            expected_bytes = if ($i -lt $bytes.Count) { [int64]$bytes[$i] } else { $null }
            library_layout = $run.library_layout
            library_strategy = $run.library_strategy
            instrument_model = $run.instrument_model
            include_in_main_download = $includeInMainDownload
            exclusion_reason = $exclusionReason
        })
    }
}

$jobsArray = @($jobs | ForEach-Object { $_ })
$jobsToProcess = $jobsArray
if ($MaxFiles -gt 0) {
    $jobsToProcess = @($jobsArray | Select-Object -First $MaxFiles)
}

$rows = New-Object System.Collections.Generic.List[object]
Write-Log "Prepared $($jobs.Count) FASTQ file jobs from $($runs.Count) runs; processing $($jobsToProcess.Count)."

foreach ($job in $jobsToProcess) {
    $out = Join-Path $outDir $job.file_name
    $fileLog = Join-Path $logDir "$($job.file_name).curl.log"
    $existingSize = if (Test-Path $out) { (Get-Item $out).Length } else { 0 }
    $completeBySize = ($job.expected_bytes -and (Test-Path $out) -and ($existingSize -eq $job.expected_bytes))
    $skipByPolicy = (-not $IncludeNonPairedArtifacts) -and (-not $job.include_in_main_download)

    if ($skipByPolicy) {
        Write-Log "SKIP policy-excluded $($job.file_name): $($job.exclusion_reason)"
    } elseif ($PlanOnly) {
        Write-Log "PLAN $($job.file_name); existing=$existingSize expected=$($job.expected_bytes)"
    } elseif ((-not $Force) -and $completeBySize) {
        Write-Log "SKIP size-complete file: $($job.file_name)"
    } else {
        Write-Log "DOWNLOAD/RESUME $($job.file_name); existing=$existingSize expected=$($job.expected_bytes)"
        & curl.exe `
            --location `
            --fail `
            --ssl-no-revoke `
            --retry 20 `
            --retry-delay 60 `
            --connect-timeout 60 `
            --speed-time 300 `
            --speed-limit 1024 `
            --continue-at - `
            --output $out `
            $job.source_url 2>&1 | Tee-Object -FilePath $fileLog -Append
        if ($LASTEXITCODE -ne 0) {
            throw "curl failed for $($job.file_name) with exit code $LASTEXITCODE"
        }
    }

    $localSize = if (Test-Path $out) { (Get-Item $out).Length } else { 0 }
    $sizeStatus = if ($skipByPolicy) { "excluded_from_main_download" } elseif ($job.expected_bytes -and ($localSize -eq $job.expected_bytes)) { "size_ok" } elseif ($job.expected_bytes) { "size_mismatch" } else { "expected_size_missing" }
    $md5 = if ($PlanOnly -or $skipByPolicy) { "" } elseif (Test-Path $out) { (Get-FileHash -Algorithm MD5 $out).Hash.ToLowerInvariant() } else { "" }
    $md5Status = if ($skipByPolicy) { "not_checked_excluded" } elseif ($PlanOnly) { "not_checked_plan_only" } elseif ($job.expected_md5 -and ($md5 -eq $job.expected_md5.ToLowerInvariant())) { "md5_ok" } elseif ($job.expected_md5) { "md5_mismatch" } else { "expected_md5_missing" }
    $sha256 = if ($PlanOnly -or $skipByPolicy) { "" } elseif (Test-Path $out) { (Get-FileHash -Algorithm SHA256 $out).Hash } else { "" }

    $rows.Add([pscustomobject]@{
        dataset_id = $dataset
        run_accession = $job.run_accession
        sample_accession = $job.sample_accession
        experiment_accession = $job.experiment_accession
        file_name = $job.file_name
        source_url = $job.source_url
        local_path = $out.Replace($ProjectRoot + "\", "")
        expected_bytes = $job.expected_bytes
        local_size_bytes = $localSize
        size_status = $sizeStatus
        expected_md5 = $job.expected_md5
        observed_md5 = $md5
        md5_status = $md5Status
        sha256 = $sha256
        library_layout = $job.library_layout
        library_strategy = $job.library_strategy
        instrument_model = $job.instrument_model
        include_in_main_download = $job.include_in_main_download
        exclusion_reason = $job.exclusion_reason
        checked_at = (Get-Date -Format "yyyy-MM-dd HH:mm:ss")
    })

    $rows | Export-Csv -Path $manifest -Delimiter "`t" -NoTypeInformation -Encoding UTF8
}

Write-Log "Wrote manifest: $manifest"
