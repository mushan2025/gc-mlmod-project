param(
    [string]$ProjectRoot = (Resolve-Path (Join-Path $PSScriptRoot "..\..")).Path,
    [string]$FileName = "SRR24947496_1.fastq.gz",
    [string]$TempRoot = (Join-Path ([System.IO.Path]::GetTempPath()) "srp444325_curl_tmp"),
    [int]$MaxAttempts = 3
)

$ErrorActionPreference = "Stop"

$dataset = "SRP444325"
$planPath = Join-Path $ProjectRoot "reports\download_resources\SRP444325_fastq_download_plan.tsv"
$outDir = Join-Path $ProjectRoot "data\public_downloads\SRP444325\fastq"
$logDir = Join-Path $ProjectRoot "logs\download_resources\SRP444325"
$reportDir = Join-Path $ProjectRoot "reports\download_resources"
$manifestPath = Join-Path $reportDir "SRP444325_fastq_download_manifest.tsv"

New-Item -ItemType Directory -Force -Path $outDir, $logDir, $reportDir, $TempRoot | Out-Null

function Write-Log {
    param([string]$Message)
    $stamp = Get-Date -Format "yyyy-MM-dd HH:mm:ss"
    $line = "[$stamp] $Message"
    Write-Host $line
    Add-Content -LiteralPath (Join-Path $logDir "curl_fallback.log") -Value $line -Encoding UTF8
}

function Get-Md5Lower {
    param([string]$Path)
    if (-not (Test-Path -LiteralPath $Path)) { return "" }
    return (Get-FileHash -Algorithm MD5 -LiteralPath $Path).Hash.ToLowerInvariant()
}

function Write-FastqManifest {
    param([object[]]$PlanRows)

    $rows = foreach ($row in $PlanRows) {
        $path = Join-Path $outDir $row.file_name
        $exists = Test-Path -LiteralPath $path
        $bytes = if ($exists) { (Get-Item -LiteralPath $path).Length } else { 0 }
        $expectedBytes = [int64]$row.expected_bytes
        $sizeStatus = if (-not $exists) {
            "missing"
        } elseif ($bytes -eq $expectedBytes) {
            "size_ok"
        } else {
            "size_mismatch"
        }
        $observedMd5 = if ($exists -and $sizeStatus -eq "size_ok") { Get-Md5Lower -Path $path } else { "" }
        $md5Status = if ([string]::IsNullOrWhiteSpace($row.expected_md5)) {
            "no_expected_md5"
        } elseif ($observedMd5 -eq "") {
            "not_checked"
        } elseif ($observedMd5 -eq $row.expected_md5) {
            "md5_ok"
        } else {
            "md5_mismatch"
        }

        [pscustomobject]@{
            dataset_id = $dataset
            run_accession = $row.run_accession
            sample_accession = $row.sample_accession
            experiment_accession = $row.experiment_accession
            file_name = $row.file_name
            source_url = $row.source_url
            local_path = $row.local_path
            expected_bytes = $row.expected_bytes
            local_size_bytes = $bytes
            size_status = $sizeStatus
            expected_md5 = $row.expected_md5
            observed_md5 = $observedMd5
            md5_status = $md5Status
            include_in_main_download = $row.include_in_main_download
            checked_at = (Get-Date -Format "yyyy-MM-dd HH:mm:ss")
        }
    }

    $rows | Export-Csv -LiteralPath $manifestPath -Delimiter "`t" -NoTypeInformation -Encoding UTF8
    return $rows
}

if (-not (Test-Path -LiteralPath $planPath)) {
    throw "Plan file not found: $planPath"
}

$planRows = Import-Csv -Delimiter "`t" -LiteralPath $planPath |
    Where-Object { $_.file_name -like "*.fastq.gz" -and $_.include_in_main_download -ne "False" }
$targetRow = $planRows | Where-Object { $_.file_name -eq $FileName } | Select-Object -First 1
if (-not $targetRow) {
    throw "FileName not found in plan: $FileName"
}

$url = $targetRow.source_url
$expectedBytes = [int64]$targetRow.expected_bytes
$expectedMd5 = $targetRow.expected_md5
$finalPath = Join-Path $outDir $FileName
$aria2Control = "$finalPath.aria2"
$tempPart = Join-Path $TempRoot "$FileName.part"
$tempDone = Join-Path $TempRoot $FileName

Write-Log "Starting curl fallback for $FileName"
Write-Log "URL: $url"
Write-Log "Temp part: $tempPart"

if (Test-Path -LiteralPath $finalPath) {
    $existingSize = (Get-Item -LiteralPath $finalPath).Length
    $existingMd5 = if ($existingSize -eq $expectedBytes) { Get-Md5Lower -Path $finalPath } else { "" }
    if ($existingSize -eq $expectedBytes -and $existingMd5 -eq $expectedMd5 -and -not (Test-Path -LiteralPath $aria2Control)) {
        Write-Log "Final file already passes size and MD5 checks."
        $manifestRows = Write-FastqManifest -PlanRows $planRows
        $okCount = ($manifestRows | Where-Object { $_.size_status -eq "size_ok" -and $_.md5_status -eq "md5_ok" }).Count
        Write-Log "Manifest written: $manifestPath; md5_ok=$okCount/$($planRows.Count)"
        exit 0
    }
    Write-Log "Removing invalid existing final file before fallback redownload."
    Remove-Item -LiteralPath $finalPath -Force
}
if (Test-Path -LiteralPath $aria2Control) {
    Remove-Item -LiteralPath $aria2Control -Force
}

$curl = Get-Command curl.exe -ErrorAction Stop
$downloadOk = $false
for ($attempt = 1; $attempt -le $MaxAttempts; $attempt++) {
    Write-Log "curl attempt $attempt/$MaxAttempts"
    & $curl.Source `
        --location `
        --fail `
        --ftp-pasv `
        --retry 30 `
        --retry-all-errors `
        --retry-delay 10 `
        --connect-timeout 60 `
        --speed-time 120 `
        --speed-limit 1024 `
        --continue-at - `
        --output $tempPart `
        $url

    $exit = $LASTEXITCODE
    $size = if (Test-Path -LiteralPath $tempPart) { (Get-Item -LiteralPath $tempPart).Length } else { 0 }
    Write-Log "curl attempt $attempt exit=$exit size=$size expected=$expectedBytes"

    if ($exit -eq 0 -and $size -eq $expectedBytes) {
        $actualMd5 = Get-Md5Lower -Path $tempPart
        Write-Log "Downloaded MD5: $actualMd5"
        if ($actualMd5 -eq $expectedMd5) {
            $downloadOk = $true
            break
        }
        Write-Log "MD5 mismatch in temp file; deleting temp part before next attempt."
        Remove-Item -LiteralPath $tempPart -Force
    } elseif ($exit -ne 0) {
        Write-Log "curl returned non-zero exit code; keeping temp part for resume."
    } else {
        Write-Log "Size mismatch after curl exit 0; keeping temp part for resume."
    }
}

if (-not $downloadOk) {
    throw "curl fallback failed for $FileName after $MaxAttempts attempts"
}

Move-Item -LiteralPath $tempPart -Destination $tempDone -Force
Move-Item -LiteralPath $tempDone -Destination $finalPath -Force
Write-Log "Moved verified file into project: $finalPath"

$manifestRows = Write-FastqManifest -PlanRows $planRows
$okRows = $manifestRows | Where-Object { $_.size_status -eq "size_ok" -and $_.md5_status -eq "md5_ok" }
Write-Log "Manifest written: $manifestPath; md5_ok=$($okRows.Count)/$($planRows.Count)"

if ($okRows.Count -ne $planRows.Count) {
    throw "Not all SRP444325 main FASTQ files pass final size+MD5 checks."
}

Write-Log "All SRP444325 main FASTQ files pass final size+MD5 checks."
