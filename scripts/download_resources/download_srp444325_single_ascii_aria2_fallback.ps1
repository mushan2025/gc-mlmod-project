param(
    [string]$ProjectRoot = (Resolve-Path (Join-Path $PSScriptRoot "..\..")).Path,
    [string]$FileName = "SRR24947496_1.fastq.gz",
    [string]$TempRoot = "C:\srp444325_aria2_tmp",
    [int]$MaxTries = 50
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
    Add-Content -LiteralPath (Join-Path $logDir "ascii_aria2_fallback.log") -Value $line -Encoding UTF8
}

function Get-Md5Lower {
    param([string]$Path)
    if (-not (Test-Path -LiteralPath $Path)) { return "" }
    return (Get-FileHash -Algorithm MD5 -LiteralPath $Path).Hash.ToLowerInvariant()
}

function Resolve-Aria2 {
    $cmd = Get-Command aria2c -ErrorAction SilentlyContinue
    if ($cmd) { return $cmd.Source }

    $wingetRoot = Join-Path $env:LOCALAPPDATA "Microsoft\WinGet\Packages"
    if (Test-Path -LiteralPath $wingetRoot) {
        $found = Get-ChildItem -Path $wingetRoot -Recurse -Filter aria2c.exe -ErrorAction SilentlyContinue | Select-Object -First 1
        if ($found) { return $found.FullName }
    }

    throw "aria2c was not found."
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

$planRows = Import-Csv -Delimiter "`t" -LiteralPath $planPath |
    Where-Object { $_.file_name -like "*.fastq.gz" -and $_.include_in_main_download -ne "False" }
$targetRow = $planRows | Where-Object { $_.file_name -eq $FileName } | Select-Object -First 1
if (-not $targetRow) {
    throw "FileName not found in plan: $FileName"
}

$url = $targetRow.source_url
if ($url -match '^ftp://') {
    $url = $url -replace '^ftp://', 'https://'
}

$expectedBytes = [int64]$targetRow.expected_bytes
$expectedMd5 = $targetRow.expected_md5
$finalPath = Join-Path $outDir $FileName
$tmpPath = Join-Path $TempRoot $FileName
$tmpControl = "$tmpPath.aria2"
$aria2Exe = Resolve-Aria2
$aria2Log = Join-Path $logDir ("SRP444325_ascii_aria2_" + (Get-Date -Format "yyyyMMdd_HHmmss") + ".log")

if (Test-Path -LiteralPath $finalPath) {
    $existingSize = (Get-Item -LiteralPath $finalPath).Length
    $existingMd5 = if ($existingSize -eq $expectedBytes) { Get-Md5Lower -Path $finalPath } else { "" }
    if ($existingSize -eq $expectedBytes -and $existingMd5 -eq $expectedMd5) {
        Write-Log "Final file already passes size and MD5 checks."
        $manifestRows = Write-FastqManifest -PlanRows $planRows
        $okCount = ($manifestRows | Where-Object { $_.size_status -eq "size_ok" -and $_.md5_status -eq "md5_ok" }).Count
        Write-Log "Manifest written: $manifestPath; md5_ok=$okCount/$($planRows.Count)"
        exit 0
    }
    Write-Log "Removing invalid final target before ASCII aria2 fallback: $finalPath"
    Remove-Item -LiteralPath $finalPath -Force
}

Write-Log "Starting ASCII aria2 fallback for $FileName"
Write-Log "URL: $url"
Write-Log "Temp path: $tmpPath"

& $aria2Exe `
    --dir=$TempRoot `
    --out=$FileName `
    --continue=true `
    --max-concurrent-downloads=1 `
    --max-connection-per-server=4 `
    --split=4 `
    --min-split-size=1M `
    --file-allocation=none `
    --auto-file-renaming=false `
    --allow-overwrite=true `
    --check-integrity=true `
    --checksum="md5=$expectedMd5" `
    --retry-wait=30 `
    --max-tries=$MaxTries `
    --timeout=60 `
    --summary-interval=60 `
    --download-result=hide `
    --console-log-level=notice `
    --log=$aria2Log `
    --log-level=notice `
    $url

if ($LASTEXITCODE -ne 0) {
    throw "ASCII aria2 fallback failed with exit code $LASTEXITCODE. See $aria2Log"
}

if (-not (Test-Path -LiteralPath $tmpPath)) {
    throw "ASCII aria2 fallback did not create expected file: $tmpPath"
}
if ((Get-Item -LiteralPath $tmpPath).Length -ne $expectedBytes) {
    throw "ASCII aria2 fallback size mismatch: $tmpPath"
}
$actualMd5 = Get-Md5Lower -Path $tmpPath
Write-Log "ASCII fallback MD5: $actualMd5"
if ($actualMd5 -ne $expectedMd5) {
    throw "ASCII aria2 fallback MD5 mismatch."
}

if (Test-Path -LiteralPath $tmpControl) {
    Remove-Item -LiteralPath $tmpControl -Force
}
Move-Item -LiteralPath $tmpPath -Destination $finalPath -Force
Write-Log "Moved verified file into project: $finalPath"

$manifestRows = Write-FastqManifest -PlanRows $planRows
$okRows = $manifestRows | Where-Object { $_.size_status -eq "size_ok" -and $_.md5_status -eq "md5_ok" }
Write-Log "Manifest written: $manifestPath; md5_ok=$($okRows.Count)/$($planRows.Count)"
if ($okRows.Count -ne $planRows.Count) {
    throw "Not all SRP444325 main FASTQ files pass final size+MD5 checks."
}

Write-Log "All SRP444325 main FASTQ files pass final size+MD5 checks."
