param(
    [string]$ProjectRoot = (Resolve-Path (Join-Path $PSScriptRoot "..\..")).Path,
    [switch]$Force,
    [switch]$SkipGzipTest,
    [switch]$PlanOnly
)

$ErrorActionPreference = "Stop"

$dataset = "GSE239676"
$baseUrl = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE239nnn/GSE239676/suppl"
$outDir = Join-Path $ProjectRoot "data\public_downloads\GSE239676"
$logDir = Join-Path $ProjectRoot "logs\download_resources\GSE239676"
$reportDir = Join-Path $ProjectRoot "reports\download_resources"
$manifest = if ($PlanOnly) {
    Join-Path $reportDir "GSE239676_download_plan.tsv"
} else {
    Join-Path $reportDir "GSE239676_download_manifest.tsv"
}

New-Item -ItemType Directory -Force -Path $outDir, $logDir, $reportDir | Out-Null

$files = @(
    "GSE239676_barcodes.tsv.gz",
    "GSE239676_features.tsv.gz",
    "GSE239676_meta.tsv.gz",
    "GSE239676_count_matrix.mtx.gz"
)

function Write-Log {
    param([string]$Message)
    $stamp = Get-Date -Format "yyyy-MM-dd HH:mm:ss"
    $line = "[$stamp] $Message"
    Write-Host $line
    Add-Content -Path (Join-Path $logDir "download.log") -Value $line -Encoding UTF8
}

function Get-RemoteSize {
    param([string]$Url)
    $head = & curl.exe -L --ssl-no-revoke --silent --show-error --head $Url
    $contentLength = ($head | Select-String -Pattern '^Content-Length:\s*(\d+)' | Select-Object -Last 1)
    if ($contentLength) {
        return [int64]$contentLength.Matches[0].Groups[1].Value
    }
    return $null
}

function Test-GzipFile {
    param([string]$Path)
    if ($SkipGzipTest) {
        return "skipped"
    }
    $fs = [System.IO.File]::OpenRead($Path)
    try {
        $b1 = $fs.ReadByte()
        $b2 = $fs.ReadByte()
    } finally {
        $fs.Close()
    }

    if (($b1 -eq 0x1f) -and ($b2 -eq 0x8b)) {
    $code = @"
import gzip
import sys
path = sys.argv[1]
with gzip.open(path, "rb") as fh:
    while fh.read(1024 * 1024):
        pass
"@
    $tmp = Join-Path $logDir "gzip_test.py"
    Set-Content -Path $tmp -Value $code -Encoding UTF8
    & python $tmp $Path
    if ($LASTEXITCODE -eq 0) { return "gzip_read_ok" }
        return "gzip_read_failed"
    }

    $textCode = @"
import sys
path = sys.argv[1]
with open(path, "rt", encoding="utf-8", errors="strict") as fh:
    for _ in range(5):
        if fh.readline() == "":
            break
"@
    $textTmp = Join-Path $logDir "plain_text_test.py"
    Set-Content -Path $textTmp -Value $textCode -Encoding UTF8
    & python $textTmp $Path
    if ($LASTEXITCODE -eq 0) { return "plain_text_read_ok_gz_suffix" }
    return "gzip_and_plain_text_read_failed"
}

$rows = New-Object System.Collections.Generic.List[object]

foreach ($file in $files) {
    $url = "$baseUrl/$file"
    $out = Join-Path $outDir $file
    $fileLog = Join-Path $logDir "$file.curl.log"
    $remoteSize = Get-RemoteSize $url
    $existingSize = if (Test-Path $out) { (Get-Item $out).Length } else { 0 }

    if ($PlanOnly) {
        Write-Log "PLAN $file; existing=$existingSize remote=$remoteSize"
    } elseif ((-not $Force) -and $remoteSize -and (Test-Path $out) -and ($existingSize -eq $remoteSize)) {
        Write-Log "SKIP complete file: $file ($existingSize bytes)"
    } else {
        Write-Log "DOWNLOAD/RESUME $file; existing=$existingSize remote=$remoteSize"
        & curl.exe `
            --location `
            --fail `
            --silent `
            --show-error `
            --ssl-no-revoke `
            --retry 10 `
            --retry-delay 30 `
            --connect-timeout 60 `
            --speed-time 180 `
            --speed-limit 1024 `
            --continue-at - `
            --output $out `
            $url 2>&1 | Tee-Object -FilePath $fileLog -Append
        if ($LASTEXITCODE -ne 0) {
            throw "curl failed for $file with exit code $LASTEXITCODE"
        }
    }

    $localSize = if (Test-Path $out) { (Get-Item $out).Length } else { 0 }
    $sizeStatus = if ($remoteSize -and ($localSize -eq $remoteSize)) { "size_ok" } elseif ($remoteSize) { "size_mismatch" } else { "remote_size_unknown" }
    $sha256 = if (Test-Path $out) { (Get-FileHash -Algorithm SHA256 $out).Hash } else { "" }
    $gzipStatus = if ($PlanOnly) { "not_checked_plan_only" } elseif (Test-Path $out) { Test-GzipFile $out } else { "missing" }

    $rows.Add([pscustomobject]@{
        dataset_id = $dataset
        file_name = $file
        source_url = $url
        local_path = $out.Replace($ProjectRoot + "\", "")
        remote_size_bytes = $remoteSize
        local_size_bytes = $localSize
        size_status = $sizeStatus
        sha256 = $sha256
        gzip_status = $gzipStatus
        checked_at = (Get-Date -Format "yyyy-MM-dd HH:mm:ss")
    })
}

$rows | Export-Csv -Path $manifest -Delimiter "`t" -NoTypeInformation -Encoding UTF8
Write-Log "Wrote manifest: $manifest"
