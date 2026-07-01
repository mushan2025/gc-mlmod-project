param(
    [string]$ProjectRoot = (Resolve-Path (Join-Path $PSScriptRoot "..\..")).Path,
    [string]$EnaRunTable = "",
    [string]$Aria2Path = "",
    [int]$MaxConcurrentDownloads = 4,
    [int]$ConnectionsPerFile = 4,
    [switch]$IncludeNonPairedArtifacts
)

$ErrorActionPreference = "Stop"

$dataset = "SRP444325"
$outDir = Join-Path $ProjectRoot "data\public_downloads\SRP444325\fastq"
$metaDir = Join-Path $ProjectRoot "data\public_downloads\SRP444325"
$logDir = Join-Path $ProjectRoot "logs\download_resources\SRP444325"
$reportDir = Join-Path $ProjectRoot "reports\download_resources"
$inputFile = Join-Path $logDir "SRP444325_aria2_input.txt"
$preflightManifest = Join-Path $reportDir "SRP444325_aria2_preflight.tsv"

New-Item -ItemType Directory -Force -Path $outDir, $metaDir, $logDir, $reportDir | Out-Null

if ([string]::IsNullOrWhiteSpace($EnaRunTable)) {
    $EnaRunTable = Join-Path $metaDir "SRP444325_ENA_read_run.tsv"
}
if (-not (Test-Path $EnaRunTable)) {
    throw "ENA run table not found: $EnaRunTable"
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
    if ($FtpPath -match '^ftp://') { return $FtpPath }
    if ($FtpPath -match '^ftp\.sra\.ebi\.ac\.uk/') { return "ftp://$FtpPath" }
    return "https://$FtpPath"
}

function Resolve-Aria2 {
    if (-not [string]::IsNullOrWhiteSpace($Aria2Path)) {
        if (Test-Path $Aria2Path) { return (Resolve-Path $Aria2Path).Path }
        throw "Aria2Path was provided but not found: $Aria2Path"
    }

    $cmd = Get-Command aria2c -ErrorAction SilentlyContinue
    if ($cmd) { return $cmd.Source }

    $wingetRoot = Join-Path $env:LOCALAPPDATA "Microsoft\WinGet\Packages"
    if (Test-Path $wingetRoot) {
        $found = Get-ChildItem -Path $wingetRoot -Recurse -Filter aria2c.exe -ErrorAction SilentlyContinue | Select-Object -First 1
        if ($found) { return $found.FullName }
    }

    throw "aria2c was not found. Install aria2 or pass -Aria2Path."
}

$aria2Exe = Resolve-Aria2
$runs = Import-Csv -Path $EnaRunTable -Delimiter "`t"
$jobs = New-Object System.Collections.Generic.List[object]

foreach ($run in $runs) {
    $urls = @($run.fastq_ftp -split ';' | Where-Object { -not [string]::IsNullOrWhiteSpace($_) })
    $md5s = @($run.fastq_md5 -split ';' | Where-Object { -not [string]::IsNullOrWhiteSpace($_) })
    $bytes = @($run.fastq_bytes -split ';' | Where-Object { $_ -match '^\d+$' })

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
            file_name = $name
            source_url = $url
            expected_md5 = if ($i -lt $md5s.Count) { $md5s[$i] } else { "" }
            expected_bytes = if ($i -lt $bytes.Count) { [int64]$bytes[$i] } else { $null }
            include_in_main_download = $includeInMainDownload
            exclusion_reason = $exclusionReason
        })
    }
}

$inputLines = New-Object System.Collections.Generic.List[string]
$rows = New-Object System.Collections.Generic.List[object]
$outDirForAria = $outDir.Replace("\", "/")
$downloadCount = 0

foreach ($job in $jobs) {
    $fileName = [string]$job.file_name
    $out = Join-Path $outDir $fileName
    $aria2Control = "$out.aria2"
    $localSize = if (Test-Path $out) { (Get-Item $out).Length } else { 0 }
    $skipByPolicy = (-not $IncludeNonPairedArtifacts) -and (-not $job.include_in_main_download)
    $complete = ($job.expected_bytes -and (Test-Path $out) -and ($localSize -eq $job.expected_bytes) -and (-not (Test-Path $aria2Control)))

    $action = if ($skipByPolicy) {
        "excluded"
    } elseif ($complete) {
        "already_complete"
    } else {
        "download_or_resume"
    }

    if ($action -eq "download_or_resume") {
        $downloadCount += 1
        $inputLines.Add([string]$job.source_url)
        $inputLines.Add("  dir=$outDirForAria")
        $inputLines.Add("  out=$fileName")
        if (-not [string]::IsNullOrWhiteSpace($job.expected_md5)) {
            $inputLines.Add("  checksum=md5=$($job.expected_md5)")
        }
    }

    $rows.Add([pscustomobject]@{
        dataset_id = $dataset
        run_accession = $job.run_accession
        file_name = $fileName
        source_url = $job.source_url
        expected_bytes = $job.expected_bytes
        local_size_bytes = $localSize
        include_in_main_download = $job.include_in_main_download
        exclusion_reason = $job.exclusion_reason
        action = $action
        aria2_control_exists = (Test-Path $aria2Control)
        checked_at = (Get-Date -Format "yyyy-MM-dd HH:mm:ss")
    })
}

$inputLines | Set-Content -Path $inputFile -Encoding UTF8
$rows | Export-Csv -Path $preflightManifest -Delimiter "`t" -NoTypeInformation -Encoding UTF8
Write-Log "Wrote aria2 input: $inputFile"
Write-Log "Wrote aria2 preflight manifest: $preflightManifest"
Write-Log "aria2 jobs to download/resume: $downloadCount; max_concurrent=$MaxConcurrentDownloads; connections_per_file=$ConnectionsPerFile"

if ($downloadCount -eq 0) {
    Write-Log "No SRP444325 FASTQ files need download/resume."
} else {
    $aria2Log = Join-Path $logDir ("SRP444325_aria2_batch_" + (Get-Date -Format "yyyyMMdd_HHmmss") + ".log")
    & $aria2Exe `
        --input-file=$inputFile `
        --continue=true `
        --max-concurrent-downloads=$MaxConcurrentDownloads `
        --max-connection-per-server=$ConnectionsPerFile `
        --split=$ConnectionsPerFile `
        --min-split-size=1M `
        --file-allocation=none `
        --auto-file-renaming=false `
        --allow-overwrite=true `
        --check-integrity=true `
        --summary-interval=60 `
        --download-result=hide `
        --console-log-level=notice `
        --retry-wait=60 `
        --max-tries=20 `
        --timeout=60 `
        --log=$aria2Log `
        --log-level=notice

    if ($LASTEXITCODE -ne 0) {
        throw "aria2 batch download failed with exit code $LASTEXITCODE. See $aria2Log"
    }
    Write-Log "aria2 batch download completed: $aria2Log"
}

$verifyScript = Join-Path $PSScriptRoot "download_srp444325_fastq.ps1"
& powershell -NoProfile -ExecutionPolicy Bypass -File $verifyScript -ProjectRoot $ProjectRoot -Aria2Connections $ConnectionsPerFile
if ($LASTEXITCODE -ne 0) {
    throw "Final SRP444325 verification script failed with exit code $LASTEXITCODE"
}
