param(
    [string]$ProjectRoot = (Resolve-Path (Join-Path $PSScriptRoot "..\..")).Path,
    [int]$MaxAttempts = 80,
    [int]$SleepSeconds = 120,
    [int]$MaxConcurrentDownloads = 3,
    [int]$ConnectionsPerFile = 2
)

$ErrorActionPreference = "Stop"

$dataset = "SRP444325"
$fastqDir = Join-Path $ProjectRoot "data\public_downloads\SRP444325\fastq"
$reportDir = Join-Path $ProjectRoot "reports\download_resources"
$logDir = Join-Path $ProjectRoot "logs\download_resources\SRP444325"
$planPath = Join-Path $reportDir "SRP444325_fastq_download_plan.tsv"
$manifestPath = Join-Path $reportDir "SRP444325_fastq_download_manifest.tsv"
$batchScript = Join-Path $PSScriptRoot "download_srp444325_fastq_aria2_batch.ps1"
$verifyScript = Join-Path $PSScriptRoot "download_srp444325_fastq.ps1"
$supervisorLog = Join-Path $logDir "SRP444325_resume_supervisor.log"

New-Item -ItemType Directory -Force -Path $fastqDir, $reportDir, $logDir | Out-Null

function Write-SupervisorLog {
    param([string]$Message)
    $stamp = Get-Date -Format "yyyy-MM-dd HH:mm:ss"
    $line = "[$stamp] $Message"
    Write-Host $line
    Add-Content -Path $supervisorLog -Value $line -Encoding UTF8
}

function Get-MainFastqRows {
    if (-not (Test-Path $planPath)) {
        throw "Download plan not found: $planPath"
    }

    return @(Import-Csv -Delimiter "`t" $planPath | Where-Object {
        $_.file_name -like "*.fastq.gz" -and $_.include_in_main_download -ne "False"
    })
}

function Get-SizeCompletionSummary {
    $rows = Get-MainFastqRows
    $complete = 0
    $partial = 0
    $missing = 0
    $completeBytes = [int64]0

    foreach ($row in $rows) {
        $file = Join-Path $fastqDir $row.file_name
        $control = "$file.aria2"
        $expected = [int64]$row.expected_bytes

        if (-not (Test-Path $file)) {
            $missing += 1
            continue
        }

        $local = (Get-Item $file).Length
        if (Test-Path $control) {
            $partial += 1
        } elseif ($local -eq $expected) {
            $complete += 1
            $completeBytes += $local
        } else {
            $partial += 1
        }
    }

    return [pscustomobject]@{
        total = $rows.Count
        complete = $complete
        partial = $partial
        missing = $missing
        complete_bytes = $completeBytes
    }
}

function Test-ManifestComplete {
    if (-not (Test-Path $manifestPath)) {
        return $false
    }
    $rows = @(Import-Csv -Delimiter "`t" $manifestPath | Where-Object {
        $_.file_name -like "*.fastq.gz" -and $_.include_in_main_download -ne "False"
    })
    if ($rows.Count -ne 30) {
        return $false
    }
    $bad = @($rows | Where-Object { $_.size_status -ne "size_ok" -or $_.md5_status -ne "md5_ok" })
    return ($bad.Count -eq 0)
}

Write-SupervisorLog "Starting $dataset resume supervisor; max_attempts=$MaxAttempts; sleep_seconds=$SleepSeconds; max_concurrent=$MaxConcurrentDownloads; connections_per_file=$ConnectionsPerFile"

for ($attempt = 1; $attempt -le $MaxAttempts; $attempt++) {
    $summary = Get-SizeCompletionSummary
    Write-SupervisorLog "Attempt $attempt/$MaxAttempts preflight: complete=$($summary.complete)/$($summary.total); partial=$($summary.partial); missing=$($summary.missing); complete_bytes=$($summary.complete_bytes)"

    if (Test-ManifestComplete) {
        Write-SupervisorLog "Manifest already complete with size_ok and md5_ok for all main FASTQ files."
        exit 0
    }

    & powershell -NoProfile -ExecutionPolicy Bypass -File $batchScript `
        -ProjectRoot $ProjectRoot `
        -MaxConcurrentDownloads $MaxConcurrentDownloads `
        -ConnectionsPerFile $ConnectionsPerFile

    $exitCode = $LASTEXITCODE
    Write-SupervisorLog "Attempt $attempt batch exit code: $exitCode"

    if ($exitCode -eq 0) {
        & powershell -NoProfile -ExecutionPolicy Bypass -File $verifyScript `
            -ProjectRoot $ProjectRoot `
            -Aria2Connections $ConnectionsPerFile
        $verifyExit = $LASTEXITCODE
        Write-SupervisorLog "Attempt $attempt verification exit code: $verifyExit"

        if ((Test-ManifestComplete) -and ($verifyExit -eq 0)) {
            $summary = Get-SizeCompletionSummary
            Write-SupervisorLog "Completed ${dataset}: complete=$($summary.complete)/$($summary.total); complete_bytes=$($summary.complete_bytes)"
            exit 0
        }
    }

    if ($attempt -lt $MaxAttempts) {
        Write-SupervisorLog "Not complete yet; sleeping $SleepSeconds seconds before retry."
        Start-Sleep -Seconds $SleepSeconds
    }
}

$finalSummary = Get-SizeCompletionSummary
Write-SupervisorLog "Stopped after MaxAttempts without complete manifest: complete=$($finalSummary.complete)/$($finalSummary.total); partial=$($finalSummary.partial); missing=$($finalSummary.missing)"
exit 1
