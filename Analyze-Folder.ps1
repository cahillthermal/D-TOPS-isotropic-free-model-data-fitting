<#
.SYNOPSIS
    Analyzes all .txt data files in a specified folder using main.py.

.DESCRIPTION
    This PowerShell script iterates through all .txt data files in a specified folder,
    invokes main.py for each file, and records fitting results in results.dat.

.PARAMETER FolderPath
    Path to the folder containing .txt data files. Defaults to current directory (".").

.PARAMETER PythonExe
    Path to Python executable. Defaults to workspace .venv if found, or 'python'.

.PARAMETER Plot
    Switch to enable plotting during analysis (plotting is disabled by default for batch runs).

.EXAMPLE
    .\Analyze-Folder.ps1 -FolderPath ".\data"
#>

[CmdletBinding()]
param (
    [Parameter(Position = 0, Mandatory = $false)]
    [string]$FolderPath = ".",

    [Parameter(Mandatory = $false)]
    [string]$PythonExe = "",

    [Parameter(Mandatory = $false)]
    [switch]$Plot
)

# Determine Python executable
if ([string]::IsNullOrWhitespace($PythonExe)) {
    $venvPython = Join-Path $PSScriptRoot ".venv\Scripts\python.exe"
    if (Test-Path $venvPython) {
        $PythonExe = $venvPython
    } else {
        $PythonExe = "python"
    }
}

# Validate folder path
if (-not (Test-Path -Path $FolderPath)) {
    Write-Error "The specified folder path does not exist: $FolderPath"
    exit 1
}

$resolvedFolder = (Get-Item -Path $FolderPath).FullName
Write-Host "Searching for .txt data files in: $resolvedFolder" -ForegroundColor Cyan

$dataFiles = Get-ChildItem -Path $resolvedFolder -Filter "*.txt" -File | Where-Object { $_.Name -notmatch '^requirements.*\.txt$' }

if ($dataFiles.Count -eq 0) {
    Write-Host "No .txt files found in $resolvedFolder" -ForegroundColor Yellow
    exit 0
}

Write-Host "Found $($dataFiles.Count) .txt file(s) to analyze." -ForegroundColor Green

$mainPyPath = Join-Path $PSScriptRoot "main.py"
if (-not (Test-Path $mainPyPath)) {
    $mainPyPath = "main.py"
}

$successful = 0
$failed = 0

foreach ($file in $dataFiles) {
    Write-Host "`n--------------------------------------------------" -ForegroundColor Yellow
    Write-Host "Processing file: $($file.Name)" -ForegroundColor Cyan
    Write-Host "--------------------------------------------------"

    $cmdArgs = @("$mainPyPath", "--file", "$($file.FullName)")
    if (-not $Plot) {
        $cmdArgs += "--no-plot"
    }

    try {
        & $PythonExe @cmdArgs
        if ($LASTEXITCODE -eq 0) {
            $successful++
            Write-Host "Successfully analyzed $($file.Name)" -ForegroundColor Green
        } else {
            $failed++
            Write-Warning "Failed analyzing $($file.Name) (Exit Code: $LASTEXITCODE)"
        }
    } catch {
        $failed++
        Write-Error "Exception encountered while processing $($file.Name): $_"
    }
}

Write-Host "`n==================================================" -ForegroundColor Cyan
Write-Host "Batch Analysis Summary" -ForegroundColor Green
Write-Host "Total files: $($dataFiles.Count)"
Write-Host "Successful:  $successful" -ForegroundColor Green
if ($failed -gt 0) {
    Write-Host "Failed:      $failed" -ForegroundColor Red
}
Write-Host "Fitted results appended to results.dat" -ForegroundColor Cyan
Write-Host "=================================================="
