<#
.SYNOPSIS
Generate a simple CHANGELOG.md from git history (simple, local-only)
#>
param(
    [string]$OutFile = "CHANGELOG.md",
    [int]$MaxEntries = 200
)

if (-not (Get-Command git -ErrorAction SilentlyContinue)) {
    Write-Error "git is not installed or not in PATH"
    exit 1
}

# Get latest commits
$commits = git --no-pager log -n $MaxEntries --pretty=format:"%h|%ad|%an|%s" --date=short
if (-not $commits) {
    Write-Host "No commits found."
    exit 0
}

$lines = @()
$lines += "# Changelog"
$lines += ""
$lines += "Generated from git log on $(Get-Date -Format o)"
$lines += ""

foreach ($c in $commits -split "`n") {
    $parts = $c -split "\|",4
    if ($parts.Length -ge 4) {
        $hash = $parts[0]
        $date = $parts[1]
        $author = $parts[2]
        $msg = $parts[3]
        $lines += "- [$date] $msg ($author) - $hash"
    }
}

$lines | Out-File -FilePath $OutFile -Encoding UTF8
Write-Host "Wrote $OutFile with $($lines.Count - 4) entries."
