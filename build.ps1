$ErrorActionPreference = 'Stop'
Push-Location $PSScriptRoot
try {
    $output = Join-Path $PSScriptRoot 'verification\build-v10'
    New-Item -ItemType Directory -Path $output -Force | Out-Null
    $latex = (Get-Command pdflatex -ErrorAction Stop).Source
    for ($pass = 1; $pass -le 3; $pass++) {
        & $latex '-no-shell-escape' '-interaction=nonstopmode' '-halt-on-error' '-file-line-error' "-output-directory=$output" 'Smale_5-CT-v10.tex'
        if ($LASTEXITCODE -ne 0) { throw "LaTeX compilation failed on pass $pass; see verification/build-v10/Smale_5-CT-v10.log." }
    }
    Copy-Item -LiteralPath (Join-Path $output 'Smale_5-CT-v10.pdf') -Destination (Join-Path $PSScriptRoot 'Smale_5-CT-v10.pdf') -Force
    Write-Host 'Created Smale_5-CT-v10.pdf'
} finally {
    Pop-Location
}
