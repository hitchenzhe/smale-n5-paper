$ErrorActionPreference = 'Stop'
Push-Location $PSScriptRoot
try {
    New-Item -ItemType Directory -Path 'build' -Force | Out-Null
    $latex = (Get-Command pdflatex -ErrorAction Stop).Source
    for ($pass = 1; $pass -le 3; $pass++) {
        & $latex '-interaction=nonstopmode' '-halt-on-error' '-file-line-error' '-output-directory=build' 'main.tex'
        if ($LASTEXITCODE -ne 0) { throw "LaTeX compilation failed on pass $pass; see build/main.log." }
    }
    Copy-Item -LiteralPath 'build/main.pdf' -Destination 'main.pdf' -Force
    Write-Host 'Created main.pdf'
} finally {
    Pop-Location
}
