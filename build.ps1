$ErrorActionPreference = 'Stop'
Push-Location $PSScriptRoot
$savedTexInputs = $env:TEXINPUTS
$savedTfmFonts = $env:TFMFONTS
$savedT1Fonts = $env:T1FONTS
try {
    $source = 'Smale_5-CT-v19.tex'
    $buildDir = Join-Path $PSScriptRoot 'build'
    $depsDir = Join-Path $buildDir 'tex-deps'
    New-Item -ItemType Directory -Path $buildDir -Force | Out-Null
    if (Test-Path -LiteralPath $depsDir) {
        $env:TEXINPUTS = ($depsDir.Replace('\','/') + '//;')
        $env:TFMFONTS = $env:TEXINPUTS
        $env:T1FONTS = $env:TEXINPUTS
        Get-ChildItem -LiteralPath (Join-Path $depsDir 'generated-fonts') -Filter '*.pk' -File | Copy-Item -Destination $buildDir -Force
    }
    Set-Location -LiteralPath $buildDir
    $latex = (Get-Command pdflatex -ErrorAction Stop).Source
    for ($pass = 1; $pass -le 3; $pass++) {
        & $latex '--disable-installer' '-no-shell-escape' '-interaction=nonstopmode' '-halt-on-error' '-file-line-error' "-include-directory=$PSScriptRoot" "-output-directory=$buildDir" (Join-Path $PSScriptRoot $source)
        if ($LASTEXITCODE -ne 0) { throw "LaTeX pass $pass failed; see build/Smale_5-CT-v19.log. Automatic installation is disabled." }
    }
    Copy-Item -LiteralPath (Join-Path $buildDir 'Smale_5-CT-v19.pdf') -Destination (Join-Path $PSScriptRoot 'Smale_5-CT-v19.pdf') -Force
    Write-Host 'Created Smale_5-CT-v19.pdf without automatic package installation.'
} finally {
    $env:TEXINPUTS = $savedTexInputs
    $env:TFMFONTS = $savedTfmFonts
    $env:T1FONTS = $savedT1Fonts
    Pop-Location
}
