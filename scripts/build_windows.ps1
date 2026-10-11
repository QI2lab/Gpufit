param([switch]$Initialized, [string]$CudaRoot, [string]$Architectures, [string]$ProjectPath,
      [ValidateSet('build', 'sync', 'rebuild')][string]$Action = 'build')
$ErrorActionPreference = 'Stop'
$projectRoot = Split-Path -Parent $PSScriptRoot
if (-not $Initialized) {
    $vswhere = Join-Path ${env:ProgramFiles(x86)} 'Microsoft Visual Studio\Installer\vswhere.exe'
    if (-not (Test-Path -LiteralPath $vswhere)) { throw 'Install VS 2022 standalone Build Tools and MSVC 14.38.' }
    $vsRoot = & $vswhere -latest -products '*' -requires Microsoft.VisualStudio.Component.VC.14.38.17.8.x86.x64 -property installationPath
    if (-not $vsRoot) { throw 'MSVC 14.38 x64/x86 component is missing.' }
    $vcvars = Join-Path $vsRoot 'VC\Auxiliary\Build\vcvarsall.bat'
    $env:GPUFIT_REQUESTED_CUDA = $CudaRoot
    $env:GPUFIT_REQUESTED_ARCHITECTURES = $Architectures
    $env:GPUFIT_REQUESTED_ACTION = $Action
    $env:GPUFIT_REQUESTED_PROJECT = $ProjectPath
    & cmd /d /c "call `"$vcvars`" x64 -vcvars_ver=14.38 >nul && powershell -NoProfile -ExecutionPolicy Bypass -File `"$PSCommandPath`" -Initialized"
    exit $LASTEXITCODE
}
if (-not $env:VCToolsVersion.StartsWith('14.38.')) { throw 'The initialized compiler must be MSVC 14.38.' }
if ($env:GPUFIT_REQUESTED_ACTION) { $Action = $env:GPUFIT_REQUESTED_ACTION }
if (-not $CudaRoot) { $CudaRoot = $env:GPUFIT_REQUESTED_CUDA }
if (-not $CudaRoot) { $CudaRoot = Join-Path $env:ProgramFiles 'NVIDIA GPU Computing Toolkit\CUDA\v12.8' }
if (-not (Test-Path -LiteralPath (Join-Path $CudaRoot 'bin\nvcc.exe'))) { throw 'Install CUDA 12.8 Update 1 nvcc and cuBLAS development components.' }
$env:CUDA_PATH = $CudaRoot
$env:CUDAToolkit_ROOT = $CudaRoot
$env:CUDACXX = Join-Path $CudaRoot 'bin\nvcc.exe'
$vsCMakeRoot = Join-Path $env:VSINSTALLDIR 'Common7\IDE\CommonExtensions\Microsoft\CMake'
$env:PATH = "$(Join-Path $vsCMakeRoot 'CMake\bin');$(Join-Path $vsCMakeRoot 'Ninja');$(Join-Path $CudaRoot 'bin');$env:PATH"
$env:UV_CACHE_DIR = Join-Path $projectRoot '.uv-cache'
$policy = (Get-Content -LiteralPath (Join-Path $projectRoot 'CMakePresets.json') -Raw | ConvertFrom-Json).configurePresets[0]
if (-not $ProjectPath) { $ProjectPath = $env:GPUFIT_REQUESTED_PROJECT }
if (-not $ProjectPath) { $ProjectPath = $projectRoot }
$ProjectPath = (Resolve-Path -LiteralPath $ProjectPath).Path
if ($Action -eq 'build' -and $ProjectPath -ne $projectRoot) { throw 'External integration projects support sync/rebuild only.' }
Set-Location -LiteralPath $ProjectPath
$settingPrefix = '--config-setting='
if ($ProjectPath -ne $projectRoot) { $settingPrefix = '--config-settings-package=pygpufit:' }
$buildSettings = @("${settingPrefix}build-dir=build/microscopy")
foreach ($setting in $policy.cacheVariables.PSObject.Properties) {
    $value = $setting.Value
    if ($setting.Name -eq 'CMAKE_CUDA_ARCHITECTURES' -and $env:GPUFIT_REQUESTED_ARCHITECTURES) { $value = $env:GPUFIT_REQUESTED_ARCHITECTURES }
    $buildSettings += "${settingPrefix}cmake.define.$($setting.Name)=$value"
}
$env:CMAKE_GENERATOR = $policy.generator
& uv --version
& cmake --version
& ninja --version
& $env:CUDACXX --version
& uv lock --check
if ($LASTEXITCODE) { exit $LASTEXITCODE }
New-Item -ItemType Directory -Force -Path (Join-Path $projectRoot 'build') | Out-Null
$ErrorActionPreference = 'Continue' # Windows PowerShell wraps native stderr, including harmless warnings.
if ($Action -in @('sync', 'rebuild')) {
    $syncOptions = @()
    if ($Action -eq 'rebuild') { $syncOptions += '--reinstall-package=pygpufit' }
    & uv sync --locked @syncOptions @buildSettings
    exit $LASTEXITCODE
}
& uv build --wheel --sdist @buildSettings 2>&1 | Tee-Object -FilePath build/build_windows.log
$buildExit = $LASTEXITCODE
exit $buildExit
