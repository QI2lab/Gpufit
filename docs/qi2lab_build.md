# QI2lab native Python builds

Root `pyproject.toml` is authoritative for `pyGpufit` (import `pygpufit`).
The tested development environment uses Python 3.12.10, uv 0.11.30,
scikit-build-core 0.11.6, MSVC 14.38.33130 (compiler 19.38.33145), Windows SDK
10.0.26100.0, CUDA 12.8 Update 1 (nvcc 12.8.93), CMake 3.31.6-msvc6 and
Ninja 1.12.1. The test device is an RTX 3090 (8.6), driver 610.88.

Install VS 2022 standalone Build Tools with C++ tools and the pinned component
`Microsoft.VisualStudio.Component.VC.14.38.17.8.x86.x64`. Keep an existing SDK and
CMake/Ninja installation. Install CUDA 12.8 Update 1's compiler, runtime, cuBLAS
development/runtime libraries, Thrust and Visual Studio integration. These are
administrator installations with substantial downloads; a restart may be requested.
Build scripts never install system tools, replace drivers or restart the machine.

Official references: [Microsoft components](https://learn.microsoft.com/en-us/visualstudio/install/workload-component-id-vs-build-tools?view=vs-2022),
[installer commands](https://learn.microsoft.com/en-us/visualstudio/install/use-command-line-parameters-to-install-visual-studio?view=vs-2022),
[CUDA archive](https://developer.nvidia.com/cuda-12-8-1-download-archive),
[CUDA Windows guide](https://docs.nvidia.com/cuda/archive/12.8.1/cuda-installation-guide-microsoft-windows/index.html).

From the repository root:

```powershell
uv lock
powershell -NoProfile -File scripts/build_windows.ps1 -Action sync -Architectures 86
uv run --locked --no-sync pytest
powershell -NoProfile -File scripts/build_windows.ps1 -Action build -Architectures 86
```

`sync` installs into this project's `.venv`. `build` produces artifacts in `dist`
without installation. The script initializes the selected compiler and locates
CUDA/CMake/Ninja in the same process that invokes uv. Global PATH is unchanged.
`CMakePresets.json` supplies the shared policy; architecture 86 is the tested
workstation selection. All microscopy builds use cuBLAS and float32.

After CUDA/C++ edits, explicitly rebuild the editable native binary:

```powershell
powershell -NoProfile -File scripts/build_windows.ps1 -Action rebuild -Architectures 86
```

Python edits are editable. Repeating `sync` preserves the fork but does not rebuild
native edits. An already initialized equivalent environment can invoke `uv build`
directly. Normal wheel installation remains standards-compliant without uv.

Wheels contain the native library, cuBLAS/cuBLASLt shared runtime dependencies,
the CUDA license and `build_manifest.json`; CUDA and MSVC runtimes are linked
statically. Tags are `py3-none-win_amd64` for the ctypes ABI. `get_capabilities()`
validates the actual binary's float32 ABI and models 15--18; stock wheels are
incompatible. Dirty source builds are identified in the manifest. The sdist
configure step includes a generated `cmake/source_identity.json` so its rebuild
retains the originating checkout revision and dirty status. This records
provenance; it is not a cryptographic attestation of an untrusted archive.

`fit_constrained` accepts C-contiguous float32 constraints `(nfits, 2*nparams)`,
interleaved lower/upper, and shared int32 `constraint_types` from
`ConstraintType.FREE/LOWER/UPPER/LOWER_UPPER`. Group differing patterns. Each fit
has its own endpoints. Active endpoints must be finite and starts feasible. Equal bounds require fixed
parameters. LSE weights are inverse variances; MLE is a pure-Poisson objective.
Inputs are host buffers with transfers, not zero-copy CuPy.

Upstream model IDs 0--14 are unchanged, including `SPLINE_4D=13` and
`SPLINE_5D=14`. The microscopy IDs were reassigned to retain upstream numbering:

| Model | Previous fork ID | Current ID | Parameters |
| --- | --- | --- | --- |
| `GAUSS_3D_ARB` | 13 | 15 | 7 |
| `GAUSS_LOR_3D_ARB` | 14 | 16 | 7 |
| `GAUSS_3D_ROT_ARB` | 15 | 17 | 11 |
| `GAUSS_3D_ASYM_ARB` | 16 | 18 | 8 |

Use `ModelID` names instead of hardcoded numbers, and rebuild the native library
after updating. Capability ABI 3 rejects earlier binaries before microscopy or
constrained calls. The fit-major constraint layout from ABI 2 is unchanged.

`tests/model_derivatives.cpp` independently checks the existing kernel derivatives
with finite differences. Compile it with `cl /EHsc /I. tests/model_derivatives.cpp`
in an initialized MSVC environment. CUDA tests are marked `gpu`; forwarding mocks
alone are not numerical validation. Linux execution is pending. The Linux script
requires an approved supported host compiler/toolkit and performs no installation.

Set `GPUFIT_MAX_CHUNK_SIZE` to a positive integer before starting Python to reduce
the native chunk capacity. It never increases the device-memory-derived limit.
The multi-chunk tests verify that each arbitrary-coordinate model reads the
correct full-batch coordinates and size table, including the final partial chunk.

An independently locked local integration project can use this same toolchain:
`powershell -File scripts/build_windows.ps1 -Action sync -Architectures 86 -ProjectPath ../localize-psf/build/local-gpu`.
The script applies native build settings specifically to `pygpufit` in that case.
Configure the documented self-hosted CI runner labels before requesting native
jobs. Local checks have run; hosted/self-hosted workflow execution is pending.

Packaging follows the [scikit-build-core configuration reference](https://scikit-build-core.readthedocs.io/en/stable/configuration/index.html),
including CMake configuration during sdist creation and a platform wheel for a
Python-independent ctypes ABI.
