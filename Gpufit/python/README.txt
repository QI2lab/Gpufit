QI2lab GPUfit Python binding

Python >=3.12 and float32 NumPy host buffers are required. The authoritative
pyproject.toml and uv.lock live at the repository root. See ../../docs/qi2lab_build.md
for approved prerequisites, uv sync, native rebuilds, wheel installation and tests.

Windows: powershell -File scripts/build_windows.ps1
Linux: sh scripts/build_linux.sh (execution remains pending on Linux hardware)

The QI2lab arbitrary-coordinate models require this fork's binary capability
manifest; a stock GPUfit binary is incompatible. gf.get_capabilities() reports
the loaded binary identity. gf.fit_constrained() exposes native box constraints.

Installation

Currently the wheel file has to be installed locally.

If NumPy is not yet installed, install it using pip from the command line

pip install numpy

Then install pyGpufit from the local folder via:

pip install --no-index --find-links=LocalPathToWheelFile pyGpufit

Examples

See project-root/examples/python folder.

Troubleshooting

A common reason for the error message 'CUDA driver version is insufficient for CUDA runtime version' is an outdated Nvidia graphics driver.
