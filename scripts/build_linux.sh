#!/usr/bin/env bash
set -euo pipefail
cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.."
command -v uv >/dev/null || { echo 'Install uv first.' >&2; exit 1; }
command -v nvcc >/dev/null || { echo 'Install a supported CUDA toolkit/host compiler first.' >&2; exit 1; }
uv lock --check
# Read the same policy used by the Windows script and standalone CMake preset.
mapfile -t settings < <(uv run --locked --no-sync python -c '
import json
p=json.load(open("CMakePresets.json"))["configurePresets"][0]
for k,v in p["cacheVariables"].items(): print("--config-setting=cmake.define."+k+"="+v)
')
export CMAKE_GENERATOR=Ninja
mkdir -p build
uv build --wheel --sdist --config-setting=build-dir=build/microscopy "${settings[@]}" 2>&1 | tee build/build_linux.log
