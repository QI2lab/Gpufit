#!/bin/sh
# Compatibility entry point for the supported uv/CMake build.
exec bash "$(dirname "$0")/scripts/build_linux.sh" "$@"
