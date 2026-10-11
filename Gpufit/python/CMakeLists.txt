# Explicit developer convenience target. Never entered by the PEP 517 build.
add_custom_target(PYTHON_WHEEL
  COMMAND uv build --wheel "${PROJECT_SOURCE_DIR}"
  COMMENT "Build the platform wheel through uv and scikit-build-core")
