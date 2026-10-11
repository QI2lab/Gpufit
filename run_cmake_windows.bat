@echo off
rem Compatibility entry point for the supported uv/CMake build.
powershell -NoProfile -ExecutionPolicy Bypass -File "%~dp0scripts\build_windows.ps1" %*
exit /b %ERRORLEVEL%
