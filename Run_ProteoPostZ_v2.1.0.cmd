@echo off
if exist "%~dp0ProteoPostZ_v2.1.0.exe" (
  start "" "%~dp0ProteoPostZ_v2.1.0.exe"
  exit /b
)
echo Cannot find ProteoPostZ_v2.1.0.exe in this folder.
pause
