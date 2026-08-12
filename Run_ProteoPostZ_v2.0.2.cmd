@echo off
if exist "%~dp0ProteoPostZ_v2.0.2.exe" (
  start "" "%~dp0ProteoPostZ_v2.0.2.exe"
  exit /b
)
echo Cannot find ProteoPostZ_v2.0.2.exe in this folder.
pause
