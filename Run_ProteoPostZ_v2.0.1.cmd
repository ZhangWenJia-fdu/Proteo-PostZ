@echo off
if exist "%~dp0ProteoPostZ_v2.0.1.exe" (
  start "" "%~dp0ProteoPostZ_v2.0.1.exe"
  exit /b
)
echo Cannot find ProteoPostZ launcher exe in this folder.
pause
