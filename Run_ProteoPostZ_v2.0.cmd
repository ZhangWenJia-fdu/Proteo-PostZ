@echo off
for %%F in ("%~dp0ProteoPostZ*.exe") do start "" "%%~fF" & exit /b
echo Cannot find ProteoPostZ launcher exe in this folder.
pause
