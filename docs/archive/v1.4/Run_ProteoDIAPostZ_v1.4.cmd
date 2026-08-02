@echo off
for %%F in ("%~dp0ProteoDIAPostZ*.exe") do start "" "%%~fF" & exit /b
echo Cannot find ProteoDIAPostZ launcher exe in this folder.
pause