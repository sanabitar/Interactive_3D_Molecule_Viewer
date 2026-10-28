@echo off
setlocal
cd /d "%~dp0"
set "LAURA_PYTHON=C:\ProgramData\miniconda3\python.exe"
if not exist "%LAURA_PYTHON%" set "LAURA_PYTHON=python"
powershell -NoProfile -Command "try { $r = Invoke-WebRequest -Uri 'http://127.0.0.1:8771/laura-viewer.html' -UseBasicParsing -TimeoutSec 2; if ($r.StatusCode -eq 200) { exit 0 } } catch {}; exit 1" >nul 2>&1
if not errorlevel 1 (
    start "" "http://127.0.0.1:8771/laura-viewer.html"
    exit /b
)
start "" /b powershell -NoProfile -WindowStyle Hidden -Command "Start-Sleep -Seconds 2; Start-Process 'http://127.0.0.1:8771/laura-viewer.html'"
echo LAURA is running at http://127.0.0.1:8771/laura-viewer.html
echo Keep this window open while using LAURA. Close it to stop the server.
echo.
"%LAURA_PYTHON%" -m http.server 8771 --bind 127.0.0.1
if errorlevel 1 pause
