@echo off
setlocal
cd /d "%~dp0"

echo === Asteroid Watch ===

if not exist "app\public\data\asteroids.json" (
  echo [data] asteroids.json missing - regenerating with Python pipeline...
  python scripts\fetch_sbdb.py || goto :fail
  python scripts\build_elements.py || goto :fail
  python scripts\trajectory.py || goto :fail
)

if not exist "app\public\data\ephemeris\planets.bin" (
  echo [data] ephemeris missing - fetching from JPL Horizons ^(slow^)...
  python scripts\fetch_ephemeris.py || goto :fail
)

cd app
if not exist node_modules (
  echo [app] installing dependencies...
  call npm install || goto :fail
)

echo [app] starting dev server on http://localhost:5174 ...
start "Asteroid Watch server" cmd /k "npm run dev -- --port 5174 --strictPort"

echo Waiting for server...
timeout /t 6 /nobreak >nul

start chrome "http://localhost:5174" 2>nul || start "" "http://localhost:5174"
echo Opened in browser. Close the "Asteroid Watch server" window to stop.
goto :eof

:fail
echo.
echo Something failed - see the messages above.
pause
exit /b 1
