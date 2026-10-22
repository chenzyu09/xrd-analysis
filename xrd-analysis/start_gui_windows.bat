@echo off
chcp 65001 >nul
cd /d "%~dp0"
set PY=
where py >NUL 2>NUL && set PY=py -3
if not defined PY (
  python --version >NUL 2>NUL && set PY=python
)
if not defined PY (
  echo [!] Python not found. Install Python 3.9+ first:
  echo     winget install -e --id Python.Python.3.12
  echo     or https://www.python.org  ^(tick "Add python.exe to PATH"^)
  pause
  exit /b 1
)
echo Using: %PY%
%PY% -m pip install -r xrd-analysis\scripts\requirements.txt streamlit
if errorlevel 1 (
  echo [!] Package install failed. If downloads are slow, run once:
  echo     %PY% -m pip config set global.index-url https://pypi.tuna.tsinghua.edu.cn/simple
  pause
  exit /b 1
)
%PY% -m streamlit run xrd-analysis\scripts\app.py
pause
