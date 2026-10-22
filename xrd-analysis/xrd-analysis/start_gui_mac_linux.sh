#!/bin/sh
# Starts the browser interface. Uses a private virtual environment (.venv) so it also works where
# the system Python refuses "pip install" (macOS Homebrew, Debian/Ubuntu: externally-managed-environment).
cd "$(dirname "$0")" || exit 1
PY=$(command -v python3 || command -v python)
if [ -z "$PY" ]; then
  echo "[!] Python 3.9+ not found. Install it from https://www.python.org (or: brew install python)"
  exit 1
fi
if [ ! -x .venv/bin/python ]; then
  "$PY" -m venv .venv || { echo "[!] could not create .venv (Debian/Ubuntu: sudo apt install python3-venv)"; exit 1; }
fi
.venv/bin/python -m pip install -q -r xrd-analysis/scripts/requirements.txt streamlit || {
  echo "[!] package install failed. If downloads are slow, try:"
  echo "    .venv/bin/python -m pip config set global.index-url https://pypi.tuna.tsinghua.edu.cn/simple"
  exit 1
}
exec .venv/bin/python -m streamlit run xrd-analysis/scripts/app.py
