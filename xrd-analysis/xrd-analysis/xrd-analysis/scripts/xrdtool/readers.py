"""Readers for common powder XRD file formats.

Every reader returns a Pattern: 2theta (deg), intensity, and metadata.
Supported:
  * Rigaku SmartLab / RAS text export (.ras, .txt with *HEADER lines)
  * Rigaku SmartLab Studio .rasx (zip archive)
  * Bruker RAW v1 ("RAW ") and v3 ("RAW1.01")
  * Generic two-column text (.xy, .xye, .txt, .dat, .csv)
  * Excel (.xlsx) with two numeric columns
"""
from __future__ import annotations

import io
import os
import re
import struct
import zipfile
from dataclasses import dataclass, field

import numpy as np

CU_KA1 = 1.540593


@dataclass
class Pattern:
    name: str
    two_theta: np.ndarray
    intensity: np.ndarray
    wavelength: float = CU_KA1
    source: str = ""
    fmt: str = ""
    meta: dict = field(default_factory=dict)

    @property
    def step(self) -> float:
        return float(np.median(np.diff(self.two_theta))) if len(self.two_theta) > 1 else 0.0


def _decode(raw: bytes) -> str:
    for enc in ("utf-8-sig", "gb18030", "latin-1"):
        try:
            return raw.decode(enc)
        except UnicodeDecodeError:
            continue
    return raw.decode("latin-1", errors="ignore")


_NUM = re.compile(r"^[\s,;]*[-+]?(\d+\.?\d*|\.\d+)([eE][-+]?\d+)?([\s,;\t]+[-+]?(\d+\.?\d*|\.\d+)([eE][-+]?\d+)?)+[\s,;]*$")


def _numeric_rows(lines):
    rows = []
    for ln in lines:
        if _NUM.match(ln):
            parts = re.split(r"[\s,;]+", ln.strip().strip(",;"))
            try:
                rows.append([float(p) for p in parts])
            except ValueError:
                continue
    return rows


def _two_columns(rows):
    if not rows:
        raise ValueError("no numeric data found")
    n = min(len(r) for r in rows)
    arr = np.array([r[:n] for r in rows], float)
    # pick the first column that looks like a monotonic 2theta axis
    for c in range(arr.shape[1] - 1):
        x = arr[:, c]
        if np.all(np.diff(x) > 0) and 0 <= x.min() and x.max() <= 180:
            return x, arr[:, c + 1]
    return arr[:, 0], arr[:, 1]


# ---------------------------------------------------------------- Rigaku text
def read_rigaku_text(text: str, name: str):
    meta = {}
    data_lines = []
    in_block = None
    for ln in text.splitlines():
        s = ln.strip()
        if s.startswith("*"):
            m = re.match(r'\*(\S+)\s+"?(.*?)"?\s*$', s)
            if m:
                meta[m.group(1)] = m.group(2)
            if s.startswith("*RAS_INT_START"):
                in_block = True
            elif s.startswith("*RAS_INT_END"):
                in_block = False
            continue
        if s.startswith("#"):
            continue
        if in_block is False:
            continue
        data_lines.append(s)
    rows = _numeric_rows(data_lines)
    arr = np.array([r[:2] for r in rows], float)
    wl = float(meta.get("HW_XG_WAVE_LENGTH_ALPHA1", CU_KA1) or CU_KA1)
    sample = meta.get("FILE_SAMPLE", "")
    return Pattern(name, arr[:, 0], arr[:, 1], wl, fmt="Rigaku RAS/TXT",
                   meta={"target": meta.get("HW_XG_TARGET_NAME", ""),
                         "sample": sample,
                         "scan_start": meta.get("MEAS_SCAN_START"),
                         "scan_stop": meta.get("MEAS_SCAN_STOP"),
                         "scan_step": meta.get("MEAS_SCAN_STEP"),
                         "scan_speed": meta.get("MEAS_SCAN_SPEED"),
                         "time": meta.get("MEAS_SCAN_START_TIME"),
                         "instrument": meta.get("FILE_SYSTEM_NAME", "")})


def read_rasx(path: str, name: str):
    with zipfile.ZipFile(path) as z:
        prof = [n for n in z.namelist() if re.search(r"Profile\d*\.txt$", n)]
        if not prof:
            raise ValueError("rasx: no Profile*.txt inside")
        text = _decode(z.read(sorted(prof)[0]))
        cond = [n for n in z.namelist() if "Conditions" in n and n.endswith(".xml")]
        xml = _decode(z.read(cond[0])) if cond else ""
    rows = _numeric_rows(text.splitlines())
    arr = np.array(rows, float)
    inten = arr[:, 1] * (arr[:, 2] if arr.shape[1] > 2 else 1.0)
    wl = CU_KA1
    m = re.search(r"WavelengthKalpha1>([\d.]+)<", xml)
    if m:
        wl = float(m.group(1))
    tgt = re.search(r"<TargetName>(\w+)<", xml)
    return Pattern(name, arr[:, 0], inten, wl, fmt="Rigaku RASX",
                   meta={"target": tgt.group(1) if tgt else ""})


# ---------------------------------------------------------------- Bruker RAW
def read_bruker_raw(b: bytes, name: str):
    if b[:4] == b"RAW ":
        steps = struct.unpack_from("<I", b, 4)[0]
        step = struct.unpack_from("<f", b, 12)[0]
        start = struct.unpack_from("<f", b, 24)[0]
        wl = struct.unpack_from("<f", b, 72)[0] or CU_KA1
        y = np.frombuffer(b, "<f4", count=steps, offset=156).astype(float)
        x = start + step * np.arange(steps)
        return Pattern(name, x, y, float(wl), fmt="Bruker RAW v1")
    if b[:7] == b"RAW1.01":
        nranges = struct.unpack_from("<I", b, 12)[0]
        off = 712
        xs, ys = [], []
        wl = CU_KA1
        for _ in range(max(nranges, 1)):
            hlen, steps = struct.unpack_from("<II", b, off)
            start = struct.unpack_from("<d", b, off + 16)[0]
            step = struct.unpack_from("<d", b, off + 176)[0]
            w = struct.unpack_from("<d", b, off + 240)[0]
            if 0.5 < w < 3:
                wl = w
            extra = struct.unpack_from("<I", b, off + 256)[0]
            doff = off + hlen + extra
            ys.append(np.frombuffer(b, "<f4", count=steps, offset=doff).astype(float))
            xs.append(start + step * np.arange(steps))
            off = doff + 4 * steps
            if off >= len(b):
                break
        return Pattern(name, np.concatenate(xs), np.concatenate(ys), wl, fmt="Bruker RAW v3")
    raise ValueError("unsupported binary .raw (only Bruker RAW v1/v3 are read; "
                     "for Rigaku .raw please use the .txt/.ras/.rasx export)")


def read_rigaku_raw(b: bytes, name: str):
    """Rigaku SmartLab binary .raw (magic 'FI'); layout found by comparison with .rasx exports."""
    start, stop, step = struct.unpack_from("<fff", b, 2962)
    if not (0 <= start < stop <= 180 and 0 < step < 1):
        raise ValueError("unrecognised Rigaku .raw layout; use the .txt/.ras/.rasx export")
    n = int(round((stop - start) / step)) + 1
    off = len(b) - 4 * n
    wl = struct.unpack_from("<d", b, 1212)[0]
    if not 0.5 < wl < 3:
        wl = CU_KA1
    y = np.frombuffer(b, "<f4", count=n, offset=off).astype(float)
    x = start + step * np.arange(n)
    return Pattern(name, x, y, wl, fmt="Rigaku RAW (binary)")


# ---------------------------------------------------------------- generic
def read_generic_text(text: str, name: str):
    rows = _numeric_rows(text.splitlines())
    x, y = _two_columns(rows)
    return Pattern(name, x, y, fmt="two-column text")


def read_excel(path: str, name: str):
    import openpyxl
    wb = openpyxl.load_workbook(path, data_only=True, read_only=True)
    ws = wb.worksheets[0]
    rows = []
    for r in ws.iter_rows(values_only=True):
        nums = [v for v in r if isinstance(v, (int, float))]
        if len(nums) >= 2:
            rows.append(nums)
    x, y = _two_columns(rows)
    return Pattern(name, x, y, fmt="Excel")


def read_pattern(path: str, name: str | None = None) -> Pattern:
    """Read any supported file; format is detected from content, not only extension."""
    name = name or os.path.splitext(os.path.basename(path))[0]
    ext = os.path.splitext(path)[1].lower()
    if ext in (".xlsx", ".xlsm"):
        p = read_excel(path, name)
    elif zipfile.is_zipfile(path) and ext == ".rasx":
        p = read_rasx(path, name)
    else:
        b = open(path, "rb").read()
        if b[:3] == b"RAW":
            p = read_bruker_raw(b, name)
        elif b[:2] == b"FI" and ext == ".raw":
            p = read_rigaku_raw(b, name)
        else:
            text = _decode(b)
            if re.search(r"^\*(FILE_|MEAS_|HW_|RAS_)", text, re.M):
                p = read_rigaku_text(text, name)
            else:
                p = read_generic_text(text, name)
    order = np.argsort(p.two_theta)
    p.two_theta, p.intensity = p.two_theta[order], p.intensity[order]
    p.source = os.path.abspath(path)
    return p


DATA_EXT = {".txt", ".ras", ".rasx", ".raw", ".xy", ".xye", ".dat", ".csv", ".xlsx"}


def find_patterns(folder: str):
    """Pick one file per measurement, preferring text exports over binary duplicates."""
    pref = [".rasx", ".ras", ".txt", ".xy", ".xye", ".dat", ".csv", ".raw", ".xlsx"]
    groups = {}
    for root, _, files in os.walk(folder):
        for f in files:
            stem, ext = os.path.splitext(f)
            if ext.lower() in DATA_EXT:
                groups.setdefault((root, stem), []).append(os.path.join(root, f))
    picked = []
    for (_, _), paths in sorted(groups.items()):
        paths.sort(key=lambda p: pref.index(os.path.splitext(p)[1].lower()))
        for p in paths:
            try:
                read_pattern(p)
                picked.append(p)
                break
            except Exception:
                continue
    return picked
