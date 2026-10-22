"""Crystallography Open Database (COD, www.crystallography.net) access.

COD entries are in the public domain, so CIFs can be downloaded, cached and redistributed.
Only the Python standard library is used.

  search_cod(elements=["Zn", "O"])             entries made of exactly these elements
  search_cod(elements=[...], subsets=True)     ... or of any subset of them
  search_cod(formula="O Zn") / text="MIL-101"  Hill formula / free text (names, minerals, papers)
  fetch_cif(9008522)                           CIF text, cached in ~/.xrd-analysis/cod
"""
from __future__ import annotations

import itertools
import json
import sys
import os
import re
import time
import urllib.error
import urllib.parse
import urllib.request

BASES = ["https://www.crystallography.net/cod"]
CACHE = os.environ.get("XRD_COD_CACHE", os.path.join(os.path.expanduser("~"), ".xrd-analysis", "cod"))
UA = {"User-Agent": "xrd-analysis (https://github.com/chenzyu09/xrd-analysis)"}


class CODUnavailable(RuntimeError):
    pass


def _log(msg):
    print(msg, file=sys.stderr, flush=True)


def _get(path, timeout=30, retries=1):
    last = None
    for base in BASES:
        url = f"{base}/{path}"
        for k in range(retries + 1):
            try:
                req = urllib.request.Request(url, headers=UA)
                with urllib.request.urlopen(req, timeout=timeout) as r:
                    data = r.read()
                txt = data.decode("utf-8", errors="replace")
                if "temporarily unavailable" in txt[:500]:
                    raise CODUnavailable("COD reports: service temporarily unavailable")
                return txt
            except (urllib.error.URLError, TimeoutError, CODUnavailable, ConnectionError) as e:
                last = e
                _log(f"  COD: no answer ({str(e)[:60]}), retrying..." if k < retries else f"  COD: no answer ({str(e)[:60]})")
                time.sleep(1.5 * (k + 1))
    raise CODUnavailable(f"cannot reach COD ({last}). Check the internet connection, or download the CIF "
                         "in a browser from www.crystallography.net and pass it with --ref.")


def _hill(formula):
    """'- Zn O -' -> 'ZnO' style compact formula for display."""
    return " ".join((formula or "").strip("- ").split())


def _query(params):
    params = dict(params, format="json")
    txt = _get("result?" + urllib.parse.urlencode(params))
    try:
        rows = json.loads(txt)
    except json.JSONDecodeError as e:
        raise CODUnavailable(f"unexpected COD answer: {txt[:200]}") from e
    out = []
    for r in rows if isinstance(rows, list) else []:
        def num(k):
            try:
                return float(r.get(k))
            except (TypeError, ValueError):
                return None
        out.append({
            "id": str(r.get("file")),
            "formula": _hill(r.get("formula")),
            "name": r.get("mineral") or r.get("commonname") or r.get("chemname") or "",
            "sg": (r.get("sg") or "").strip(),
            "cell": [num(k) for k in ("a", "b", "c", "alpha", "beta", "gamma")],
            "volume": num("vol"),
            "year": r.get("year"),
            "title": (r.get("title") or "")[:160],
            "journal": r.get("journal"),
        })
    return out


def search_cod(elements=None, subsets=False, formula=None, text=None, exclude=None, max_results=500):
    """Search COD. elements: exact element set (or subsets of it with subsets=True)."""
    results = {}
    queries = []
    if elements:
        els = [e.strip().capitalize() for e in elements if e.strip()]
        sets = [els]
        if subsets:
            sets = [list(c) for n in range(1, len(els) + 1) for c in itertools.combinations(els, n)]
        for s in sets:
            q = {f"el{i + 1}": e for i, e in enumerate(s[:8])}
            q.update(strictmin=len(s), strictmax=len(s))
            queries.append(q)
    if formula:
        queries.append({"formula": formula})
    if text:
        queries = [dict(q, text=text) for q in queries] if queries else [{"text": text}]
    for n, q in enumerate(queries, 1):
        els = [v for k, v in q.items() if k.startswith("el")]
        _log(f"COD: searching {'-'.join(els) or q.get('formula') or q.get('text')} ({n}/{len(queries)})...")
        for i, e in enumerate(exclude or []):
            q[f"nel{i + 1}"] = e
        for row in _query(q):
            results[row["id"]] = row
        if len(results) >= max_results:
            break
    return list(results.values())[:max_results]


def deduplicate(entries):
    """One entry per (formula, space group, rounded cell) - COD has many redeterminations."""
    seen = {}
    for e in entries:
        cell = tuple(round(v, 1) if v else None for v in e["cell"][:3])
        key = (e["formula"], e["sg"].replace(" ", ""), cell)
        if key not in seen or (e.get("year") or 0) > (seen[key].get("year") or 0):
            seen[key] = e
    return list(seen.values())


def cif_path(cod_id):
    return os.path.join(CACHE, f"{cod_id}.cif")


def fetch_cif(cod_id, cache=True):
    cod_id = str(cod_id).strip()
    if not re.fullmatch(r"\d{7}", cod_id):
        raise ValueError(f"COD ids have 7 digits, got {cod_id!r}")
    p = cif_path(cod_id)
    if cache and os.path.exists(p) and os.path.getsize(p) > 100:
        return p
    txt = _get(f"{cod_id}.cif")
    if "_cell_length_a" not in txt:
        raise CODUnavailable(f"COD {cod_id}: no CIF returned")
    os.makedirs(CACHE, exist_ok=True)
    with open(p, "w", encoding="utf-8") as f:
        f.write(txt)
    return p
