"""Candidate phases for search-match: the built-in library and COD."""
from __future__ import annotations

import hashlib
import json
import os

import numpy as np

from .reference import Reference, reference_from_cif

LIB_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))),
                       "references", "library")
PATTERN_CACHE = os.path.join(os.path.expanduser("~"), ".xrd-analysis", "patterns")


def _cached_reference(cif, wl, tmin, tmax, name=None, extra_meta=None):
    """reference_from_cif with an on-disk cache (big MOF cells take seconds to compute)."""
    h = hashlib.sha1(open(cif, "rb").read()).hexdigest()[:16]
    key = f"v2_{h}_{wl:.5f}_{tmin:.1f}_{tmax:.1f}.json"
    p = os.path.join(PATTERN_CACHE, key)
    if os.path.exists(p):
        try:
            d = json.load(open(p, encoding="utf-8"))
            ref = Reference(name or d["name"], np.array(d["t"]), np.array(d["i"]), d["hkl"], np.array(d["d"]),
                            source=cif, meta=d["meta"],
                            hkl_int=np.array(d["hkl_int"], int).reshape(-1, 3) if d.get("hkl_int") else None)
            ref.meta.update(extra_meta or {})
            return ref
        except Exception:  # noqa: BLE001
            pass
    ref = reference_from_cif(cif, wl, tmin, tmax, name=name)
    ref.meta.update(extra_meta or {})
    try:
        os.makedirs(PATTERN_CACHE, exist_ok=True)
        json.dump({"name": ref.name, "t": ref.two_theta.tolist(), "i": ref.intensity.tolist(), "hkl": ref.hkl,
                   "d": ref.d.tolist(), "meta": ref.meta,
                   "hkl_int": ref.hkl_int.tolist() if ref.hkl_int is not None else None}, open(p, "w", encoding="utf-8"))
    except OSError:
        pass
    return ref


def library_index():
    p = os.path.join(LIB_DIR, "index.json")
    if not os.path.exists(p):
        return []
    return json.load(open(p, encoding="utf-8"))


def library_candidates(wl, tmin, tmax, groups=None):
    out, log = [], []
    for e in library_index():
        if groups and not set(groups) & set(e.get("groups", [])):
            continue
        cif = os.path.join(LIB_DIR, e["file"])
        try:
            meta = {k: e[k] for k in ("formula", "source", "cod_id", "groups") if k in e}
            out.append(_cached_reference(cif, wl, tmin, tmax, name=e["name"], extra_meta=meta))
        except Exception as ex:  # noqa: BLE001
            log.append(f"library {e['file']}: {ex}")
    log.insert(0, f"built-in library: {len(out)} phases compared")
    return out, log


def cod_candidates(elements, subsets, text, wl, tmin, tmax, max_candidates=80):
    from .cod import deduplicate, fetch_cif, search_cod
    if not elements and not text:
        raise ValueError("COD search needs elements (e.g. --search Zn,O,S,H) or a keyword")
    entries = deduplicate(search_cod(elements=elements, subsets=subsets, text=text or None))
    log = [f"COD: {len(entries)} distinct entries for {'-'.join(elements or [])} {text or ''}".strip()]
    # prefer smaller cells first (faster, and most common inorganic phases), keep the cap
    entries.sort(key=lambda e: (e["volume"] or 1e9))
    if len(entries) > max_candidates:
        log.append(f"COD: only the {max_candidates} entries with the smallest cells are compared "
                   "(raise search.max_candidates or narrow the elements)")
        entries = entries[:max_candidates]
    import sys
    from concurrent.futures import ThreadPoolExecutor
    print(f"COD: {len(entries)} structures to compare (downloaded once, then cached)", file=sys.stderr, flush=True)
    done = [0]

    def get(e):
        try:
            return e, fetch_cif(e["id"]), None
        except Exception as ex:  # noqa: BLE001
            return e, None, ex
        finally:
            done[0] += 1
            if done[0] % 10 == 0 or done[0] == len(entries):
                print(f"  downloaded {done[0]}/{len(entries)}", file=sys.stderr, flush=True)

    with ThreadPoolExecutor(max_workers=6) as pool:  # polite parallelism
        fetched = list(pool.map(get, entries))
    print("  calculating patterns...", file=sys.stderr, flush=True)
    out = []
    for e, p, ex in fetched:
        if ex is not None:
            log.append(f"COD {e['id']}: skipped ({str(ex)[:80]})")
            continue
        try:
            nm = e["name"] if e["name"] and len(e["name"]) <= 30 else e["formula"]
            out.append(_cached_reference(p, wl, tmin, tmax, name=f"{nm}", extra_meta={
                "formula": e["formula"], "cod_id": e["id"], "source": "COD"}))
        except Exception as ex:  # noqa: BLE001
            log.append(f"COD {e['id']}: skipped ({str(ex)[:80]})")
    return out, log
