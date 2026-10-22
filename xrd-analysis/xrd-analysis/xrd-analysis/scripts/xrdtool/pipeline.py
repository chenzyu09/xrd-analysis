"""End-to-end analysis: read -> process -> peaks -> reference matching -> figures/tables/report."""
from __future__ import annotations

import copy
import json
import os
import re

import numpy as np

from . import processing as P
from .plotting import PALETTE, REF_COLORS, label_y, peaks_plot, stacked_plot
from .readers import read_pattern
from .reference import load_reference, match_reference, texture_coefficients

DEFAULTS = {
    "_help": "Edit and re-run with --config this_file.json. Lists follow the order of input files.",
    "samples": [],              # filled on first run: [{"file", "label", "color", "show"}]
    "references": [],           # [{"file": "x.cif" or stick list, "label": "...", "show": true}]
    "x_range": None,            # [min, max] in 2theta; null = common range of the data
    "background": "snip",       # "snip" | "none"
    "background_width": 1.5,    # deg; ~ width of the broadest real peak (raise for broad humps)
    "subtract_background_in_plot": False,
    "strip_ka2": False,
    "smooth": 0.0,              # deg window for Savitzky-Golay in the plot (0 = off)
    "normalize": "max",         # "max" | "none"
    "intensity_scale": "linear",  # "linear" | "sqrt" | "log" (sqrt/log reveal weak peaks next to huge ones)
    "offset": 1.1,              # vertical spacing between stacked curves (in normalised units)
    "stack_order": "bottom_up", # first file at the bottom; "top_down" puts it on top
    "reference_style": "sticks",  # "sticks" (panel under the data) | "profile" (calc. pattern) | "none"
    "reference_fwhm": 0.1,      # deg, for "profile"
    "reference_scale": 1.0,
    "label_hkl": True,
    "hkl_min_intensity": 8,     # only label sticks above this %
    "guide_lines": False,       # dashed vertical lines at the first reference's strong peaks
    "guide_min_intensity": 20,
    "label_x": 0.98,            # label position as fraction of x range
    "fig_size": [6.0, 6.5],
    "stack_panel_ratio": 4,
    "font_size": 14,
    "frame_width": 1.5,
    "line_width": 1.2,
    "formats": ["png", "pdf"],
    "title": "",
    "peak_min_snr": 5.0,
    "peak_min_rel": 0.02,       # ignore peaks below this fraction of the strongest
    "match_tolerance": 0.2,     # deg 2theta
    "scherrer_K": 0.9,
    "instrument_fwhm": 0.0,     # deg; measure with a Si/LaB6 standard, 0 = not corrected
    "ratios": [],               # e.g. ["003/104"] (hkl of first reference) or ["18.7/44.5"] (2theta)
    "texture": "auto",          # "auto" (if a reference is given) | "off"
    "search": {                 # automatic phase identification (search-match)
        "source": "auto",       # "auto" (built-in library if no reference given) | "library" | "cod" | "none"
        "elements": [],         # for "cod": elements present, e.g. ["Zn", "O", "S", "H"]
        "subsets": True,        # also phases made of a subset of those elements (e.g. ZnO for Zn-O-S-H)
        "text": "",             # for "cod": extra keyword (name, mineral)
        "cod_ids": [],          # COD entries to use directly as references, e.g. ["9008522"]
        "library_groups": [],   # restrict the built-in library, e.g. ["battery", "mof"]; [] = all
        "add_to_plot": None,    # null = add accepted phases only if no reference was given
        "max_candidates": 80,
        "tolerance": 0.15,
        "max_shift": 0.3,
    },
    "origin_font_size": 36,
    "origin_frame_pt": 3,
    "origin_line_pt": 2,
}

# Common substrate / current-collector lines (Cu K-alpha, strongest first)
SUBSTRATES = {
    "Cu foil": [43.30, 50.43, 74.13],
    "Al foil": [38.47, 44.74, 65.13, 78.23],
    "Ni (foam/foil)": [44.51, 51.85, 76.37],
    "Ti foil": [40.17, 35.09, 38.42, 53.00],
    "Zn foil": [43.23, 36.30, 38.99, 54.34],
    "graphite / carbon paper/cloth": [26.5, 54.6],
    "stainless steel (austenitic)": [43.6, 50.8, 74.7],
    "Si (wafer/standard)": [28.44, 47.30, 56.12],
    "Kapton/Al holder": [38.47, 44.74],
}


def merge_config(user):
    cfg = copy.deepcopy(DEFAULTS)
    for k, v in (user or {}).items():
        cfg[k] = v
    return cfg


def _hkl_key(h):
    return re.sub(r"[()\s,]", "", str(h)).split("/")[0]


def _ratio(spec, peaks, ref, tol):
    a, b = spec.split("/")
    vals = []
    for part in (a, b):
        part = part.strip()
        t = None
        try:
            t = float(part)
            if t < 2.0 or "." not in part:  # looks like an hkl such as 003
                t = None
        except ValueError:
            pass
        if t is None and ref is not None:
            for tt, h in zip(ref.two_theta, ref.hkl):
                if _hkl_key(h) == part or part in str(h).split("/"):
                    t = tt
                    break
        if t is None:
            return {"label": spec, "value": None, "note": f"'{part}' not found in reference"}
        cand = [p for p in peaks if abs(p["center"] - t) <= max(tol, 0.3)]
        if not cand:
            return {"label": spec, "value": None, "note": f"no peak near {t:.2f} deg"}
        p = max(cand, key=lambda q: q["height"])
        vals.append((p["area"] if p.get("fitted") and np.isfinite(p["area"]) else p["height"], p["center"]))
    note = f"I({a})@{vals[0][1]:.2f} / I({b})@{vals[1][1]:.2f}; areas if fitted else heights"
    return {"label": spec, "value": vals[0][0] / vals[1][0] if vals[1][0] else None, "note": note}


def analyse_one(pat, refs, cfg):
    x, y = pat.two_theta, pat.intensity.astype(float)
    step = pat.step
    checks = []
    if cfg["strip_ka2"]:
        y = P.strip_ka2(x, y, wl1=pat.wavelength)
    bg = P.snip_background(y, cfg["background_width"], step) if cfg["background"] == "snip" else np.zeros_like(y)
    net = np.clip(y - bg, 0, None)
    is_cu = abs(pat.wavelength - 1.5406) < 0.002
    peaks, noise, n_ka2 = P.find_xrd_peaks(x, net, cfg["peak_min_snr"], cfg["peak_min_rel"],
                                           ka2=is_cu and not cfg["strip_ka2"])
    for p in peaks:
        p["d"] = float(P.d_spacing(p["center"], pat.wavelength))
        p["size_nm"] = float(P.scherrer(p["fwhm"], p["center"], pat.wavelength, cfg["scherrer_K"],
                                        cfg["instrument_fwhm"])) if p.get("fitted") else np.nan
    # ---- data quality checks
    top = net.max() if len(net) else 0
    snr = top / noise if noise > 0 else np.inf
    if snr < 20:
        checks.append(f"weak signal: strongest peak is only ~{snr:.0f}x the noise; consider longer counting time")
    tot = np.trapezoid(y) if hasattr(np, "trapezoid") else np.trapz(y)
    bga = np.trapezoid(bg) if hasattr(np, "trapezoid") else np.trapz(bg)
    bg_frac = bga / tot if tot > 0 else 0
    if bg_frac > 0.85 and cfg["background"] == "snip":
        checks.append(f"background is {100 * bg_frac:.0f}% of total intensity: low crystallinity, amorphous "
                      "content, or fluorescence (Fe/Co/Mn-rich samples under a Cu tube)")
    hump = _amorphous_hump(x, bg)
    if hump:
        checks.append(f"broad hump centred near {hump:.1f} deg in the background "
                      "(amorphous phase, glass holder or polymer binder/separator)")
    fitted = [p for p in peaks if p.get("fitted")]
    if fitted:
        narrow = min(p["fwhm"] for p in fitted)
        if narrow < 3 * step:
            checks.append(f"step size {step:.3f} deg is coarse for the narrowest peak (FWHM {narrow:.3f} deg); "
                          "FWHM/size values for sharp peaks are unreliable")
    if np.any(y < 0):
        checks.append("negative intensities in file (data already background-subtracted?)")
    if len(x) > 1 and np.any(np.diff(x) <= 0):
        checks.append("2theta axis not strictly increasing (merged scans?)")
    # ---- reference matching
    matches, texture, assigned = [], {}, {}
    for r in refs:
        m = match_reference(peaks, r, cfg["match_tolerance"], x[0], x[-1])
        matches.append(m)
        for ln in m["lines"]:
            if ln["found"]:
                k = ln["peak_index"]
                if k not in assigned or ln["ref_I"] > assigned[k][2]:
                    assigned[k] = (ln["hkl"], r.name, ln["ref_I"])
        if np.isfinite(m["median_shift"]) and abs(m["median_shift"]) > 0.08 and m["score"] > 0.5:
            sign = "higher" if m["median_shift"] > 0 else "lower"
            why = (m.get("shift_model") or {}).get("note") or (
                "lattice change (doping / solid solution / lithiation) or sample-height error")
            checks.append(f"peaks sit {abs(m['median_shift']):.2f} deg {sign} than {r.name}: {why}")
        if cfg["texture"] != "off":
            tc = texture_coefficients(peaks, r, cfg["match_tolerance"], tmin=x[0], tmax=x[-1])
            if tc:
                texture[r.name] = tc
    for k, (h, rn, _) in assigned.items():
        peaks[k]["hkl"], peaks[k]["phase"] = h, rn
    # unexplained peaks vs common substrates
    unexplained = [p for i, p in enumerate(peaks) if i not in assigned]
    subs = []
    ref_lines = np.concatenate([r.two_theta[r.intensity >= 20] for r in refs]) if refs else np.array([])
    for sname, lines in SUBSTRATES.items():
        if len(ref_lines) and np.min(np.abs(ref_lines - lines[0])) < 0.3:
            continue  # the substrate is (or coincides with) one of the reference phases
        hit = [t for t in lines[:3] if any(abs(p["center"] - t) < 0.25 for p in (unexplained if refs else peaks))]
        if len(hit) >= 2 or (len(hit) == 1 and hit[0] == lines[0] and len(lines) <= 2):
            subs.append((sname, hit))
    if subs and refs:
        checks.append("unassigned peaks match " + "; ".join(f"{s} ({', '.join(f'{t:.1f}' for t in h)})"
                                                         for s, h in subs) + " - substrate/current collector/holder?")
    ratios = [_ratio(s, peaks, refs[0] if refs else None, cfg["match_tolerance"]) for s in cfg["ratios"]]
    return {"name": pat.name, "file": pat.source, "format": pat.fmt, "wavelength": pat.wavelength,
            "range": [float(x[0]), float(x[-1])], "step": step, "n_points": len(x),
            "x": x, "y": y, "bg": bg, "net": net, "peaks": peaks, "noise": noise, "snr": snr,
            "ka2_doublets_merged": n_ka2, "background_fraction": bg_frac,
            "checks": checks, "matches": matches, "texture": texture, "ratios": ratios,
            "unassigned": [p["center"] for p in unexplained], "substrate_hits": subs, "meta": pat.meta}


def _amorphous_hump(x, bg):
    if len(x) < 50:
        return None
    k = max(5, len(x) // 50)
    b = np.convolve(bg, np.ones(k) / k, mode="same")
    i = int(np.argmax(b[k:-k])) + k
    edge = np.median(np.r_[b[k:2 * k], b[-2 * k:-k]])
    if b[i] > 1.6 * edge and 0.1 * len(x) < i < 0.9 * len(x):
        return float(x[i])
    return None


def run(files, cfg, out_dir, names=None):
    os.makedirs(out_dir, exist_ok=True)
    pats = [read_pattern(f) for f in files]
    seen = {}
    for p in pats:  # make sample names unique (same file name in different folders)
        n = seen.get(p.name, 0)
        seen[p.name] = n + 1
        if n:
            p.name = f"{p.name}_{n + 1}"
    if not cfg["samples"]:
        cfg["samples"] = [{"file": os.path.abspath(f), "label": (names[i] if names and i < len(names) else p.name),
                           "color": PALETTE[i % len(PALETTE)], "show": True}
                          for i, (f, p) in enumerate(zip(files, pats))]
    wl = pats[0].wavelength if pats else 1.540593
    lo = min(p.two_theta[0] for p in pats)
    hi = max(p.two_theta[-1] for p in pats)
    if not cfg["x_range"]:
        cfg["x_range"] = [float(np.floor(max(p.two_theta[0] for p in pats))),
                          float(np.ceil(min(p.two_theta[-1] for p in pats)))]
    refs = []
    tmin_ref, tmax_ref = max(1.0, lo - 1), min(150, hi + 1)
    for r in cfg["references"]:
        ref = load_reference(r["file"], wl, tmin_ref, tmax_ref, name=r.get("label"))
        ref.show = r.get("show", True)
        refs.append(ref)
    sc = dict(DEFAULTS["search"], **(cfg.get("search") or {}))
    cfg["search"] = sc
    search_log = []
    if sc["cod_ids"]:
        from .cod import fetch_cif
        for cid in sc["cod_ids"]:
            try:
                ref = load_reference(fetch_cif(cid), wl, tmin_ref, tmax_ref)
                ref.show = True
                refs.append(ref)
            except Exception as e:  # noqa: BLE001
                search_log.append(f"COD {cid}: {e}")
    print(f"Reading and analysing {len(pats)} pattern(s)...", flush=True)
    results = [analyse_one(p, refs, cfg) for p in pats]
    source = sc["source"]
    if source == "auto":
        source = "library" if not refs else "none"
    if source in ("library", "cod"):
        from .candidates import cod_candidates, library_candidates
        from .searchmatch import search_match
        try:
            if source == "library":
                cands, log = library_candidates(wl, tmin_ref, tmax_ref, sc["library_groups"])
            else:
                cands, log = cod_candidates(sc["elements"], sc["subsets"], sc["text"], wl, tmin_ref, tmax_ref,
                                            sc["max_candidates"])
            search_log += log
        except Exception as e:  # noqa: BLE001
            cands = []
            search_log.append(f"search failed: {e}")
            if source == "cod":
                cands, log = library_candidates(wl, tmin_ref, tmax_ref)
                search_log += ["falling back to the built-in library"] + log
        print(f"Phase search: comparing {len(cands)} candidate phases...", flush=True)
        accepted_all = []
        for res in results:
            sm = search_match(res["peaks"], cands, res["range"][0], res["range"][1], sc["tolerance"],
                              sc["max_shift"])
            res["search"] = sm
            for a in sm["accepted"]:
                if all(a["ref"] is not b for b in accepted_all):
                    accepted_all.append(a["ref"])
        add = sc["add_to_plot"]
        if add is None:
            add = not refs
        if add and accepted_all:
            known = {(r.meta.get("formula"), r.meta.get("spacegroup")) for r in refs}
            for r in accepted_all:
                if (r.meta.get("formula"), r.meta.get("spacegroup")) not in known:
                    r.show = True
                    refs.append(r)
                    # remember it, so editing config.json and re-running keeps the same references
                    cfg["references"].append({"file": r.source, "label": r.name, "show": True,
                                              "found_by": "search-match"})
            results_new = [analyse_one(p, refs, cfg) for p in pats]
            for old, new in zip(results, results_new):
                new["search"] = old["search"]
            results = results_new
    for res in results:
        res["search_log"] = search_log
    for msg in search_log:
        if "failed" in msg or "falling back" in msg:
            print("note:", msg)
    # ---- curves for the stacked plot
    curves = []
    for res, s in zip(results, cfg["samples"]):
        if not s.get("show", True):
            continue
        y = res["net"] if cfg["subtract_background_in_plot"] else res["y"]
        y = P.smooth(y, cfg["smooth"], res["step"])
        if cfg["intensity_scale"] == "sqrt":
            y = np.sqrt(np.clip(y, 0, None))
        elif cfg["intensity_scale"] == "log":
            y = np.log10(np.clip(y, 1, None))
        xmin, xmax = cfg["x_range"]
        m = (res["x"] >= xmin) & (res["x"] <= xmax)
        if cfg["normalize"] == "max" and m.any():
            seg = y[m]
            y = (y - seg.min()) / ((seg.max() - seg.min()) or 1)
        curves.append({"x": res["x"], "y": y, "label": s["label"], "color": s["color"]})
    shown_refs = [r for r in refs if getattr(r, "show", True)]
    print("Writing figures, Excel, Origin script and report...", flush=True)
    paths = stacked_plot(curves, shown_refs, cfg, os.path.join(out_dir, "xrd_stacked"))
    # offsets used in the figure (mirrors stacked_plot) -> for Excel / Origin
    off = 0.0
    order = list(range(len(curves)))
    if cfg["stack_order"] == "top_down":
        order = order[::-1]
    xmin, xmax = cfg["x_range"]
    for i in order:
        c = curves[i]
        m = (c["x"] >= xmin) & (c["x"] <= xmax)
        c["y_plot"] = c["y"] + off
        c["label_y"] = off + label_y(c["x"][m], c["y"][m], c["label"], cfg)
        off += cfg["offset"] * max(1.0, float(np.nanmax(c["y"][m])) if m.any() else 1.0)
    diag = os.path.join(out_dir, "diagnostics")
    os.makedirs(diag, exist_ok=True)
    for res in results:
        best = None
        if refs and res["matches"]:
            best = refs[int(np.argmax([m["score"] for m in res["matches"]]))]
        res["diagnostic_png"] = peaks_plot(res["x"], res["y"], res["bg"], res["peaks"], best, res["name"],
                                           os.path.join(diag, _safe(res["name"])), cfg)
    from .export import write_excel, write_origin_script
    from .report import write_report, write_summary_json
    write_excel(os.path.join(out_dir, "xrd_data.xlsx"), curves, refs, results)
    write_origin_script(os.path.join(out_dir, "origin_plot.py"), curves, shown_refs, cfg, off + 0.1, REF_COLORS)
    write_report(os.path.join(out_dir, "report.md"), results, refs, cfg)
    write_summary_json(os.path.join(out_dir, "summary.json"), results, refs, cfg)
    cfg_path = os.path.join(out_dir, "config.json")
    with open(cfg_path, "w", encoding="utf-8") as f:
        json.dump(cfg, f, ensure_ascii=False, indent=2)
    return {"figures": paths, "out_dir": out_dir, "config": cfg_path, "results": results}


def _safe(s):
    return re.sub(r'[\\/:*?"<>|\s]+', "_", s)
