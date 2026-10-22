"""Human-readable report (Markdown) and machine-readable summary (JSON) for AI agents."""
from __future__ import annotations

import json

import numpy as np


def _f(v, n=3):
    return "-" if v is None or (isinstance(v, float) and not np.isfinite(v)) else f"{v:.{n}f}"


def write_report(path, results, refs, cfg):
    L = ["# XRD analysis report", ""]
    L.append(f"Samples: {len(results)} | References: {', '.join(r.name for r in refs) or 'none'} | "
             f"Plot range: {cfg['x_range'][0]}-{cfg['x_range'][1]} deg 2theta")
    if cfg["instrument_fwhm"] == 0:
        L.append("")
        L.append("> Crystallite sizes are NOT corrected for instrument broadening (instrument_fwhm = 0). "
                 "Treat them as lower bounds and compare samples only relative to each other.")
    for res in results:
        L += ["", f"## {res['name']}", "",
              f"- File: `{res['file']}` ({res['format']})",
              f"- Range {res['range'][0]:.2f}-{res['range'][1]:.2f} deg, step {res['step']:.4f} deg, "
              f"{res['n_points']} points, wavelength {res['wavelength']:.5f} A",
              f"- Peaks found: {len(res['peaks'])}; signal-to-noise of strongest peak ~{res['snr']:.0f}"]
        if res["checks"]:
            L += ["", "**Checks**", ""] + [f"- {c}" for c in res["checks"]]
        for m in res["matches"]:
            found = sum(l["found"] for l in m["lines"])
            L += ["", f"**Match with {m['reference']}**: {found}/{len(m['lines'])} strong lines found, "
                      f"weighted score {m['score']:.2f}, median shift {_f(m['median_shift'])} deg", "",
                  "| hkl | ref 2theta | ref I | found | obs 2theta | shift |", "|---|---|---|---|---|---|"]
            for l in m["lines"]:
                L.append(f"| {l['hkl']} | {l['ref_2theta']:.3f} | {l['ref_I']:.0f} | {'yes' if l['found'] else 'no'} | "
                         f"{_f(l['obs_2theta'])} | {_f(l['shift'])} |")
        sm = res.get("search")
        if sm:
            from .searchmatch import describe
            L += ["", "**Phase search (search-match)**", ""]
            for msg in res.get("search_log", [])[:6]:
                L.append(f"- {msg}")
            if sm["accepted"]:
                L += ["", "Proposed phase assemblage:", ""]
                for a in sm["accepted"]:
                    miss = (f"; strong lines missing at {', '.join(f'{t:.1f}' for t in a['strong_missing'])}"
                            if a["strong_missing"] else "")
                    lat = ""
                    if a.get("strain_axes"):
                        lat = ", lattice vs reference " + "/".join(f"{v * 100:+.1f}%" for v in a["strain_axes"]) + " (by axis)"
                    elif abs(a.get("strain", 0)) >= 0.002:
                        lat = f", lattice vs reference {a['strain'] * 100:+.1f}%"
                    L.append(f"- **{describe(a)}** - FOM {a['fom']:.2f} (lines found {a['n_hit']}/{a['n_lines']}, "
                             f"completeness {a['completeness']:.2f}, explains {a['coverage'] * 100:.0f}% of peak "
                             f"intensity, zero shift {a['shift']:+.2f} deg{lat}){miss}")
            else:
                L.append("- no candidate reached the acceptance threshold (FOM 0.45)")
            if sm["ranking"]:
                L += ["", "| rank | candidate | FOM | completeness | coverage | shift |", "|---|---|---|---|---|---|"]
                for k, r in enumerate(sm["ranking"][:10], 1):
                    L.append(f"| {k} | {describe(r)} | {r['fom']:.2f} | {r['completeness']:.2f} | "
                             f"{r['coverage']:.2f} | {r['shift']:+.2f} |")
            if sm["unexplained"]:
                L.append("")
                L.append("Peaks not explained by the proposed phases (2theta): "
                         + ", ".join(f"{t:.2f}" for t in sm["unexplained"][:20]))
            L.append("")
            L.append("_Search-match ranks candidates by peak positions; confirm with chemistry, other techniques "
                     "or refinement before claiming a phase._")
        for rn, tc in res["texture"].items():
            L += ["", f"**Texture coefficients vs {rn}** (TC > 1 = preferred orientation)", "",
                  "| hkl | 2theta | TC |", "|---|---|---|"]
            L += [f"| {t['hkl']} | {t['obs_2theta']:.2f} | {t['TC']:.2f} |" for t in tc]
        if res["ratios"]:
            L += ["", "**Intensity ratios**", ""]
            L += [f"- {r['label']}: {_f(r['value'])} ({r['note']})" for r in res["ratios"]]
        L += ["", "**Peak list** (strongest 15)", "",
              "| 2theta | d (A) | rel I | FWHM | size (nm) | assignment |", "|---|---|---|---|---|---|"]
        for p in sorted(res["peaks"], key=lambda q: -q["height"])[:15]:
            asg = f"{p.get('phase', '')} ({p['hkl']})" if p.get("hkl") else ""
            L.append(f"| {p['center']:.3f} | {p['d']:.4f} | {p.get('rel_intensity', 0):.0f} | {p['fwhm']:.3f} | "
                     f"{_f(p.get('size_nm'), 1)} | {asg} |")
        if refs and res["unassigned"]:
            L += ["", "Unassigned peaks (2theta): " + ", ".join(f"{t:.2f}" for t in res["unassigned"][:20])]
    L += ["", "---", "Notes: intensities are from background-subtracted data (SNIP). Phase scores only check "
          "the references you supplied; a peak list match is qualitative, not a Rietveld refinement."]
    with open(path, "w", encoding="utf-8") as f:
        f.write("\n".join(L) + "\n")


def _clean(o):
    from .reference import Reference
    if isinstance(o, Reference):
        return {"name": o.name, **{k: v for k, v in o.meta.items()}}
    if isinstance(o, dict):
        return {k: _clean(v) for k, v in o.items() if k not in ("x", "y", "bg", "net")}
    if isinstance(o, (list, tuple)):
        return [_clean(v) for v in o]
    if isinstance(o, (np.floating, float)):
        return None if not np.isfinite(o) else round(float(o), 5)
    if isinstance(o, np.integer):
        return int(o)
    if isinstance(o, np.bool_):
        return bool(o)
    return o


def write_summary_json(path, results, refs, cfg):
    data = {"references": [{"name": r.name, "source": r.source, **r.meta} for r in refs],
            "x_range": cfg["x_range"], "samples": _clean(results)}
    with open(path, "w", encoding="utf-8") as f:
        json.dump(data, f, ensure_ascii=False, indent=1)
