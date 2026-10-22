"""Reference patterns: calculated from a CIF, or read from a stick list (CSV / Jade-style text).

A reference is a list of sticks: two_theta, intensity (0-100), hkl label, d.
"""
from __future__ import annotations

import os
import re
from dataclasses import dataclass, field

import numpy as np

from .processing import d_spacing, two_theta_from_d


@dataclass
class Reference:
    name: str
    two_theta: np.ndarray
    intensity: np.ndarray
    hkl: list
    d: np.ndarray
    source: str = ""
    meta: dict = field(default_factory=dict)
    hkl_int: np.ndarray = None  # (n, 3) integer indices of each line (CIF references only)

    def strongest(self, n=10, tmin=None, tmax=None):
        m = np.ones(len(self.two_theta), bool)
        if tmin is not None:
            m &= self.two_theta >= tmin
        if tmax is not None:
            m &= self.two_theta <= tmax
        idx = np.where(m)[0]
        return idx[np.argsort(-self.intensity[idx])][:n]

    def profile(self, x, fwhm=0.1, eta=0.5):
        """Continuous 'calculated pattern' (pseudo-Voigt sticks), max = 1."""
        y = np.zeros_like(x, dtype=float)
        for t, i in zip(self.two_theta, self.intensity):
            if x[0] - 5 * fwhm < t < x[-1] + 5 * fwhm:
                dx2 = 4 * (x - t) ** 2 / fwhm ** 2
                y += i * (eta / (1 + dx2) + (1 - eta) * np.exp(-np.log(2) * dx2))
        return y / y.max() if y.max() > 0 else y


# ---------------------------------------------------------------- from CIF
def _multiplicity(hkl, ops):
    eq = set()
    for op in ops:
        h = tuple(op.apply_to_hkl(list(hkl)))
        eq.add(h)
        eq.add(tuple(-v for v in h))
    return len(eq)


def reference_from_cif(path, wavelength=1.540593, tmin=2.0, tmax=90.0, B=1.0, min_rel=0.1,
                       name=None):
    """Kinematic powder pattern from a CIF (needs the 'gemmi' package).

    I = m * |F|^2 * Lp, with IT92 form factors and an isotropic Debye-Waller factor B (A^2).
    """
    try:
        import gemmi
    except ImportError as e:
        raise ImportError("calculating patterns from CIF needs: pip install gemmi") from e
    # read the text ourselves: gemmi's file opening can fail on Windows paths with non-ASCII characters
    small = gemmi.make_small_structure_from_block(_cif_block(path))
    # occupancies in CIFs are chemical (1.0 on special positions); do not rescale them
    sg = small.find_spacegroup() if hasattr(small, "find_spacegroup") else None
    if sg is None:
        sg = gemmi.find_spacegroup_by_name(small.spacegroup_hm) if small.spacegroup_hm else None
    if sg is None:
        sg = gemmi.SpaceGroup("P 1")
    cell = small.cell
    dmin = wavelength / (2 * np.sin(np.radians(tmax / 2)))
    hkls = gemmi.make_miller_array(cell, sg, dmin, 0, True)
    hkls = np.array(hkls, int).reshape(-1, 3)
    hkls = hkls[np.any(hkls != 0, axis=1)]
    sites = small.get_all_unit_cell_sites()
    frac = np.array([[s.fract.x, s.fract.y, s.fract.z] for s in sites])
    occ = np.array([s.occ for s in sites])
    elems = [s.element for s in sites]
    d = np.array([cell.calculate_d(list(h)) for h in hkls])
    stol2 = 1 / (4 * d ** 2)
    F = np.zeros(len(hkls), complex)
    for el_name in sorted({e.name for e in elems}):
        idx = [i for i, e in enumerate(elems) if e.name == el_name]
        el = gemmi.Element(el_name)
        f = np.array([el.it92.calculate_sf(s) for s in stol2]) if el.it92 else np.full(len(hkls), el.atomic_number, float)
        for chunk in range(0, len(hkls), 4000):
            hh = hkls[chunk:chunk + 4000]
            ph = np.exp(2j * np.pi * hh @ frac[idx].T) @ occ[idx]
            F[chunk:chunk + 4000] += f[chunk:chunk + 4000] * ph
    F *= np.exp(-B * stol2)
    ops = list(sg.operations().sym_ops)
    mult = np.array([_multiplicity(h, ops) for h in hkls])
    tt = two_theta_from_d(d, wavelength)
    th = np.radians(tt / 2)
    lp = (1 + np.cos(2 * th) ** 2) / (np.sin(th) ** 2 * np.cos(th))
    I = mult * np.abs(F) ** 2 * lp
    m = np.isfinite(tt) & (tt >= tmin) & (tt <= tmax)
    tt, I, d, hkls = tt[m], I[m], d[m], hkls[m]
    # merge reflections that coincide (same d)
    order = np.argsort(tt)
    tt, I, d, hkls = tt[order], I[order], d[order], hkls[order]
    out_t, out_i, out_d, out_h, out_hi = [], [], [], [], []
    for t, i, dd, h in zip(tt, I, d, hkls):
        lab = "".join(str(v) if v >= 0 else f"-{-v}" for v in h)
        if out_t and abs(t - out_t[-1]) < 1e-4:
            out_i[-1] += i
            out_h[-1] = out_h[-1] + "/" + lab if out_i[-1] < 2 * i else out_h[-1]
        else:
            out_t.append(t), out_i.append(i), out_d.append(dd), out_h.append(lab), out_hi.append(list(h))
    I = np.array(out_i)
    I = 100 * I / I.max() if len(I) else I
    keep = I >= min_rel
    info = cif_info(path)
    return Reference(name or info.get("name") or os.path.splitext(os.path.basename(path))[0],
                     np.array(out_t)[keep], I[keep], [h for h, k in zip(out_h, keep) if k],
                     np.array(out_d)[keep], source=os.path.abspath(path),
                     meta=dict(info, spacegroup=sg.xhm(),
                               cell=[cell.a, cell.b, cell.c, cell.alpha, cell.beta, cell.gamma]),
                     hkl_int=np.array(out_hi, int).reshape(-1, 3)[keep])


def _cif_block(path):
    import gemmi
    from .readers import _decode
    text = _decode(open(path, "rb").read())
    return gemmi.cif.read_string(text).sole_block()


def cif_info(path):
    """Formula, a short display name and COD id from the CIF text (best effort)."""
    import gemmi
    out = {}
    try:
        b = _cif_block(path)
    except Exception:
        return out
    def val(*tags):
        for t in tags:
            v = b.find_value(t)
            if v and v not in ("?", "."):
                return gemmi.cif.as_string(v).strip()
        return ""
    f = val("_chemical_formula_sum")
    if f:
        out["formula"] = re.sub(r"\s+", "", f)
    nm = val("_chemical_name_mineral", "_chemical_name_common", "_chemical_name_systematic")
    nm = nm.split("\n")[0].strip()
    if nm and len(nm) <= 40:
        out["name"] = nm
    elif f:
        out["name"] = out["formula"]
    cid = val("_cod_database_code")
    if cid:
        out["cod_id"] = cid
    return out


# ---------------------------------------------------------------- from stick list
def reference_from_table(path, wavelength=1.540593, name=None):
    """Stick list text/CSV: columns 2theta (or d), intensity, optional hkl.

    Header words decide meaning: a column named d / d-spacing is treated as d.
    Works with Jade/HighScore 'peak list' exports and hand-written CSVs.
    """
    from .readers import _decode
    text = _decode(open(path, "rb").read())
    lines = [l for l in text.splitlines() if l.strip()]
    is_d = False
    rows = []
    for l in lines:
        parts = re.split(r"[\s,;\t]+", l.strip())
        nums = []
        for p in parts:
            try:
                nums.append(float(p))
            except ValueError:
                nums.append(None)
        if nums and nums[0] is None and not rows:
            hdr = [p.lower() for p in parts]
            if hdr and hdr[0].startswith("d"):
                is_d = True
            continue
        if len(nums) >= 2 and nums[0] is not None and nums[1] is not None:
            hkl = ""
            if len(parts) >= 5 and all(n is not None for n in nums[2:5]):
                hkl = "".join(str(int(v)) for v in nums[2:5])
            elif len(parts) >= 3:
                hkl = parts[2].strip("()")
            rows.append((nums[0], nums[1], hkl))
    if not rows:
        raise ValueError(f"no sticks found in {path}")
    a = np.array([[r[0], r[1]] for r in rows], float)
    if is_d or a[:, 0].max() < 15 and a[:, 0].min() > 0.5 and np.all(np.diff(a[:, 0]) <= 0):
        d = a[:, 0]
        tt = two_theta_from_d(d, wavelength)
    else:
        tt = a[:, 0]
        d = d_spacing(tt, wavelength)
    I = 100 * a[:, 1] / a[:, 1].max()
    o = np.argsort(tt)
    return Reference(name or os.path.splitext(os.path.basename(path))[0], tt[o], I[o],
                     [rows[i][2] for i in o], d[o], source=os.path.abspath(path))


def load_reference(path, wavelength=1.540593, tmin=2.0, tmax=90.0, name=None):
    if path.lower().endswith(".cif"):
        return reference_from_cif(path, wavelength, tmin, tmax, name=name)
    return reference_from_table(path, wavelength, name=name)


# ---------------------------------------------------------------- matching
def match_reference(peaks, ref: Reference, tol=0.2, tmin=None, tmax=None, n_strong=10, max_shift=0.3):
    """Compare observed peaks with a reference.

    Returns per-stick matches for the strongest reference lines in range, the fraction of those
    found, and the median 2theta offset (sample - reference).
    """
    obs = np.array([p["center"] for p in peaks]) if peaks else np.array([])
    rows = []
    idx_strong = ref.strongest(n_strong, tmin, tmax)
    # common offset (zero error / sample height) is found first, so a uniformly shifted pattern
    # still matches; the reported shifts are the raw differences
    off = 0.0
    if len(obs) and len(idx_strong):
        best = -1
        for s in np.arange(-max_shift, max_shift + 1e-9, 0.01):
            dd = np.min(np.abs(obs[None, :] - (ref.two_theta[idx_strong][:, None] + s)), axis=1)
            w = ref.intensity[idx_strong][dd <= tol / 2].sum() - 0.001 * abs(s)
            if w > best:
                best, off = w, float(s)
    for i in idx_strong:
        t = ref.two_theta[i]
        if len(obs):
            j = int(np.argmin(np.abs(obs - (t + off))))
            ok = abs(obs[j] - t - off) <= tol
        else:
            ok = False
        rows.append({"hkl": ref.hkl[i], "ref_2theta": float(t), "ref_I": float(ref.intensity[i]),
                     "found": bool(ok), "obs_2theta": float(obs[j]) if ok else np.nan,
                     "shift": float(obs[j] - t) if ok else np.nan,
                     "obs_I": peaks[j].get("rel_intensity", np.nan) if ok else np.nan,
                     "peak_index": j if ok else None})
    found = [r for r in rows if r["found"]]
    weight = sum(r["ref_I"] for r in rows) or 1
    score = sum(r["ref_I"] for r in found) / weight
    shift = float(np.median([r["shift"] for r in found])) if found else np.nan
    unexplained = []
    matched_idx = {r["peak_index"] for r in found}
    for k, p in enumerate(peaks):
        if k in matched_idx:
            continue
        if len(ref.two_theta) and np.min(np.abs(ref.two_theta - p["center"])) <= tol:
            continue
        if (tmin is None or p["center"] >= tmin) and (tmax is None or p["center"] <= tmax):
            unexplained.append(p)
    return {"reference": ref.name, "lines": rows, "score": score, "median_shift": shift,
            "shift_model": diagnose_shift(found),
            "unexplained": unexplained}


def texture_coefficients(peaks, ref: Reference, tol=0.2, n=None, tmin=None, tmax=None):
    """Harris texture coefficient TC(hkl) = (I/I0) / mean(I/I0) over matched reflections.

    TC > 1 means the plane is preferentially oriented (e.g. Zn(002) on zinc anodes).
    Uses fitted peak areas when available, otherwise heights.
    """
    m = match_reference(peaks, ref, tol, tmin, tmax, n_strong=n or len(ref.two_theta))
    rows = [r for r in m["lines"] if r["found"]]
    if len(rows) < 2:
        return []
    vals = []
    for r in rows:
        p = peaks[r["peak_index"]]
        I = p["area"] if p.get("fitted") and np.isfinite(p.get("area", np.nan)) else p["height"]
        vals.append(I / r["ref_I"])
    vals = np.array(vals)
    tc = vals / vals.mean()
    return [{"hkl": r["hkl"], "ref_2theta": r["ref_2theta"], "obs_2theta": r["obs_2theta"],
             "TC": float(t)} for r, t in zip(rows, tc)]


def diagnose_shift(found):
    """Which systematic error best explains the peak shifts (needs >= 4 matched lines)?

    zero offset      : d2theta = const
    sample height    : d2theta ~ cos(theta)   (specimen displacement)
    lattice change   : d2theta ~ tan(theta)   (uniform strain / expansion / contraction)
    """
    if len(found) < 4:
        return None
    t = np.array([r["ref_2theta"] for r in found])
    d = np.array([r["shift"] for r in found])
    th = np.radians(t / 2)
    if np.max(np.abs(d)) < 0.05:
        return {"model": "none", "note": "shifts below 0.05 deg"}
    out = {}
    for name, basis in (("zero offset", np.ones_like(th)), ("sample height", np.cos(th)),
                        ("lattice change", np.tan(th))):
        k = float(np.dot(basis, d) / np.dot(basis, basis))
        out[name] = (float(np.sqrt(np.mean((d - k * basis) ** 2))), k)
    best = min(out, key=lambda n: out[n][0])
    rms = {n: round(v[0], 4) for n, v in out.items()}
    res = {"model": best, "rms": rms}
    if best == "lattice change":
        res["strain"] = float(-out[best][1] * np.pi / 180 / 2)  # delta d / d = -cot(theta) d(theta)
        res["note"] = (f"shifts scale with tan(theta): uniform lattice {'expansion' if res['strain'] > 0 else 'contraction'}"
                       f" of about {abs(res['strain']) * 100:.2f}% (delta d/d)")
    elif best == "sample height":
        res["note"] = "shifts scale with cos(theta): most likely sample surface above/below the focusing plane"
    else:
        res["note"] = "constant shift at all angles: most likely a goniometer zero offset"
    vals = sorted(v[0] for v in out.values())
    if vals[1] < 1.3 * vals[0]:
        res["note"] += " (other models fit almost as well; cannot be distinguished with these peaks)"
    return res
