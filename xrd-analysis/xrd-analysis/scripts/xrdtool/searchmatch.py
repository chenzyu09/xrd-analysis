"""Search-match: rank candidate phases (CIF / stick lists) against the observed peaks.

Figure of merit per candidate (0-1), evaluated at the best constant 2theta shift within +-max_shift:
  completeness  intensity-weighted fraction of the candidate's lines (I >= 5 %, in range) that have
                an observed peak within tol            -> "are the phase's strong lines all there?"
  coverage      fraction of the observed (still unexplained) intensity those lines account for
                                                       -> "does the phase explain the pattern?"
  FOM = completeness * (0.4 + 0.6 * coverage) * accuracy, accuracy = 1 - 0.5 * mean|delta|/tol
A missing strong line costs a lot; extra observed peaks only lower coverage, so minor phases in a
mixture still rank. Phases are then accepted greedily: best phase first, its peaks marked as
explained, candidates re-scored on what is left.

Without intensity information from the sample (texture!) this is deliberately position-driven,
like the first stage of Jade/HighScore search-match. Results are candidates, not proof.
"""
from __future__ import annotations

import numpy as np


def _strained(two_theta, eps):
    """2theta after scaling all d-spacings by (1 + eps) (isotropic lattice expansion eps > 0)."""
    s = np.sin(np.radians(two_theta / 2)) / (1 + eps)
    return np.degrees(2 * np.arcsin(np.clip(s, -1, 1)))


def _metric_groups(ref, mask):
    """Split 1/d^2 of each line into contributions that scale together for the crystal system.

    cubic: 1 group; tetragonal / hexagonal / rhombohedral(hex axes): (a-terms, c-terms);
    orthorhombic: (a, b, c). Lower symmetry -> None (isotropic scaling only).
    """
    if ref.hkl_int is None or not ref.meta.get("cell"):
        return None
    a, b, c, al, be, ga = ref.meta["cell"]
    ang = np.radians([al, be, ga])
    if not np.allclose([al, be], 90, atol=0.01):
        return None
    G = np.array([[a * a, a * b * np.cos(ang[2]), a * c * np.cos(ang[1])],
                  [a * b * np.cos(ang[2]), b * b, b * c * np.cos(ang[0])],
                  [a * c * np.cos(ang[1]), b * c * np.cos(ang[0]), c * c]])
    Gs = np.linalg.inv(G)
    h = ref.hkl_int[mask].astype(float)
    terms = {(i, j): (2 - (i == j)) * Gs[i, j] * h[:, i] * h[:, j] for i in range(3) for j in range(i, 3)}
    if abs(ga - 90) < 0.01:
        if np.isclose(a, b, rtol=1e-4) and np.isclose(b, c, rtol=1e-4):
            groups = [list(terms)]
        elif np.isclose(a, b, rtol=1e-4):
            groups = [[(0, 0), (1, 1), (0, 1)], [(2, 2), (0, 2), (1, 2)]]
        else:
            groups = [[(0, 0)], [(1, 1)], [(2, 2)], [(0, 1), (0, 2), (1, 2)]]
            groups = groups[:3] if np.allclose([terms[k] for k in groups[3]], 0) else None
            if groups is None:
                return None
    elif abs(ga - 120) < 0.01 and np.isclose(a, b, rtol=1e-4):
        groups = [[(0, 0), (1, 1), (0, 1)], [(2, 2), (0, 2), (1, 2)]]
    else:
        return None
    return np.stack([sum(terms[k] for k in g) for g in groups], axis=1)


def _evaluate(rt, ri, obs_t, obs_w, active, tot_obs, tol):
    d = np.abs(obs_t[None, :] - rt[:, None])  # lines x peaks
    j = np.argmin(d, axis=1)
    dmin = d[np.arange(len(rt)), j]
    hit = dmin <= tol
    comp = ri[hit].sum() / ri.sum()
    used = np.zeros(len(obs_t), bool)
    used[j[hit]] = True
    cov = obs_w[used & active].sum() / tot_obs
    acc = 1 - 0.5 * (dmin[hit].mean() / tol) if hit.any() else 0
    return comp * (0.4 + 0.6 * cov) * acc, comp, cov, hit, used


def _refine_cell(Gm, ri, rt_start, obs_t, wl, shift, max_strain=0.03):
    """Least-squares lattice refinement on matched lines (weights sqrt(I)), tolerances 0.35->0.15 deg."""
    q_ref = Gm.sum(axis=1)
    lo, hi = 1 / (1 + max_strain) ** 2, 1 / (1 - max_strain) ** 2
    rt = rt_start.copy()
    s = np.ones(Gm.shape[1])
    for tol in (0.35, 0.25, 0.15):
        d = np.abs(obs_t[None, :] - rt[:, None])
        j = np.argmin(d, axis=1)
        ok = d[np.arange(len(rt)), j] <= tol
        if ok.sum() < Gm.shape[1] + 1:
            return None
        th = np.radians((obs_t[j[ok]] - shift) / 2)
        q_obs = (2 * np.sin(th) / wl) ** 2
        w = np.sqrt(ri[ok])
        s, *_ = np.linalg.lstsq(Gm[ok] * w[:, None], q_obs * w, rcond=None)
        s = np.clip(s, lo, hi)
        q = Gm @ s
        rt = np.degrees(2 * np.arcsin(np.clip(wl * np.sqrt(q) / 2, -1, 1))) + shift
    return rt, 1 / np.sqrt(s) - 1  # refined positions, strain per axis group


def _score(ref, obs_t, obs_w, active, tol, max_shift, tmin, tmax, min_I=5.0, max_strain=0.02):
    """Best FOM over a constant 2theta offset (zero / height error) and an isotropic lattice scale
    (doping, solid solution, lithiation, temperature, DFT-relaxed reference cells), followed by a
    symmetry-constrained cell refinement (a and c separately for tetragonal/hexagonal, etc.)."""
    m = (ref.two_theta >= tmin) & (ref.two_theta <= tmax) & (ref.intensity >= min_I)
    if m.sum() == 0:
        return None
    rt0, ri = ref.two_theta[m], ref.intensity[m]
    tot_obs = obs_w[active].sum() or 1.0
    best = None
    strains = np.arange(-max_strain, max_strain + 1e-12, 0.002) if max_strain > 0 else [0.0]
    shifts = np.arange(-max_shift, max_shift + 1e-9, 0.03)
    for eps in strains:
        rt_e = _strained(rt0, eps)
        for shift in shifts:
            fom, comp, cov, hit, used = _evaluate(rt_e + shift, ri, obs_t, obs_w, active, tot_obs, tol)
            # prefer the simplest explanation: small offsets and small lattice changes
            fom *= (1 - 0.1 * abs(shift) / max(max_shift, 1e-9)) * (1 - 0.15 * abs(eps) / max(max_strain, 1e-9))
            if best is None or fom > best["fom"]:
                best = {"fom": float(fom), "completeness": float(comp), "coverage": float(cov),
                        "shift": float(shift), "strain": float(eps), "strain_axes": None,
                        "hit": hit, "used": used, "rt": rt_e + shift}
    # anisotropic refinement
    Gm = _metric_groups(ref, m)
    if Gm is not None and Gm.shape[1] > 1 and best["completeness"] >= 0.3:
        wl = float(np.median(2 * ref.d[m] * np.sin(np.radians(rt0 / 2))))
        r = _refine_cell(Gm, ri, best["rt"], obs_t, wl, best["shift"], max_strain=0.03)
        if r is not None:
            rt, st = r
            fom, comp, cov, hit, used = _evaluate(rt, ri, obs_t, obs_w, active, tot_obs, tol)
            fom *= (1 - 0.1 * abs(best["shift"]) / max(max_shift, 1e-9)) * (1 - 0.15 * min(1, np.max(np.abs(st)) / 0.03))
            if fom > best["fom"]:
                best.update(fom=float(fom), completeness=float(comp), coverage=float(cov),
                            strain=float(np.mean(st)), strain_axes=[float(v) for v in st],
                            hit=hit, used=used, rt=rt)
    hit = best.pop("hit")
    used = best.pop("used")
    best.pop("rt")
    best.update(n_lines=int(m.sum()), n_hit=int(hit.sum()), explains=np.where(used)[0].tolist(),
                strong_missing=[float(t) for t, i, h in zip(rt0, ri, hit) if not h and i >= 30],
                strongest_hit=bool(hit[np.argmax(ri)]))
    return best


def search_match(peaks, candidates, tmin, tmax, tol=0.15, max_shift=0.3, max_phases=3,
                 min_fom=0.45, top=15):
    """Rank candidates and pick a phase assemblage.

    peaks: list of dicts with 'center' and 'height' (from find_xrd_peaks)
    candidates: list of Reference
    Returns {'ranking': [...], 'accepted': [...], 'unexplained': [2theta,...]}
    """
    if not peaks or not candidates:
        return {"ranking": [], "accepted": [], "unexplained": [p["center"] for p in peaks or []]}
    obs_t = np.array([p["center"] for p in peaks])
    obs_w = np.array([max(p["height"], 0) for p in peaks], float)
    obs_w = np.sqrt(obs_w)  # compress dynamic range: one huge textured peak must not dominate
    active = np.ones(len(obs_t), bool)
    ranking = []
    for ref in candidates:
        s = _score(ref, obs_t, obs_w, active, tol, max_shift, tmin, tmax)
        if s:
            ranking.append(dict(s, ref=ref))
    ranking.sort(key=lambda r: -r["fom"])
    accepted = []
    pool = list(ranking)
    for _ in range(max_phases):
        if not active.any() or not pool:
            break
        rescored = []
        for r in pool:
            s = _score(r["ref"], obs_t, obs_w, active, tol, max_shift, tmin, tmax)
            if s:
                rescored.append(dict(s, ref=r["ref"]))
        rescored.sort(key=lambda r: -r["fom"])
        if not rescored:
            break
        b = None
        for cand in rescored:
            new = [k for k in cand["explains"] if active[k]]
            if not accepted:
                ok = cand["fom"] >= min_fom
            else:
                # a further phase: its strong lines must be present, it must explain peaks the
                # accepted phases do not, and its strongest line must be among them or present
                new_share = obs_w[new].sum() / (obs_w.sum() or 1)
                ok = (cand["completeness"] >= 0.7 and cand["coverage"] >= 0.12 and cand["fom"] >= 0.35
                      and cand["strongest_hit"] and len(new) >= 2 and new_share >= 0.05)
            if ok:
                b = cand
                break
            if not accepted:
                break
        if b is None:
            break
        new = [k for k in b["explains"] if active[k]]
        accepted.append(b)
        active[new] = False
        pool = [r for r in rescored[1:] if not _same_phase(r["ref"], b["ref"])]
    return {"ranking": ranking[:top], "accepted": accepted,
            "unexplained": [float(t) for t, a in zip(obs_t, active) if a]}


def _same_phase(a, b):
    fa, fb = a.meta.get("formula"), b.meta.get("formula")
    return bool(fa) and fa == fb and a.meta.get("spacegroup") == b.meta.get("spacegroup")


def describe(r):
    ref = r["ref"]
    meta = ref.meta
    parts = [ref.name]
    if meta.get("formula") and meta.get("formula") != ref.name:
        parts.append(meta["formula"])
    if meta.get("spacegroup"):
        parts.append(meta["spacegroup"])
    if meta.get("cod_id"):
        parts.append(f"COD {meta['cod_id']}")
    return ", ".join(parts)
