"""Signal processing and peak analysis for powder XRD patterns."""
from __future__ import annotations

import numpy as np
from scipy.optimize import curve_fit
from scipy.signal import find_peaks, savgol_filter

CU_KA2_KA1 = 1.544414 / 1.540593


# ---------------------------------------------------------------- basics
def smooth(y, window_deg=0.0, step=0.01, order=2):
    """Savitzky-Golay smoothing; window given in degrees 2theta (0 = off)."""
    if window_deg <= 0:
        return y.copy()
    n = max(5, int(round(window_deg / step)) | 1)
    if n >= len(y):
        return y.copy()
    return savgol_filter(y, n, order)


def snip_background(y, width_deg=1.5, step=0.01, iterations=None):
    """SNIP background (Ryan et al. 1988) on LLS-transformed data. width ~ widest peak."""
    m = max(2, int(round(width_deg / step)))
    yy = np.log(np.log(np.sqrt(np.clip(y, 0, None) + 1) + 1) + 1)
    n = len(yy)
    for p in range(1, (iterations or m) + 1):
        a = np.empty_like(yy)
        a[:p] = yy[:p]
        a[n - p:] = yy[n - p:]
        a[p:n - p] = np.minimum(yy[p:n - p], 0.5 * (yy[:n - 2 * p] + yy[2 * p:]))
        yy = a
    return (np.exp(np.exp(yy) - 1) - 1) ** 2 - 1


def strip_ka2(x, y, ratio=0.5, wl1=1.540593, wl2=1.544414):
    """Rachinger correction: remove Cu Kalpha2 contribution."""
    out = y.astype(float).copy()
    th = np.radians(x / 2)
    shift = np.degrees(2 * np.arcsin(np.clip(np.sin(th) * wl2 / wl1, -1, 1))) - x
    for i in range(len(x)):
        x2 = x[i] - shift[i]
        if x2 < x[0]:
            continue
        out[i] -= ratio * np.interp(x2, x, out)
    return np.clip(out, 0, None)


def normalize(y, mode="max"):
    if mode == "max":
        r = np.nanmax(y) - np.nanmin(y)
        return (y - np.nanmin(y)) / (r if r else 1)
    if mode == "area":
        a = np.trapezoid(np.abs(y)) if hasattr(np, "trapezoid") else np.trapz(np.abs(y))
        return y / (a if a else 1)
    return y


def d_spacing(two_theta, wl=1.540593):
    return wl / (2 * np.sin(np.radians(np.asarray(two_theta) / 2)))


def two_theta_from_d(d, wl=1.540593):
    s = wl / (2 * np.asarray(d, float))
    return np.where(s < 1, np.degrees(2 * np.arcsin(np.clip(s, -1, 1))), np.nan)


def noise_level(y):
    d = np.diff(y)
    return 1.4826 * np.median(np.abs(d - np.median(d))) / np.sqrt(2)


# ---------------------------------------------------------------- peaks
def pseudo_voigt(x, h, c, w, eta, b0, b1):
    g = np.exp(-4 * np.log(2) * (x - c) ** 2 / w ** 2)
    lo = 1 / (1 + 4 * (x - c) ** 2 / w ** 2)
    return h * (eta * lo + (1 - eta) * g) + b0 + b1 * (x - c)


def ka2_shift(two_theta, wl1=1.540593, wl2=1.544414):
    """2theta(Ka2) - 2theta(Ka1) in degrees."""
    s = np.sin(np.radians(np.asarray(two_theta) / 2)) * wl2 / wl1
    return np.degrees(2 * np.arcsin(np.clip(s, -1, 1))) - np.asarray(two_theta)


def doublet(x, h, c, w, eta, b0, b1, dc=0.0, ratio=0.5):
    return pseudo_voigt(x, h, c, w, eta, b0, b1) + pseudo_voigt(x, ratio * h, c + dc, w, eta, 0, 0)


def fit_peak(x, y, center, guess_fwhm=0.15, window=None, ka2=False):
    """Pseudo-Voigt fit; with ka2=True a Cu Ka1/Ka2 doublet (fixed 2:1 ratio and splitting)."""
    dc = float(ka2_shift(center)) if ka2 else 0.0
    window = window or max(4 * guess_fwhm + dc, 0.6)
    m = (x > center - window / 2) & (x < center + window / 2 + dc)
    xs, ys = x[m], y[m]
    if len(xs) < 7:
        return None
    b0 = float(min(ys[0], ys[-1]))
    h = float(ys.max() - b0)
    p0 = [h, center, max(guess_fwhm, 2 * float(np.median(np.diff(xs)))), 0.5, b0, 0.0]
    step = float(np.median(np.diff(xs)))
    lo = [0, center - window / 4, 1.5 * step, 0, -np.inf, -np.inf]
    hi = [np.inf, center + window / 4, window, 1, np.inf, np.inf]
    model = (lambda xx, *q: doublet(xx, *q, dc=dc)) if ka2 else pseudo_voigt
    try:
        p, cov = curve_fit(model, xs, ys, p0=p0, bounds=(lo, hi), maxfev=4000)
    except Exception:
        return None
    h, c, w, eta, *_ = p
    resid = ys - model(xs, *p)
    r2 = 1 - np.sum(resid ** 2) / max(np.sum((ys - ys.mean()) ** 2), 1e-12)
    # area of pseudo-Voigt
    area = h * w * (eta * np.pi / 2 + (1 - eta) * np.sqrt(np.pi / (4 * np.log(2))))
    return {"center": float(c), "height": float(h), "fwhm": float(w), "eta": float(eta),
            "area": float(area), "r2": float(r2)}


def drop_ka2_companions(x, ys, idx, prom):
    """Remove detected maxima that are the Cu Ka2 partner of a stronger peak."""
    keep = np.ones(len(idx), bool)
    for a in range(len(idx)):
        if not keep[a]:
            continue
        t, h = x[idx[a]], ys[idx[a]]
        dc = float(ka2_shift(t))
        for b in range(a + 1, len(idx)):
            tb = x[idx[b]]
            if tb > t + 1.5 * dc:
                break
            if abs(tb - t - dc) < 0.35 * dc + 0.01 and 0.2 * h < ys[idx[b]] < 0.75 * h:
                keep[b] = False
    return idx[keep], prom[keep], int((~keep).sum())


def typical_fwhm(x, y_net, noise):
    """Median FWHM of the clearly visible peaks, found on heavily smoothed data."""
    step = float(np.median(np.diff(x)))
    ys = smooth(y_net, max(0.3, 5 * step), step)
    idx, props = find_peaks(ys, prominence=max(10 * noise, 0.1 * (ys.max() - np.median(ys))))
    if not len(idx):
        return 0.2
    _, _, l, r = __import__("scipy.signal", fromlist=["peak_widths"]).peak_widths(ys, idx, rel_height=0.5)
    w = (r - l) * step
    # heavy smoothing broadens narrow peaks: subtract the smoothing width in quadrature
    w = np.sqrt(np.clip(w ** 2 - 0.3 ** 2 / 3, (2 * step) ** 2, None))
    return float(np.median(w))


def find_xrd_peaks(x, y_net, min_snr=5.0, min_rel=0.02, min_sep_deg=0.08, fit=True, ka2=True):
    """Detect peaks on a background-subtracted pattern and refine them with pseudo-Voigt fits.

    ka2=True (Cu tube, Ka2 not stripped): Ka2 partner maxima are discarded and every peak is
    fitted as a Ka1/Ka2 doublet, so positions/FWHM refer to Ka1.
    """
    step = float(np.median(np.diff(x)))
    noise = max(noise_level(y_net), 1e-9)
    width = typical_fwhm(x, y_net, noise)
    win = float(np.clip(0.4 * width, 5 * step, 0.6))
    ys = smooth(y_net, win, step)
    # noise left after smoothing (SG window of n points reduces white noise ~ sqrt(n/2))
    n_win = max(5, int(round(win / step)))
    noise_s = noise / np.sqrt(max(1.0, n_win / 2.0))
    prom = max(min_snr * max(noise_s, noise / 4), min_rel * (ys.max() - np.median(ys)))
    sep = max(min_sep_deg, 0.5 * width)
    idx, props = find_peaks(ys, prominence=prom, distance=max(1, int(sep / step)))
    proms = props["prominences"]
    n_ka2 = 0
    if ka2 and len(idx):
        idx, proms, n_ka2 = drop_ka2_companions(x, ys, idx, proms)
    peaks = []
    for i, pr in zip(idx, proms):
        half = ys[i] - pr / 2
        l = i
        while l > 0 and ys[l] > half:
            l -= 1
        r = i
        while r < len(ys) - 1 and ys[r] > half:
            r += 1
        fw0 = max((x[r] - x[l]), 2 * step)
        pk = {"center": float(x[i]), "height": float(ys[i]), "fwhm": float(fw0), "eta": np.nan,
              "area": np.nan, "r2": np.nan, "prominence": float(pr), "fitted": False}
        if fit:
            # doublet fit only when the Ka1/Ka2 splitting is resolvable relative to the width
            use2 = ka2 and float(ka2_shift(x[i])) > 0.25 * fw0
            f = fit_peak(x, y_net, x[i], fw0, ka2=use2)
            if f and abs(f["center"] - x[i]) < fw0 and f["r2"] > 0.8 and f["fwhm"] < 5 * fw0:
                pk.update(f)
                pk["fitted"] = True
        # local significance: peak height vs noise in a +-2 deg neighbourhood, data edges excluded
        if x[i] - x[0] < 0.3 or x[-1] - x[i] < 0.3:
            continue
        loc = (x > x[i] - 2) & (x < x[i] + 2)
        ln = max(noise_level(y_net[loc]), noise * 0.5) if loc.sum() > 20 else noise
        if pk["height"] < min_snr * ln:
            continue
        pk["snr"] = float(pk["height"] / ln)
        peaks.append(pk)
    if peaks:
        hmax = max(p["height"] for p in peaks)
        for p in peaks:
            p["rel_intensity"] = 100 * p["height"] / hmax
    return peaks, noise, n_ka2


def scherrer(fwhm_deg, two_theta, wl=1.540593, K=0.9, instrument_fwhm=0.0):
    """Crystallite size in nm. Instrument broadening removed in quadrature (Gaussian approx.)."""
    b2 = np.asarray(fwhm_deg, float) ** 2 - instrument_fwhm ** 2
    b = np.radians(np.sqrt(np.clip(b2, 1e-12, None)))
    th = np.radians(np.asarray(two_theta) / 2)
    return K * wl / (b * np.cos(th)) / 10.0
