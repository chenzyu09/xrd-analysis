"""Publication-style figures."""
from __future__ import annotations

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

PALETTE = ["#1f1f1f", "#c0392b", "#2e86c1", "#28b463", "#8e44ad", "#d68910", "#17a589",
           "#a04000", "#5d6d7e", "#b03a2e", "#1a5276", "#7d6608"]
REF_COLORS = ["#555555", "#c0392b", "#2e86c1", "#28b463", "#8e44ad", "#d68910"]


def _fonts():
    """Latin font first, then a CJK font so Chinese/Japanese sample names render."""
    from matplotlib import font_manager
    have = {f.name for f in font_manager.fontManager.ttflist}
    latin = next((f for f in ["Arial", "Liberation Sans", "Helvetica", "DejaVu Sans"] if f in have), "DejaVu Sans")
    cjk = [f for f in ["Microsoft YaHei", "SimHei", "PingFang SC", "Heiti SC", "Noto Sans CJK SC",
                       "Noto Sans SC", "Source Han Sans SC", "WenQuanYi Micro Hei", "Noto Sans CJK JP"]
           if f in have]
    return [latin] + cjk[:1]


def label_y(x, y, label, cfg):
    """Height for a curve label: just above the curve's highest point under the label text."""
    xmin, xmax = cfg["x_range"]
    xr = xmin + cfg["label_x"] * (xmax - xmin)
    width = 0.024 * len(label) * (xmax - xmin) * (cfg["font_size"] / 14) * (6.0 / cfg["fig_size"][0])
    m = (x >= xr - width) & (x <= xr)
    seg = y[m] if m.any() else y[-max(1, len(y) // 20):]
    return float(np.nanmax(seg)) + 0.06 if len(seg) else 0.06


def style(base=14, frame=1.5):
    plt.rcParams.update({
        "font.family": _fonts(), "axes.unicode_minus": False, "font.size": base, "axes.linewidth": frame,
        "xtick.major.width": frame, "ytick.major.width": frame, "xtick.direction": "in",
        "ytick.direction": "in", "xtick.major.size": 5, "xtick.minor.size": 3,
        "xtick.minor.width": frame * 0.7, "mathtext.default": "regular",
        "savefig.dpi": 600, "savefig.bbox": "tight",
    })


def stacked_plot(curves, refs, cfg, out_base):
    """curves: list of dict(x, y, label, color); refs: list of Reference.

    Layout follows the common paper style: samples stacked with an offset, label above the
    right end of each curve, reference sticks (with optional hkl) in a panel at the bottom.
    """
    style(cfg["font_size"], cfg["frame_width"])
    xmin, xmax = cfg["x_range"]
    n_ref = len(refs) if cfg["reference_style"] != "none" else 0
    ref_mode = cfg["reference_style"]
    stick_panel = n_ref and ref_mode == "sticks"
    if stick_panel:
        fig, axes = plt.subplots(1 + n_ref, 1, sharex=True, figsize=cfg["fig_size"],
                                 gridspec_kw={"height_ratios": [cfg["stack_panel_ratio"]] + [1] * n_ref,
                                              "hspace": 0})
        ax = axes[0]
    else:
        fig, ax = plt.subplots(figsize=cfg["fig_size"])
        axes = [ax]
    off = 0.0
    offs = []
    bottom_profiles = []
    if n_ref and ref_mode == "profile":
        for k, r in enumerate(refs):
            xx = np.linspace(xmin, xmax, 4000)
            yy = r.profile(xx, cfg["reference_fwhm"])
            bottom_profiles.append((xx, yy * cfg["reference_scale"], r.name, REF_COLORS[k % len(REF_COLORS)]))
    for xx, yy, name, col in bottom_profiles:
        ax.plot(xx, yy + off, color=col, lw=cfg["line_width"])
        ax.text(cfg["label_x"] * (xmax - xmin) + xmin, off + 0.12 * cfg["reference_scale"],
                name, color=col, ha="right", va="bottom", fontsize=cfg["font_size"] * 0.9)
        off += cfg["reference_scale"] * cfg["offset"]
    order = list(range(len(curves)))
    if cfg["stack_order"] == "top_down":
        order = order[::-1]
    for i in order:
        c = curves[i]
        m = (c["x"] >= xmin) & (c["x"] <= xmax)
        y = c["y"][m]
        ax.plot(c["x"][m], y + off, color=c["color"], lw=cfg["line_width"])
        ytxt = off + label_y(c["x"][m], y, c["label"], cfg)
        ax.text(cfg["label_x"] * (xmax - xmin) + xmin, ytxt, c["label"], color=c["color"],
                ha="right", va="bottom", fontsize=cfg["font_size"] * 0.9)
        offs.append(off)
        top = off + (float(np.nanmax(y)) if len(y) else 1.0)
        off += cfg["offset"] * max(1.0, float(np.nanmax(y)) if len(y) else 1.0)
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(-0.05, top + 0.08 if curves else off)
    ax.set_yticks([])
    ax.set_ylabel("Intensity (a.u.)")
    if cfg["guide_lines"] and refs:
        r = refs[0]
        for t, I in zip(r.two_theta, r.intensity):
            if xmin <= t <= xmax and I >= cfg["guide_min_intensity"]:
                ax.axvline(t, color="#999999", lw=0.8, ls="--", zorder=0)
    if stick_panel:
        for k, r in enumerate(refs):
            a = axes[k + 1]
            col = REF_COLORS[k % len(REF_COLORS)]
            m = (r.two_theta >= xmin) & (r.two_theta <= xmax)
            hs = 0.55 if cfg["label_hkl"] else 0.9
            a.vlines(r.two_theta[m], 0, hs * r.intensity[m], color=col, lw=cfg["line_width"] * 1.1)
            a.set_ylim(0, 100)
            a.set_yticks([])
            a.text(0.99, 0.88, r.name, transform=a.transAxes, ha="right", va="top", color=col,
                   fontsize=cfg["font_size"] * 0.85)
            if cfg["label_hkl"]:
                # strongest first; a label is skipped if it would collide with one already placed or
                # with the phase name in the top-right corner
                min_gap = 0.018 * (xmax - xmin) * (cfg["font_size"] / 14) * (6.0 / cfg["fig_size"][0])
                name_w = 0.022 * len(r.name) * (xmax - xmin) * (cfg["font_size"] * 0.85 / 14) * (6.0 / cfg["fig_size"][0])
                name_x0 = xmin + 0.99 * (xmax - xmin) - name_w - min_gap
                placed = []
                sel = [k for k in np.argsort(-r.intensity) if m[k] and r.intensity[k] >= cfg["hkl_min_intensity"] and r.hkl[k]]
                for k in sel:
                    t = r.two_theta[k]
                    lab = "(" + r.hkl[k].split("/")[0] + ")"
                    top = hs * min(r.intensity[k], 100) + 3 + 2.2 * len(lab) * (cfg["font_size"] / 14)
                    if any(abs(p[0] - t) < min_gap for p in placed) or top > 98:
                        continue
                    if t > name_x0 and top > 60:
                        continue
                    placed.append((t, lab, r.intensity[k]))
                for t, lab, I in placed:
                    a.text(t, hs * min(I, 100) + 3, lab, rotation=90, ha="center",
                           va="bottom", fontsize=cfg["font_size"] * 0.55, color=col)
        axes[-1].set_xlabel(r"2$\theta$ (°)")
    else:
        ax.set_xlabel(r"2$\theta$ (°)")
    axes[-1].minorticks_on()
    for a in axes:
        a.tick_params(axis="y", which="both", left=False, right=False)
    for a in axes[:-1]:
        a.tick_params(axis="x", which="both", bottom=False)
    if cfg.get("title"):
        axes[0].set_title(cfg["title"])
    paths = []
    for ext in cfg["formats"]:
        p = f"{out_base}.{ext}"
        fig.savefig(p)
        paths.append(p)
    plt.close(fig)
    return paths


def peaks_plot(x, y_raw, bg, peaks, ref, name, out_base, cfg):
    """Diagnostic figure: raw data, background, detected peaks, reference positions."""
    style(11, 1.0)
    fig, ax = plt.subplots(figsize=(9, 4))
    ax.plot(x, y_raw, color="#333333", lw=0.8, label="raw")
    if bg is not None:
        ax.plot(x, bg, color="#e67e22", lw=1.0, label="background")
    for p in peaks:
        ax.plot(p["center"], np.interp(p["center"], x, y_raw), "v", color="#c0392b", ms=5)
    if ref is not None:
        top = y_raw.max()
        for t, I in zip(ref.two_theta, ref.intensity):
            if x[0] <= t <= x[-1]:
                ax.vlines(t, 0, top * 0.15 * I / 100, color="#2e86c1", lw=1)
    ax.set_xlim(x[0], x[-1])
    ax.set_xlabel(r"2$\theta$ (°)")
    ax.set_ylabel("Counts")
    ax.set_title(f"{name}: {len(peaks)} peaks" + (f" | ref {ref.name}" if ref is not None else ""))
    ax.legend(frameon=False, fontsize=9)
    p = f"{out_base}.png"
    fig.savefig(p, dpi=200)
    plt.close(fig)
    return p
