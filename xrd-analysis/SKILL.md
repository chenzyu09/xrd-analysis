---
name: xrd-analysis
description: Analyse powder XRD data files (Rigaku .txt/.ras/.rasx/.raw, Bruker .raw, .xy, .xlsx) - peak finding, phase matching against CIF or stick-list references, crystallite size, texture coefficients, intensity ratios such as I(003)/I(104) - and produce publication-style stacked plots, an Origin script and Excel tables. Use when the user shares XRD files or asks to plot, compare or interpret XRD patterns.
license: MIT
---

# XRD analysis

All numbers come from the Python tool in `scripts/`. Do not estimate peak positions, widths or
ratios yourself; run the tool and read its outputs. Your job is to choose the options, run it,
check the figures, and explain the results.

## 1. Run

Requirements: Python 3.9+, `pip install -r scripts/requirements.txt` (numpy, scipy, matplotlib,
openpyxl; `gemmi` only for CIF references).

```bash
python scripts/xrd.py <files or folders> [options]
```

| Option | Use |
|---|---|
| `--names A B C` | Display names, same order as the files |
| `--ref phase.cif` (repeatable) `--ref-label "LiCoO2"` | Reference: CIF (pattern calculated) or stick list (2theta/d, I, hkl) |
| `--range 10 80` | 2theta plot range |
| `--ratio 003/104` | Intensity ratio by hkl of the first reference, or by 2theta (`18.7/44.5`) |
| `--style sticks / profile / none` | Reference as stick panel (default) or as a calculated pattern curve |
| `--subtract-bg`, `--smooth 0.05`, `--ka2` | Plot background-subtracted, smoothed, or Ka2-stripped data |
| `--instrument-fwhm 0.06` | Instrumental broadening for Scherrer correction (from a Si/LaB6 standard) |
| `--search Zn,O,S,H` | Identify phases from COD online: downloads every structure made of these elements (or a subset), calculates patterns, ranks them against the peaks, adds the accepted phases to the plot. `--exact-elements`, `--search-text MIL-101` narrow it |
| `--cod 1234567` (repeatable) | Use a COD entry as reference (downloaded and cached) |
| `--library [groups]` | Identify phases from the built-in library (34 common phases: metals/current collectors, Li/Zn battery materials, catalyst oxides); `--list-library` shows them |
| `--no-search` | Skip phase identification |
| `--config out/config.json` | Re-run after editing labels, colours, order, fonts, offsets |
| `--calc-ref phase.cif` | Print the calculated stick list of a CIF and exit |

A folder input picks one file per measurement (text exports are preferred over binary duplicates).

Outputs (default `xrd_output/` next to the first data file):
- `xrd_stacked.png/.pdf` - publication figure
- `diagnostics/<sample>.png` - raw data, background, detected peaks, reference positions
- `report.md` - checks, phase matching, texture, ratios, peak lists
- `summary.json` - the same, machine readable (read this to answer questions)
- `xrd_data.xlsx` - plot data (3 header rows: long name / unit / comment, Origin-ready), peaks, matches
- `origin_plot.py` - Origin 2021+ embedded-Python script with the data inside
- `config.json` - every plotting/analysis setting; edit and re-run with `--config`

## 2. Workflow

1. Ask (only if unknown): what the samples are and which elements they contain. References, in order
   of preference: the user's own CIF / card peak list; `--cod ID` if they know the entry; `--search
   <elements>` (needs internet); otherwise the built-in library is searched automatically when no
   reference is given. Never invent a reference pattern from memory.
2. Run the tool. Look at `xrd_stacked.png` and every `diagnostics/*.png` before reporting.
3. Fix presentation on request through `config.json` (labels, colours, `stack_order`,
   `intensity_scale: sqrt` for very strong textured peaks, `x_range`, `font_size`, `fig_size`).
4. Report per sample: phase match (score, missing strong lines, unexplained peaks), shifts and their
   diagnosed cause, crystallite size (with caveats), texture/ratios, and every item under "Checks".

## 3. Interpreting the outputs

- **Search-match FOM** (0-1): completeness of the phase's strong lines x share of the observed intensity
  it explains x position accuracy, after a common shift of up to +-0.3 deg. The report lists the
  proposed assemblage, a ranking table and the peaks left unexplained. It is position-based (textured
  samples are fine); treat it as candidates. Chemistry must agree: reject a phase whose elements are not
  in the sample. "No confident match" with the built-in library usually means the phase is not in it:
  suggest `--search <elements>`.

- **Match score** = intensity-weighted fraction of the strongest reference lines found within
  `match_tolerance` (0.2 deg). >0.8 good, 0.5-0.8 partial (texture, overlap, or minor phase), <0.5 doubtful.
  Missing a strong line matters more than extra weak peaks.
- **Unassigned peaks**: second phase, impurity, precursor, substrate. The tool flags lines of Cu/Al/Ni/Ti/Zn
  foils, carbon paper, stainless steel and Si. Suggest candidates but say they need a reference to confirm.
- **Shift model** (needs >= 4 matched lines): constant shift = zero offset; ~cos(theta) = sample height
  error (common for foils and electrodes); ~tan(theta) = real lattice change (report delta d/d).
  Only call it lattice expansion/contraction when the tan(theta) model wins.
- **Crystallite size** (Scherrer, K = 0.9): lower bound unless `instrument_fwhm` is set; ignores strain;
  meaningless above ~100-150 nm or when the report says the step is too coarse. Compare samples relatively.
- **Texture coefficient TC(hkl)** (Harris): 1 = random, >1 preferred. For Zn anodes TC(002) and I(002)/I(101)
  are the usual descriptors of (002)-textured deposition. Rolled foils are already textured; compare with the
  pristine foil, not with the powder reference alone.
- **Layered cathodes** (R-3m, LiCoO2/NCM/LNO): `--ratio 003/104`; >1.2 usually read as low Li/Ni mixing.
  Clear splitting of (006)/(012) and (018)/(110) indicates good layering.
- **Broad hump** in the background: amorphous content, glass holder, polymer binder/separator, or
  low-crystallinity carbon (~20-26 deg).
- **High background fraction** with Fe/Co/Mn samples under Cu radiation: fluorescence, not amorphous phase.
- **MOFs**: low-angle peaks (2-10 deg) carry the framework information; check sharp extra peaks against
  the linker (e.g. terephthalic acid ~17.3, 25.1, 27.9 deg) and metal salts.

Keep conclusions qualitative where the method is qualitative. Phase fractions, lattice parameters and
site occupancies need Rietveld refinement (GSAS-II, FullProf, TOPAS); say so if asked.

## 4. Getting references

- CIF: COD (crystallography.net, free), Materials Project, ICSD/CCDC (licensed), or the paper's SI.
- Stick list: export the PDF card peak list from Jade/HighScore as text (2theta or d, I, h k l) and pass it
  with `--ref`; the header word `d` marks d-spacing columns.
- COD: `--search` / `--cod` download automatically (cache in `~/.xrd-analysis/cod`). If COD is unreachable
  the tool says so and falls back to the built-in library.
- `references/library/` holds the built-in idealised structures (positions and strongest lines are
  self-tested by `tools/build_builtin_library.py`); intensities are calculated, not measured.
