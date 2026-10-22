"""Benchmark search-match on the measured patterns published with XRD-AutoAnalyzer
(Szymanski et al., npj Comput. Mater. 7, 2021; https://github.com/njszym/XRD-AutoAnalyzer).

The true phases are encoded in each spectrum's file name (e.g. Mn2O3+LiF.xy); candidates are the
56 reference CIFs of that repository. Usage:
  git clone https://github.com/njszym/XRD-AutoAnalyzer
  python tools/benchmark_autoanalyzer.py XRD-AutoAnalyzer/Example
"""
import glob, os, re, sys, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "xrd-analysis", "scripts"))
from xrdtool.readers import read_pattern
from xrdtool.processing import snip_background, find_xrd_peaks
from xrdtool.reference import reference_from_cif
from xrdtool.searchmatch import search_match
import numpy as np


def norm_formula(f):
    return re.sub(r"[()\s]", "", f)


def main(root):
    cands = []
    for p in sorted(glob.glob(os.path.join(root, "References", "*.cif"))):
        stem = os.path.splitext(os.path.basename(p))[0]
        try:
            r = reference_from_cif(p, tmin=5, tmax=90, name=stem)
            r.meta["formula"] = stem
            cands.append(r)
        except Exception as e:  # noqa: BLE001
            print("skip", stem, e)
    print(f"{len(cands)} candidate phases")
    tot = hit = fp = 0
    rows = []
    for sp in sorted(glob.glob(os.path.join(root, "Spectra", "*.xy"))):
        truth = [norm_formula(t) for t in os.path.splitext(os.path.basename(sp))[0].split("+")]
        pat = read_pattern(sp)
        y = pat.intensity.astype(float)
        bg = snip_background(y, 1.5, pat.step)
        peaks, noise, _ = find_xrd_peaks(pat.two_theta, np.clip(y - bg, 0, None), ka2=True)
        t0 = time.time()
        sm = search_match(peaks, cands, pat.two_theta[0], pat.two_theta[-1])
        found = [a["ref"].name for a in sm["accepted"]]
        found_f = [norm_formula(f.split("_")[0]) for f in found]
        ok = [t for t in truth if t in found_f]
        extra = [f for f, ff in zip(found, found_f) if ff not in truth]
        tot += len(truth); hit += len(ok); fp += len(extra)
        top = [r["ref"].name for r in sm["ranking"][:5]]
        rows.append((os.path.basename(sp), truth, found, ok, extra, top))
        print(f"{os.path.basename(sp):26s} truth={'+'.join(truth):22s} found={' + '.join(found) or '-':40s} "
              f"{'OK' if len(ok) == len(truth) and not extra else 'PARTIAL' if ok else 'MISS'}"
              f"{'  extra: ' + ', '.join(extra) if extra else ''}")
        if len(ok) < len(truth):
            print(f"{'':26s} top5: {', '.join(top)}")
    print(f"\nphases found {hit}/{tot}; wrong extra phases {fp}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "XRD-AutoAnalyzer/Example")
