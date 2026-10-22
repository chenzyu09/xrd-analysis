#!/usr/bin/env python3
"""xrd-analysis command line.

Examples
  python xrd.py sample1.txt sample2.rasx --ref Zn.cif
  python xrd.py data_folder/ --ref LiCoO2.cif --ratio 003/104 --range 10 80
  python xrd.py --config xrd_output/config.json          # re-plot after editing the config
  python xrd.py --calc-ref Zn.cif                          # print the calculated stick list
"""
import argparse
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from xrdtool.pipeline import merge_config, run  # noqa: E402
from xrdtool.readers import find_patterns  # noqa: E402


def main(argv=None):
    ap = argparse.ArgumentParser(description="Powder XRD analysis and publication plots")
    ap.add_argument("inputs", nargs="*", help="data files or folders")
    ap.add_argument("--ref", action="append", default=[], help="reference: .cif or stick list (repeatable)")
    ap.add_argument("--ref-label", action="append", default=[], help="display name for each --ref")
    ap.add_argument("--names", nargs="+", help="display names for the samples, in input order")
    ap.add_argument("--range", nargs=2, type=float, metavar=("MIN", "MAX"), help="2theta plot range")
    ap.add_argument("--ratio", action="append", default=[], help='intensity ratio, e.g. "003/104" or "18.7/44.5"')
    ap.add_argument("--style", choices=["sticks", "profile", "none"], help="how references are drawn")
    ap.add_argument("--subtract-bg", action="store_true", help="plot background-subtracted curves")
    ap.add_argument("--ka2", action="store_true", help="strip Cu K-alpha2 before analysis")
    ap.add_argument("--smooth", type=float, help="smoothing window for the plot, deg")
    ap.add_argument("--instrument-fwhm", type=float, help="instrumental FWHM (deg) for Scherrer correction")
    ap.add_argument("--search", metavar="ELEMENTS",
                    help="identify phases from COD (online) using these elements, e.g. Zn,O,S,H")
    ap.add_argument("--exact-elements", action="store_true",
                    help="with --search: only phases containing all the given elements (default: any subset)")
    ap.add_argument("--search-text", default="", help="with --search: extra keyword, e.g. a mineral or MOF name")
    ap.add_argument("--cod", action="append", default=[], metavar="ID", help="use COD entry as reference (repeatable)")
    ap.add_argument("--library", nargs="*", metavar="GROUP",
                    help="identify phases from the built-in library (optionally only groups: battery cathode "
                         "anode metal oxide catalyst substrate carbon sei salt sulfide)")
    ap.add_argument("--no-search", action="store_true", help="do not run automatic phase identification")
    ap.add_argument("--list-library", action="store_true", help="list the built-in phases and exit")
    ap.add_argument("--config", help="config.json from a previous run")
    ap.add_argument("--out", default=None, help="output folder (default: xrd_output next to the data)")
    ap.add_argument("--calc-ref", help="only print the calculated pattern of this CIF and exit")
    a = ap.parse_args(argv)

    if a.list_library:
        from xrdtool.candidates import library_index
        for e in library_index():
            print(f"{e['name']:30s} {e['formula']:12s} {', '.join(e['groups'])}")
        return 0

    if a.calc_ref:
        from xrdtool.reference import reference_from_cif
        r = reference_from_cif(a.calc_ref, tmin=(a.range or [5, 90])[0], tmax=(a.range or [5, 90])[1])
        print(f"# {r.name}  {r.meta}")
        print("2theta\tI\td\thkl")
        for t, i, d, h in zip(r.two_theta, r.intensity, r.d, r.hkl):
            print(f"{t:.3f}\t{i:.1f}\t{d:.4f}\t{h}")
        return 0

    user = {}
    if a.config:
        with open(a.config, encoding="utf-8") as f:
            user = json.load(f)
    cfg = merge_config(user)
    files = []
    if a.inputs:
        cfg["samples"] = [] if not a.config else cfg["samples"]
        for p in a.inputs:
            files += find_patterns(p) if os.path.isdir(p) else [p]
    elif cfg["samples"]:
        files = [s["file"] for s in cfg["samples"]]
    if not files:
        ap.error("no input data (give files/folders or --config)")
    if a.ref:
        cfg["references"] = [{"file": os.path.abspath(r),
                              "label": a.ref_label[i] if i < len(a.ref_label) else None, "show": True}
                             for i, r in enumerate(a.ref)]
    if a.range:
        cfg["x_range"] = list(a.range)
    sc = cfg.setdefault("search", {})
    if a.cod:
        sc["cod_ids"] = a.cod
    if a.search:
        sc.update(source="cod", elements=[e for e in a.search.replace(" ", ",").split(",") if e],
                  subsets=not a.exact_elements, text=a.search_text, add_to_plot=True)
    elif a.library is not None:
        sc.update(source="library", library_groups=a.library, add_to_plot=True)
    if a.no_search:
        sc["source"] = "none"
    if a.ratio:
        cfg["ratios"] = a.ratio
    if a.style:
        cfg["reference_style"] = a.style
    if a.subtract_bg:
        cfg["subtract_background_in_plot"] = True
    if a.ka2:
        cfg["strip_ka2"] = True
    if a.smooth is not None:
        cfg["smooth"] = a.smooth
    if a.instrument_fwhm is not None:
        cfg["instrument_fwhm"] = a.instrument_fwhm
    out = a.out or (os.path.dirname(os.path.abspath(a.config)) if a.config else
                    os.path.join(os.path.dirname(os.path.abspath(files[0])), "xrd_output"))
    res = run(files, cfg, out, a.names)
    print(f"Analysed {len(files)} pattern(s):")
    for r in res["results"]:
        flag = f"  [{len(r['checks'])} check(s)]" if r["checks"] else ""
        best = ""
        if r.get("search") is not None:
            acc = r["search"]["accepted"]
            best = (", phases: " + " + ".join(f"{x['ref'].name} (FOM {x['fom']:.2f})" for x in acc)
                    if acc else ", phases: no confident match")
        elif r["matches"]:
            m = max(r["matches"], key=lambda m: m["score"])
            best = f", match {m['reference']} score {m['score']:.2f}"
        print(f"  {r['name']}: {len(r['peaks'])} peaks{best}{flag}")
    print(f"Outputs in {res['out_dir']}:")
    for f in sorted(os.listdir(res["out_dir"])):
        print("  ", f)
    print("Edit config.json and run:  python xrd.py --config", res["config"])
    return 0


if __name__ == "__main__":
    sys.exit(main())
