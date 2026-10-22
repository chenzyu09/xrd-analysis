"""Write idealised CIFs for common phases from standard crystallographic data and check them.

Each entry: name, formula, space group (gemmi notation incl. setting), cell, sites, groups,
and the expected positions (Cu Ka1, deg 2theta) of its strongest lines from the PDF/literature,
used as a self-test. Run: python tools/build_builtin_library.py
"""
import json, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
LIB = os.path.join(HERE, "..", "xrd-analysis", "references", "library")
sys.path.insert(0, os.path.join(HERE, "..", "xrd-analysis", "scripts"))

E = [
 # metals / current collectors / catalysts
 ("Zn", "Zn", "P 63/m m c", (2.6649, 2.6649, 4.9468, 90, 90, 120), [("Zn", "Zn", 1/3, 2/3, 0.25, 1)], ["battery", "metal"], [36.30, 38.99, 43.23]),
 ("Cu", "Cu", "F m -3 m", (3.6149,)*3 + (90,)*3, [("Cu", "Cu", 0, 0, 0, 1)], ["metal", "substrate"], [43.30, 50.43, 74.13]),
 ("Al", "Al", "F m -3 m", (4.0495,)*3 + (90,)*3, [("Al", "Al", 0, 0, 0, 1)], ["metal", "substrate"], [38.47, 44.74, 65.13]),
 ("Ni", "Ni", "F m -3 m", (3.5238,)*3 + (90,)*3, [("Ni", "Ni", 0, 0, 0, 1)], ["metal", "substrate", "catalyst"], [44.51, 51.85, 76.37]),
 ("Fe (alpha)", "Fe", "I m -3 m", (2.8665,)*3 + (90,)*3, [("Fe", "Fe", 0, 0, 0, 1)], ["metal", "catalyst"], [44.67, 65.02, 82.33]),
 ("Co (hcp)", "Co", "P 63/m m c", (2.5071, 2.5071, 4.0695, 90, 90, 120), [("Co", "Co", 1/3, 2/3, 0.25, 1)], ["metal", "catalyst"], [47.4]),
 ("Ru", "Ru", "P 63/m m c", (2.7058, 2.7058, 4.2816, 90, 90, 120), [("Ru", "Ru", 1/3, 2/3, 0.25, 1)], ["metal", "catalyst"], [38.39, 42.15, 44.01]),
 ("Pt", "Pt", "F m -3 m", (3.9231,)*3 + (90,)*3, [("Pt", "Pt", 0, 0, 0, 1)], ["metal", "catalyst"], [39.76, 46.24, 67.45]),
 ("Au", "Au", "F m -3 m", (4.0786,)*3 + (90,)*3, [("Au", "Au", 0, 0, 0, 1)], ["metal", "catalyst"], [38.18, 44.39, 64.58]),
 ("Ag", "Ag", "F m -3 m", (4.0862,)*3 + (90,)*3, [("Ag", "Ag", 0, 0, 0, 1)], ["metal", "catalyst"], [38.12, 44.30, 64.44]),
 ("Ti (alpha)", "Ti", "P 63/m m c", (2.9508, 2.9508, 4.6855, 90, 90, 120), [("Ti", "Ti", 1/3, 2/3, 0.25, 1)], ["metal", "substrate"], [35.09, 38.42, 40.17]),
 ("Li", "Li", "I m -3 m", (3.5093,)*3 + (90,)*3, [("Li", "Li", 0, 0, 0, 1)], ["battery", "metal"], [36.18, 52.05]),
 ("Si", "Si", "F d -3 m :1", (5.4309,)*3 + (90,)*3, [("Si", "Si", 0, 0, 0, 1)], ["battery", "substrate"], [28.44, 47.30, 56.12]),
 ("Graphite (2H)", "C", "P 63/m m c", (2.464, 2.464, 6.711, 90, 90, 120), [("C1", "C", 0, 0, 0.25, 1), ("C2", "C", 1/3, 2/3, 0.25, 1)], ["battery", "carbon"], [26.54, 44.56]),
 # Zn battery
 ("ZnO", "ZnO", "P 63 m c", (3.2498, 3.2498, 5.2066, 90, 90, 120), [("Zn", "Zn", 1/3, 2/3, 0, 1), ("O", "O", 1/3, 2/3, 0.3819, 1)], ["battery", "oxide"], [31.77, 34.42, 36.25, 56.60]),
 ("beta-MnO2 (pyrolusite)", "MnO2", "P 42/m n m", (4.3983, 4.3983, 2.8730, 90, 90, 90), [("Mn", "Mn", 0, 0, 0, 1), ("O", "O", 0.3047, 0.3047, 0, 1)], ["battery", "oxide"], [28.68, 37.33, 56.65]),
 # Li-ion electrodes
 ("LiCoO2", "LiCoO2", "R -3 m :H", (2.8156, 2.8156, 14.0542, 90, 90, 120), [("Co", "Co", 0, 0, 0, 1), ("Li", "Li", 0, 0, 0.5, 1), ("O", "O", 0, 0, 0.2604, 1)], ["battery", "cathode"], [18.92, 37.39, 45.25]),
 ("LiNiO2 (layered, NCM-type)", "LiNiO2", "R -3 m :H", (2.8780, 2.8780, 14.1900, 90, 90, 120), [("Ni", "Ni", 0, 0, 0, 1), ("Li", "Li", 0, 0, 0.5, 1), ("O", "O", 0, 0, 0.2580, 1)], ["battery", "cathode"], [18.74, 36.65, 44.40]),
 ("LiFePO4", "LiFePO4", "P n m a", (10.3290, 6.0065, 4.6908, 90, 90, 90), [("Li", "Li", 0, 0, 0, 1), ("Fe", "Fe", 0.2822, 0.25, 0.9747, 1), ("P", "P", 0.0948, 0.25, 0.4182, 1), ("O1", "O", 0.0968, 0.25, 0.7428, 1), ("O2", "O", 0.4571, 0.25, 0.2060, 1), ("O3", "O", 0.1656, 0.0464, 0.2848, 1)], ["battery", "cathode"], [17.15, 20.78, 25.56, 29.72, 32.21, 35.59]),
 ("LiMn2O4", "LiMn2O4", "F d -3 m :2", (8.2476,)*3 + (90,)*3, [("Li", "Li", 0.125, 0.125, 0.125, 1), ("Mn", "Mn", 0.5, 0.5, 0.5, 1), ("O", "O", 0.2631, 0.2631, 0.2631, 1)], ["battery", "cathode"], [18.61, 36.09, 43.87, 63.79]),
 ("Li4Ti5O12", "Li4Ti5O12", "F d -3 m :2", (8.3595,)*3 + (90,)*3, [("Li1", "Li", 0.125, 0.125, 0.125, 1), ("Li2", "Li", 0.5, 0.5, 0.5, 1/6), ("Ti", "Ti", 0.5, 0.5, 0.5, 5/6), ("O", "O", 0.2625, 0.2625, 0.2625, 1)], ["battery", "anode"], [18.36, 35.57, 43.25, 62.83]),
 ("LiF", "LiF", "F m -3 m", (4.0270,)*3 + (90,)*3, [("Li", "Li", 0, 0, 0, 1), ("F", "F", 0.5, 0.5, 0.5, 1)], ["battery", "sei"], [38.68, 44.99, 65.52]),
 ("Li2O", "Li2O", "F m -3 m", (4.6114,)*3 + (90,)*3, [("O", "O", 0, 0, 0, 1), ("Li", "Li", 0.25, 0.25, 0.25, 1)], ["battery", "sei"], [33.59, 56.36]),
 ("NaCl", "NaCl", "F m -3 m", (5.6402,)*3 + (90,)*3, [("Na", "Na", 0, 0, 0, 1), ("Cl", "Cl", 0.5, 0.5, 0.5, 1)], ["salt"], [31.69, 45.45, 56.47]),
 # oxides common in electrocatalysis / precursors
 ("NiO", "NiO", "F m -3 m", (4.1771,)*3 + (90,)*3, [("Ni", "Ni", 0, 0, 0, 1), ("O", "O", 0.5, 0.5, 0.5, 1)], ["oxide", "catalyst"], [37.25, 43.28, 62.86]),
 ("Cu2O", "Cu2O", "P n -3 m :2", (4.2696,)*3 + (90,)*3, [("O", "O", 0.25, 0.25, 0.25, 1), ("Cu", "Cu", 0, 0, 0, 1)], ["oxide", "catalyst"], [36.42, 42.30, 61.38]),
 ("CuO (tenorite)", "CuO", "C 1 2/c 1", (4.6837, 3.4226, 5.1288, 90, 99.54, 90), [("Cu", "Cu", 0.25, 0.25, 0, 1), ("O", "O", 0, 0.4184, 0.25, 1)], ["oxide", "catalyst"], [35.54, 38.71, 48.72]),
 ("Co3O4", "Co3O4", "F d -3 m :2", (8.0837,)*3 + (90,)*3, [("Co1", "Co", 0.125, 0.125, 0.125, 1), ("Co2", "Co", 0.5, 0.5, 0.5, 1), ("O", "O", 0.2631, 0.2631, 0.2631, 1)], ["oxide", "catalyst"], [31.27, 36.85, 59.36, 65.24]),
 ("Fe3O4 (magnetite)", "Fe3O4", "F d -3 m :2", (8.3960,)*3 + (90,)*3, [("Fe1", "Fe", 0.125, 0.125, 0.125, 1), ("Fe2", "Fe", 0.5, 0.5, 0.5, 1), ("O", "O", 0.2549, 0.2549, 0.2549, 1)], ["oxide", "catalyst"], [30.10, 35.45, 43.09, 62.57]),
 ("alpha-Fe2O3 (hematite)", "Fe2O3", "R -3 c :H", (5.0356, 5.0356, 13.7489, 90, 90, 120), [("Fe", "Fe", 0, 0, 0.35530, 1), ("O", "O", 0.3059, 0, 0.25, 1)], ["oxide", "catalyst"], [24.14, 33.15, 35.61, 54.09]),
 ("TiO2 (anatase)", "TiO2", "I 41/a m d :2", (3.7852, 3.7852, 9.5139, 90, 90, 90), [("Ti", "Ti", 0, 0.75, 0.125, 1), ("O", "O", 0, 0.25, 0.0816, 1)], ["oxide"], [25.28, 37.80, 48.05]),
 ("TiO2 (rutile)", "TiO2", "P 42/m n m", (4.5937, 4.5937, 2.9587, 90, 90, 90), [("Ti", "Ti", 0, 0, 0, 1), ("O", "O", 0.3048, 0.3048, 0, 1)], ["oxide"], [27.45, 36.09, 54.32]),
 ("ZnS (sphalerite)", "ZnS", "F -4 3 m", (5.4093,)*3 + (90,)*3, [("Zn", "Zn", 0, 0, 0, 1), ("S", "S", 0.25, 0.25, 0.25, 1)], ["battery", "sulfide"], [28.56, 47.52, 56.37]),
 ("CeO2", "CeO2", "F m -3 m", (5.4113,)*3 + (90,)*3, [("Ce", "Ce", 0, 0, 0, 1), ("O", "O", 0.25, 0.25, 0.25, 1)], ["oxide", "catalyst"], [28.55, 33.08, 47.48, 56.33]),
]

CIF = """data_{tag}
_chemical_name_common '{name}'
_chemical_formula_sum '{formula}'
_xrd_analysis_source 'idealised structure from standard crystallographic data (xrd-analysis built-in library)'
_cell_length_a {a}
_cell_length_b {b}
_cell_length_c {c}
_cell_angle_alpha {al}
_cell_angle_beta {be}
_cell_angle_gamma {ga}
_space_group_name_H-M_alt '{sg}'
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_occupancy
{sites}
"""


TOP = {'Zn': 43.23, 'Cu': 43.3, 'Al': 38.47, 'Ni': 44.5, 'Fe (alpha)': 44.67, 'Co (hcp)': 47.4, 'Ru': 44.0, 'Pt': 39.76, 'Au': 38.18, 'Ag': 38.1, 'Ti (alpha)': 40.17, 'Li': 36.2, 'Si': 28.44, 'Graphite (2H)': 26.54, 'ZnO': 36.25, 'beta-MnO2 (pyrolusite)': 28.68, 'LiCoO2': 18.92, 'LiNiO2 (layered, NCM-type)': 18.74, 'LiMn2O4': 18.61, 'Li4Ti5O12': 18.36, 'LiF': 44.99, 'Li2O': 33.6, 'NaCl': 31.7, 'NiO': 43.28, 'Cu2O': 36.42, 'Co3O4': 36.85, 'Fe3O4 (magnetite)': 35.45, 'alpha-Fe2O3 (hematite)': 33.15, 'TiO2 (anatase)': 25.28, 'TiO2 (rutile)': 27.45, 'ZnS (sphalerite)': 28.56, 'CeO2': 28.55}


def main():
    from xrdtool.reference import reference_from_cif
    os.makedirs(LIB, exist_ok=True)
    index, bad = [], 0
    for name, formula, sg, cell, sites, groups, expect in E:
        tag = "".join(ch for ch in name if ch.isalnum() or ch in "-_")
        fn = f"{tag}.cif"
        txt = CIF.format(tag=tag, name=name, formula=formula, a=cell[0], b=cell[1], c=cell[2], al=cell[3],
                         be=cell[4], ga=cell[5], sg=sg,
                         sites="\n".join(f"{l} {el} {x:.5f} {y:.5f} {z:.5f} {o:.4f}" for l, el, x, y, z, o in sites))
        open(os.path.join(LIB, fn), "w", encoding="utf-8").write(txt)
        r = reference_from_cif(os.path.join(LIB, fn), tmin=5, tmax=90)
        strong = r.two_theta[r.intensity >= 10]
        miss = [t for t in expect if np.min(np.abs(strong - t)) > 0.12]
        top = r.two_theta[np.argmax(r.intensity)]
        exp_top = TOP.get(name)
        top_bad = exp_top is not None and abs(top - exp_top) > 0.15
        ok = not miss and not top_bad
        bad += not ok
        top3 = r.two_theta[np.argsort(-r.intensity)[:3]]
        print(f"{'OK ' if ok else 'BAD'} {name:28s} {r.meta['spacegroup']:14s} top3 {" ".join(f"{t:.2f}" for t in top3)}  lines>=10%: "
              + " ".join(f"{t:.2f}" for t in strong[:7]) + (f"   MISSING {miss}" if miss else "") + (f"   STRONGEST should be {exp_top}" if top_bad else ""))
        index.append({"file": fn, "name": name, "formula": formula, "groups": groups,
                      "source": "idealised (standard crystallographic data)"})
    json.dump(index, open(os.path.join(LIB, "index.json"), "w", encoding="utf-8"), indent=1, ensure_ascii=False)
    print(f"{len(index)} phases, {bad} failing")


if __name__ == "__main__":
    main()
