"""Run cells 1-6 on an existing mock folder with a given TP53 mode and proliferation setting."""
import json, os, sys, matplotlib
matplotlib.use("Agg")
folder, mode, prolif = sys.argv[1], sys.argv[2], sys.argv[3] == "1"
os.environ["ROUTE_A_DATA"] = folder
nb = json.load(open("RouteA_ploidy_redox.ipynb"))
codes = ["".join(c["source"]) for c in nb["cells"] if c["cell_type"] == "code"]
codes[0] = codes[0].replace('TP53_MODE = "All lines"', f'TP53_MODE = "{mode}"').replace(
    "ADJUST_PROLIFERATION = True", f"ADJUST_PROLIFERATION = {prolif}")
ns = {}
for c in codes:
    exec(c.replace("display(", "(lambda *a: None)("), ns)
r = ns["rna"]["result"].set_index("programme")
print(f"{mode:24s} prolif={prolif!s:5s} lines={len(ns['rna']['table'])}  " + "  ".join(
    f"{k.split(': ')[1][:12]}={r.loc[k, 'partial_rho']:+.2f}" for k in ["REDOX: NRF2 targets", "REDOX: glutathione", "CONTROL: DNA damage"]))
