import json, os, sys, time, matplotlib
matplotlib.use("Agg")
import pandas as pd
pd.set_option("display.width", 220); pd.set_option("display.max_columns", 20)
os.environ["ROUTE_CR_DATA"] = sys.argv[1]; os.environ["ROUTE_CR_GEO"] = sys.argv[2]
nb = json.load(open("RouteC_rat_GSE132409.ipynb")); ns = {}
for i, c in enumerate([c for c in nb["cells"] if c["cell_type"] == "code"]):
    t = time.time(); exec("".join(c["source"]), ns); print(f"--- cell {i+1} ok {time.time()-t:.0f}s", flush=True)
