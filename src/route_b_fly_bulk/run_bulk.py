import json, os, sys, matplotlib
matplotlib.use("Agg")
import pandas as pd
pd.set_option("display.width", 220); pd.set_option("display.max_columns", 20); pd.set_option("display.max_rows", 60)
os.environ["ROUTE_BB_DATA"] = sys.argv[1]; os.environ["ROUTE_BB_GTF"] = sys.argv[1] + "/Drosophila_melanogaster.gtf.gz"
nb = json.load(open("RouteB_fly_bulk_GSE153165.ipynb"))
codes = ["".join(c["source"]) for c in nb["cells"] if c["cell_type"] == "code"]
ns = {}
for i, c in enumerate(codes):
    exec(c, ns); print(f"--- cell {i+1} ok", flush=True)
