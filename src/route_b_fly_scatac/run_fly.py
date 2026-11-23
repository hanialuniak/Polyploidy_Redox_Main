import json, os, sys, time, matplotlib
matplotlib.use("Agg")
import pandas as pd
pd.set_option("display.width", 200); pd.set_option("display.max_columns", 20)
os.environ["ROUTE_B_DATA"] = sys.argv[1]
nb = json.load(open("RouteB_fly_ploidy_redox.ipynb"))
codes = ["".join(c["source"]) for c in nb["cells"] if c["cell_type"] == "code"]
if len(sys.argv) > 2:
    codes[0] = codes[0].replace("N_CLUSTERS = 20", f"N_CLUSTERS = {sys.argv[2]}")
codes[0] = codes[0].replace("N_RANDOM_SETS = 500", "N_RANDOM_SETS = 100").replace("N_PERMUTATIONS = 2000", "N_PERMUTATIONS = 500")
codes[0] = codes[0].replace(codes[0].split("RUNS = ")[1].split("  #@param")[0], "\"adult_brain\"")
ns = {}
for i, c in enumerate(codes):
    t = time.time(); exec(c, ns); print(f"--- cell {i+1} ok {time.time()-t:.0f}s", flush=True)
