import json, os, sys, time, matplotlib
matplotlib.use("Agg")
os.environ["ROUTE_E_DATA"], os.environ["ROUTE_E_GEO"] = sys.argv[1], sys.argv[2]
nb = json.load(open("RouteE_ovary_endocycle_redox.ipynb")); ns = {}
for i, c in enumerate([c for c in nb["cells"] if c["cell_type"] == "code"]):
    t = time.time(); exec("".join(c["source"]), ns); print(f"--- cell {i+1} ok {time.time()-t:.0f}s", flush=True)
