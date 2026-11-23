import json, os, sys, time, matplotlib
matplotlib.use("Agg")
os.environ["ROUTE_D_DATA"], os.environ["ROUTE_D_GEO"], os.environ["ROUTE_D_GTF"], gse = sys.argv[1:5]
nb = json.load(open("RouteD_salivary_endocycle_redox.ipynb")); ns = {}
for i, c in enumerate([c for c in nb["cells"] if c["cell_type"] == "code"]):
    src = "".join(c["source"])
    if i == 0: src = src.replace('GSE = "GSE121956"', f'GSE = "{gse}"')
    t = time.time(); exec(src, ns); print(f"--- cell {i+1} ok {time.time()-t:.0f}s", flush=True)
