import json, os, sys, types
here = os.path.dirname(os.path.abspath(__file__))
nb = json.load(open(os.path.join(here, "RouteF_embryo_salivary_redox.ipynb")))
g = {"__name__": "__main__"}
for i, c in enumerate([c for c in nb["cells"] if c["cell_type"] == "code"], start=1):
    src = "".join(c["source"])
    print(f"\n{'=' * 24} CELL {i} {'=' * 24}")
    try:
        exec(compile(src, f"cell{i}", "exec"), g)
    except Exception as e:
        print(f"!! CELL {i} FAILED: {type(e).__name__}: {e}")
        import traceback; traceback.print_exc()
        sys.exit(1)
print("\nAll cells ran.")
