import json, os
here = os.path.dirname(os.path.abspath(__file__))
intro = ("# Route B (fly): download single-cell ATAC files from the Aerts lab\n\n"
         "Run Cells 1 and 2 to see what is available. Then type the runs you want into **RUNS** in Cell 1, "
         "run Cell 1 again, and run Cell 3. The files go into the folder the main notebook reads.")
cells = [{"cell_type": "markdown", "metadata": {}, "source": intro}]
for n in ("c1.py", "c2.py", "c3.py"):
    cells.append({"cell_type": "code", "metadata": {"cellView": "form"}, "execution_count": None, "outputs": [],
                  "source": open(os.path.join(here, n)).read().rstrip("\n")})
for c in cells:
    parts = c["source"].split("\n"); c["source"] = [s + "\n" for s in parts[:-1]] + [parts[-1]]
json.dump({"nbformat": 4, "nbformat_minor": 0, "cells": cells,
           "metadata": {"colab": {"provenance": []}, "kernelspec": {"name": "python3", "display_name": "Python 3"}}},
          open(os.path.join(here, "RouteB_fly_download_Aerts.ipynb"), "w"), indent=1)
print("built")
