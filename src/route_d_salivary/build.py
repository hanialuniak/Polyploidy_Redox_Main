import glob, json, os
here = os.path.dirname(os.path.abspath(__file__))
cells = [{"cell_type": "markdown", "metadata": {}, "source": "# Route D: redox in endocycling fly cells (salivary glands GSE121956; induced endocycles in S2 cells GSE121955)\n\nPick the dataset in Cell 1, then **Runtime → Run all**. Data download from GEO into *Colab Notebooks/route_d_salivary*. Run once per dataset."}]
for p in sorted(glob.glob(os.path.join(here, "cells", "c*.py"))):
    cells.append({"cell_type": "code", "metadata": {"cellView": "form"}, "execution_count": None, "outputs": [], "source": open(p).read().rstrip("\n")})
for c in cells:
    s = c["source"].split("\n"); c["source"] = [x + "\n" for x in s[:-1]] + [s[-1]]
json.dump({"nbformat": 4, "nbformat_minor": 0, "cells": cells, "metadata": {"kernelspec": {"name": "python3", "display_name": "Python 3"}}},
          open(os.path.join(here, "RouteD_salivary_endocycle_redox.ipynb"), "w"), indent=1)
print("built")
