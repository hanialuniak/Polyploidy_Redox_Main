import glob, json, os
here = os.path.dirname(os.path.abspath(__file__))
cells = [{"cell_type": "markdown", "metadata": {}, "source": "# Route C (rat): redox enzyme genes in hepatocytes sorted by ploidy (GSE132409, Katsuda et al. 2020)\n\nChoose **Runtime → Run all**. The data download from GEO automatically into *Colab Notebooks/route_c_rat*."}]
for p in sorted(glob.glob(os.path.join(here, "cells", "c*.py"))):
    cells.append({"cell_type": "code", "metadata": {"cellView": "form"}, "execution_count": None, "outputs": [], "source": open(p).read().rstrip("\n")})
for c in cells:
    s = c["source"].split("\n"); c["source"] = [x + "\n" for x in s[:-1]] + [s[-1]]
json.dump({"nbformat": 4, "nbformat_minor": 0, "cells": cells, "metadata": {"kernelspec": {"name": "python3", "display_name": "Python 3"}}},
          open(os.path.join(here, "RouteC_rat_GSE132409.ipynb"), "w"), indent=1)
print("built")
