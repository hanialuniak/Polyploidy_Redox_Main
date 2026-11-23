import glob, json, os
here = os.path.dirname(os.path.abspath(__file__))
cells = [{"cell_type": "markdown", "metadata": {}, "source": "# Route E: redox across the mitotic-to-endocycle switch in the fly ovary\n\nFollicle cells divide until stage 6, then endocycle. Both are the same cell type in the same tissue, so this is a cleaner test than the salivary gland. Pick the dataset in Cell 1, then **Runtime → Run all**. Data download from GEO into *Colab Notebooks/route_e_ovary*."}]
for p in sorted(glob.glob(os.path.join(here, "cells", "c*.py"))):
    cells.append({"cell_type": "code", "metadata": {"cellView": "form"}, "execution_count": None, "outputs": [], "source": open(p).read().rstrip("\n")})
for c in cells:
    s = c["source"].split("\n"); c["source"] = [x + "\n" for x in s[:-1]] + [s[-1]]
json.dump({"nbformat": 4, "nbformat_minor": 0, "cells": cells, "metadata": {"kernelspec": {"name": "python3", "display_name": "Python 3"}}},
          open(os.path.join(here, "RouteE_ovary_endocycle_redox.ipynb"), "w"), indent=1)
print("built", len(cells) - 1, "code cells")
