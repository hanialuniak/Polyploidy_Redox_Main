import glob, json, os
here = os.path.dirname(os.path.abspath(__file__))
cells = [{"cell_type": "markdown", "metadata": {}, "source": "# Route F: redox in the embryonic salivary gland, stages 10-12 against 13-16\n\nThe nuclear HyPer7 signal in salivary glands is reported at stage 14. This notebook asks what changes in that tissue across exactly that window, using the embryo atlas of Peng, Jackson et al. 2024 (*Development*), which covers both windows in one study.\n\nEvery change is measured **relative to all other embryonic cells**, so anything happening across the whole embryo is removed. Choose **Runtime → Run all**."}]
for p in sorted(glob.glob(os.path.join(here, "cells", "c*.py"))):
    cells.append({"cell_type": "code", "metadata": {"cellView": "form"}, "execution_count": None, "outputs": [], "source": open(p).read().rstrip("\n")})
for c in cells:
    s = c["source"].split("\n"); c["source"] = [x + "\n" for x in s[:-1]] + [s[-1]]
json.dump({"nbformat": 4, "nbformat_minor": 0, "cells": cells, "metadata": {"kernelspec": {"name": "python3", "display_name": "Python 3"}}},
          open(os.path.join(here, "RouteF_embryo_salivary_redox.ipynb"), "w"), indent=1)
print("built", len(cells) - 1, "code cells")
