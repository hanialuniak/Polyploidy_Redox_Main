import glob, json, os
here = os.path.dirname(os.path.abspath(__file__))
intro = ("# Route B (fly): polyploidy and redox in the Drosophila brain, from single-cell ATAC fragments\n\n"
         "Choose **Runtime → Run all**. The step-by-step guide (PDF) explains every cell and screen.")
cells = [{"cell_type": "markdown", "metadata": {}, "source": intro}]
for p in sorted(glob.glob(os.path.join(here, "cells", "c*.py"))):
    cells.append({"cell_type": "code", "metadata": {"cellView": "form"}, "execution_count": None, "outputs": [],
                  "source": open(p).read().rstrip("\n")})
for c in cells:
    parts = c["source"].split("\n"); c["source"] = [s + "\n" for s in parts[:-1]] + [parts[-1]]
json.dump({"nbformat": 4, "nbformat_minor": 0, "cells": cells,
           "metadata": {"colab": {"provenance": []}, "kernelspec": {"name": "python3", "display_name": "Python 3"},
                        "language_info": {"name": "python"}}}, open(os.path.join(here, "RouteB_fly_ploidy_redox.ipynb"), "w"), indent=1)
print(len(cells) - 1, "code cells")
