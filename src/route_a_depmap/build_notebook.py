"""Builds RouteA_ploidy_redox.ipynb from cells/c*.py."""
import glob, json, os
here = os.path.dirname(os.path.abspath(__file__))
intro = ("# Route A: an exploratory screen of ploidy, redox and Geminin in cancer cell lines\n\n"
         "Choose **Runtime \u2192 Run all**. The step-by-step guide (PDF) explains every cell and screen.")
cells = [{"cell_type": "markdown", "metadata": {}, "source": intro}]
for path in sorted(glob.glob(os.path.join(here, "cells", "c*.py"))):
    cells.append({"cell_type": "code", "metadata": {"cellView": "form"}, "execution_count": None, "outputs": [],
                  "source": open(path).read().rstrip("\n")})
for cell in cells:
    parts = cell["source"].split("\n")
    cell["source"] = [p + "\n" for p in parts[:-1]] + [parts[-1]]
notebook = {"nbformat": 4, "nbformat_minor": 0, "cells": cells,
            "metadata": {"colab": {"provenance": []}, "kernelspec": {"name": "python3", "display_name": "Python 3"},
                         "language_info": {"name": "python"}}}
json.dump(notebook, open(os.path.join(here, "RouteA_ploidy_redox.ipynb"), "w"), indent=1)
print(len(cells) - 1, "code cells written")
