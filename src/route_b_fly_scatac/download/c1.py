#@title Cell 1 - Settings and connect Google Drive
import os, re, subprocess, urllib.parse
import pandas as pd
try:
    from IPython.display import display
except ImportError:
    display = print

BASE_URL = "https://flybrain.aertslab.org/downloads/cellRanger/"  #@param {type:"string"}
DATA_DIR = "/content/drive/MyDrive/Colab Notebooks/route_b_fly"  #@param {type:"string"}
#@markdown **RUNS** is pre-filled with the two larval brain runs (about 660 MB), which carry the known polyploid
#@markdown surface glia used as a positive control. Just run all three cells.
#@markdown To choose different runs, clear RUNS, run Cells 1 and 2 to see the full list, then type the names here.
RUNS = "DFB__303662__10x_scATAC_DG2_larval_brain, DFB__96a6e8__10x_scATAC_DG2_larval_brain"  #@param {type:"string"}
DATA_DIR = os.environ.get("ROUTE_B_DATA", DATA_DIR)
if DATA_DIR.startswith("/content/drive"):
    from google.colab import drive
    drive.mount("/content/drive")
os.makedirs(DATA_DIR, exist_ok=True)


class StopHere(Exception):
    """Stops the notebook with a plain-language message."""


WANTED = ("fragments.tsv.gz", "singlecell.csv")
print("Listing:", BASE_URL, "| saving to:", DATA_DIR)
