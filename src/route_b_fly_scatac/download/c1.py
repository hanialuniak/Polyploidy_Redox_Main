#@title Cell 1 - Settings and connect Google Drive
import os, re, subprocess, urllib.parse
import pandas as pd
try:
    from IPython.display import display
except ImportError:
    display = print

BASE_URL = "https://flybrain.aertslab.org/downloads/cellRanger/"  #@param {type:"string"}
DATA_DIR = "/content/drive/MyDrive/Colab Notebooks/route_b_fly"  #@param {type:"string"}
#@markdown Leave **RUNS** empty the first time: Cell 2 lists what is available. Then type the run names you want
#@markdown (comma-separated) and run Cell 3.
RUNS = ""  #@param {type:"string"}
DATA_DIR = os.environ.get("ROUTE_B_DATA", DATA_DIR)
if DATA_DIR.startswith("/content/drive"):
    from google.colab import drive
    drive.mount("/content/drive")
os.makedirs(DATA_DIR, exist_ok=True)


class StopHere(Exception):
    """Stops the notebook with a plain-language message."""


WANTED = ("fragments.tsv.gz", "singlecell.csv")
print("Listing:", BASE_URL, "| saving to:", DATA_DIR)
