#@title Cell 1 - Setup, settings and connect Google Drive
import glob, gzip, io, os, re, subprocess, tarfile, warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
warnings.filterwarnings("ignore")
try:
    from IPython.display import display
except ImportError:
    display = print

DATA_DIR = "/content/drive/MyDrive/Colab Notebooks/route_b_fly_bulk"  #@param {type:"string"}
GTF_PATH = "/content/drive/MyDrive/Colab Notebooks/route_b_fly/Drosophila_melanogaster.gtf.gz"  #@param {type:"string"}
DATA_DIR = os.environ.get("ROUTE_BB_DATA", DATA_DIR)
GTF_PATH = os.environ.get("ROUTE_BB_GTF", GTF_PATH)
if DATA_DIR.startswith("/content/drive"):
    from google.colab import drive
    drive.mount("/content/drive")
os.makedirs(DATA_DIR, exist_ok=True)


class StopHere(Exception):
    """Stops the notebook with a plain-language message."""


GSE = "GSE153165"
N_RANDOM_SETS = 500  #@param {type:"integer"}
SEED = 1  #@param {type:"integer"}
FDR_MAX, RANDOM_P_MAX, CONTROL_P_MAX = 0.05, 0.05, 0.05
rng = np.random.default_rng(SEED)
print("Data folder:", DATA_DIR)
