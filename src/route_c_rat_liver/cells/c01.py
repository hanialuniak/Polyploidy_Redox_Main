#@title Cell 1 - Settings and connect Google Drive
import gzip, io, os, re, subprocess, warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
warnings.filterwarnings("ignore")
try:
    from IPython.display import display
except ImportError:
    display = print

DATA_DIR = "/content/drive/MyDrive/Colab Notebooks/route_c_rat"  #@param {type:"string"}
GSE = "GSE132409"  #@param {type:"string"}
N_RANDOM_SETS = 2000  #@param {type:"integer"}
SEED = 1  #@param {type:"integer"}
DATA_DIR = os.environ.get("ROUTE_CR_DATA", DATA_DIR)
if DATA_DIR.startswith("/content/drive"):
    from google.colab import drive
    drive.mount("/content/drive")
os.makedirs(DATA_DIR, exist_ok=True)
rng = np.random.default_rng(SEED)


class StopHere(Exception):
    """Stops the notebook with a plain-language message."""


print("Dataset:", GSE, "| folder:", DATA_DIR)
