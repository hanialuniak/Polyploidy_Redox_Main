#@title Cell 1 - Settings and connect Google Drive
import gzip, io, os, re, subprocess, tarfile, warnings
import numpy as np
import pandas as pd
import scipy.sparse as sp
import matplotlib.pyplot as plt
warnings.filterwarnings("ignore")
try:
    from IPython.display import display
except ImportError:
    display = print

#@markdown **GSE146040**: adult ovary single cells, Jevitt et al. 2020 (Deng lab), ~35,000 cells.
#@markdown **GSE136162**: adult ovary single cells, Rust et al. 2020 (Nystul lab).
#@markdown Both hold mitotic follicle cells (stages 1-6) and endocycling follicle cells (stages 7-10) from the same tissue.
GSE = "GSE146040"  #@param ["GSE146040", "GSE136162"]
DATA_DIR = "/content/drive/MyDrive/Colab Notebooks/route_e_ovary"  #@param {type:"string"}
GTF_PATH = "/content/drive/MyDrive/Colab Notebooks/route_b_fly/Drosophila_melanogaster.gtf.gz"  #@param {type:"string"}
N_CLUSTERS = 18  #@param {type:"integer"}
MIN_GENES = 400  #@param {type:"integer"}
MIN_COUNTS = 1000  #@param {type:"integer"}
N_RANDOM_SETS = 2000  #@param {type:"integer"}
SEED = 1  #@param {type:"integer"}
DATA_DIR = os.environ.get("ROUTE_E_DATA", DATA_DIR)
GTF_PATH = os.environ.get("ROUTE_E_GTF", GTF_PATH)
if DATA_DIR.startswith("/content/drive"):
    from google.colab import drive
    drive.mount("/content/drive")
GDIR = os.path.join(DATA_DIR, GSE)
os.makedirs(GDIR, exist_ok=True)
rng = np.random.default_rng(SEED)
pd.set_option("display.width", 230); pd.set_option("display.max_columns", 16); pd.set_option("display.max_rows", 220)


class StopHere(Exception):
    """Stops the notebook with a plain-language message."""


print("Dataset:", GSE, "| folder:", GDIR)
