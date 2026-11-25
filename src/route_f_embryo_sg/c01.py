#@title Cell 1 - Settings and connect Google Drive
import gzip, os, re, subprocess, warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
warnings.filterwarnings("ignore")
try:
    from IPython.display import display
except ImportError:
    display = print
pd.set_option("display.width", 230); pd.set_option("display.max_columns", None); pd.set_option("display.max_rows", 220)

#@markdown **GSE234602** (Peng, Jackson et al. 2024, *Development*): single-cell RNA-seq of the *Drosophila*
#@markdown embryo across two windows, **stages 10-12** and **stages 13-16**, with the salivary gland annotated.
#@markdown Both windows come from the same study, so the early-to-late comparison carries no batch confound.
#@markdown The HyPer7 nuclear signal in salivary glands is reported at stage 14, in the late window.
GSE = "GSE234602"  #@param {type:"string"}
DATA_DIR = "/content/drive/MyDrive/Colab Notebooks/route_f_embryo_sg"  #@param {type:"string"}
#@markdown Which annotated cell type counts as salivary gland (a regular expression, case-insensitive).
TISSUE_PATTERN = "saliv"  #@param {type:"string"}
#@markdown The gene annotation converts FlyBase IDs to symbols. It is the same file Route B and Route D use.
GTF_PATH = "/content/drive/MyDrive/Colab Notebooks/route_b_fly/Drosophila_melanogaster.gtf.gz"  #@param {type:"string"}
GTF_PATH = os.environ.get("ROUTE_F_GTF", GTF_PATH)
#@markdown Minimum raw counts a gene needs in EVERY pseudobulk group. The salivary gland groups hold only a
#@markdown few hundred cells, so a high value here silently removes well-measured genes.
MIN_COUNTS_PER_GROUP = 5  #@param {type:"integer"}
#@markdown Minimum counts summed across all four groups, which keeps the ratios stable.
MIN_TOTAL_COUNTS = 50  #@param {type:"integer"}
N_RANDOM_SETS = 2000  #@param {type:"integer"}
SEED = 1  #@param {type:"integer"}

DATA_DIR = os.environ.get("ROUTE_F_DATA", DATA_DIR)
if DATA_DIR.startswith("/content/drive"):
    from google.colab import drive
    drive.mount("/content/drive")
GDIR = os.path.join(DATA_DIR, GSE)
os.makedirs(GDIR, exist_ok=True)
rng = np.random.default_rng(SEED)


class StopHere(Exception):
    """Stops the notebook with a plain-language message."""


print("Dataset:", GSE, "| folder:", GDIR)
print("Salivary gland cells will be those whose annotated cell type matches:", repr(TISSUE_PATTERN))
