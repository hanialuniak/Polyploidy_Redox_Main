#@title Cell 1 - Setup, settings and connect Google Drive
import gzip, os, re, warnings
import numpy as np
import pandas as pd
import scipy.sparse as sp
import matplotlib.pyplot as plt
from scipy.special import betaln, gammaln, logsumexp
warnings.filterwarnings("ignore")
try:
    from IPython.display import display
except ImportError:
    display = print

DATA_DIR = "/content/drive/MyDrive/Colab Notebooks/route_b_fly"  #@param {type:"string"}
DATA_DIR = os.environ.get("ROUTE_B_DATA", DATA_DIR)
if DATA_DIR.startswith("/content/drive"):
    from google.colab import drive
    drive.mount("/content/drive")


class StopHere(Exception):
    """Stops the notebook with a plain-language message."""


if not os.path.isdir(DATA_DIR):
    raise StopHere(f"Folder not found: {DATA_DIR}. Check the folder name in Google Drive and fix DATA_DIR above.")

#@markdown **Runs**: one sub-folder per sequencing run, each holding fragments.tsv.gz and singlecell.csv.
#@markdown Prefer runs that include the optic lobes and flies at least 1-3 weeks old (see the guide).
RUNS = "FDM__f2501a__10x_scATAC_mix_DGRP_adult_brain, FDM__1d60cb__10x_scATAC_mix_DGRP_adult_central_brain, FDM__9c338b__10x_scATAC_mix_DGRP_adult_central_brain"  #@param {type:"string"}
GTF_FILE = "Drosophila_melanogaster.gtf.gz"  #@param {type:"string"}
N_CLUSTERS = 20  #@param {type:"integer"}
MIN_FRAGMENTS = 2000  #@param {type:"integer"}
MAX_FRAGMENTS = 60000  #@param {type:"integer"}
GENES_OF_INTEREST = "geminin, dup, Cdc6, Orc1"  #@param {type:"string"}
N_RANDOM_SETS = 500  #@param {type:"integer"}
N_PERMUTATIONS = 2000  #@param {type:"integer"}
SEED = 1  #@param {type:"integer"}

FDR_MAX, RANDOM_P_MAX, CONTROL_P_MAX = 0.05, 0.05, 0.05
PLOIDIES = (2, 4, 8, 16)
AUTOSOMES = {"2L", "2R", "3L", "3R"}
BIN = 500            # base pairs per bin for gene accessibility
PROMOTER = 1000      # base pairs upstream of each gene included
rng = np.random.default_rng(SEED)
RUN_LIST = [r.strip() for r in RUNS.split(",") if r.strip()]
GOI = [g.strip() for g in re.split(r"[,\s;]+", GENES_OF_INTEREST) if g.strip()]
print("Runs:", RUN_LIST, "| genes of interest:", GOI, "| data folder:", DATA_DIR)
