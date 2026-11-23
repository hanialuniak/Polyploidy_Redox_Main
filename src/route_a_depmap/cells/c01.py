#@title Cell 1 - Setup, settings and connect Google Drive
import io, os, re, warnings
import numpy as np
import pandas as pd
import requests
import matplotlib.pyplot as plt
warnings.filterwarnings("ignore")
try:
    from IPython.display import display
except ImportError:
    display = print

# Your data folder in Google Drive (where you put the 7 downloaded files)
DATA_DIR = "/content/drive/MyDrive/Colab Notebooks/route_a_data"  #@param {type:"string"}
DATA_DIR = os.environ.get("ROUTE_A_DATA", DATA_DIR)
if DATA_DIR.startswith("/content/drive"):
    from google.colab import drive
    drive.mount("/content/drive")


class StopHere(Exception):
    """Stops the notebook with a plain-language message."""


if not os.path.isdir(DATA_DIR):
    raise StopHere(f"Folder not found: {DATA_DIR}. Check the folder name in Google Drive and fix DATA_DIR above.")

# Settings (Colab shows these as boxes on the right)
LINEAGE = "All lineages"  #@param ["All lineages", "CNS/Brain only"]
USE_PROTEINS = True  #@param {type:"boolean"}
EXCLUDE_NRF2_MUTANTS = True  #@param {type:"boolean"}
MY_GENES = ""  #@param {type:"string"}
N_RANDOM_SETS = 500  #@param {type:"integer"}
N_PERMUTATIONS = 2000  #@param {type:"integer"}
MIN_LINES_PER_TISSUE = 10  #@param {type:"integer"}
DEPMAP_RELEASE = "DepMap Public 24Q4"  #@param {type:"string"}
SEED = 1  #@param {type:"integer"}
GENES_OF_INTEREST = "GMNN, CDT1"  #@param {type:"string"}
USE_CRISPR = True  #@param {type:"boolean"}
#@markdown **TP53**: p53 status affects both ploidy tolerance and some redox genes.
TP53_MODE = "All lines"  #@param ["All lines", "TP53-normal lines only", "Adjust for TP53 status"]
#@markdown **Proliferation adjustment**: untick to see how much it changes the results.
ADJUST_PROLIFERATION = True  #@param {type:"boolean"}

# Reporting thresholds (conventions for labelling findings, not hypotheses)
FDR_MAX = 0.05            # false discovery rate
RANDOM_P_MAX = 0.05       # must beat expression-matched random gene sets
CONTROL_P_MAX = 0.05      # must beat every non-redox stress programme

rng = np.random.default_rng(SEED)
GOI = [g.strip().upper() for g in re.split(r"[,\s;]+", GENES_OF_INTEREST) if g.strip()]
print("Genes of interest:", GOI, "| TP53:", TP53_MODE, "| proliferation adjustment:", ADJUST_PROLIFERATION)
print("Data folder:", DATA_DIR)
