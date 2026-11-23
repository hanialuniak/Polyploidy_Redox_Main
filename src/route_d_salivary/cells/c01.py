#@title Cell 1 - Settings and connect Google Drive
import gzip, io, os, re, subprocess, tarfile, warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
warnings.filterwarnings("ignore")
try:
    from IPython.display import display
except ImportError:
    display = print

#@markdown **GSE121956**: salivary glands (endocycling) vs brain and imaginal discs (dividing), feeding 3rd-instar larvae.
#@markdown **GSE121955**: S2 cells pushed into endocycles by Cyclin A or Myb knockdown vs GFP-knockdown controls.
GSE = "GSE121956"  #@param ["GSE121956", "GSE121955"]
DATA_DIR = "/content/drive/MyDrive/Colab Notebooks/route_d_salivary"  #@param {type:"string"}
GTF_PATH = "/content/drive/MyDrive/Colab Notebooks/route_b_fly/Drosophila_melanogaster.gtf.gz"  #@param {type:"string"}
N_RANDOM_SETS = 2000  #@param {type:"integer"}
SEED = 1  #@param {type:"integer"}
DATA_DIR = os.environ.get("ROUTE_D_DATA", DATA_DIR)
GTF_PATH = os.environ.get("ROUTE_D_GTF", GTF_PATH)
if DATA_DIR.startswith("/content/drive"):
    from google.colab import drive
    drive.mount("/content/drive")
GDIR = os.path.join(DATA_DIR, GSE)
os.makedirs(GDIR, exist_ok=True)
rng = np.random.default_rng(SEED)
pd.set_option("display.width", 220); pd.set_option("display.max_columns", 14); pd.set_option("display.max_rows", 200)


class StopHere(Exception):
    """Stops the notebook with a plain-language message."""


# Which samples are endocycling ("test") and which are dividing ("reference"), recognised from the GEO descriptions
GROUP_RULES = {
    "GSE121956": [("salivary gland", r"(?i:salivary)|(?<![A-Za-z])SG(?![a-z])", "test"),
                  ("brain + imaginal disc", r"(?i:brain|disc|disk|wing|eye)|(?<![A-Za-z])BD(?![a-z])", "reference")],
    "GSE121955": [("CycA knockdown", r"(?i:cyc\s*-?\s*a|cyclin\s*a)", "test"), ("Myb knockdown", r"(?i:myb)", "test"),
                  ("GFP control", r"(?i:gfp|control)", "reference")],
}
print("Dataset:", GSE, "| folder:", GDIR)
