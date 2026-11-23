#@title Cell 4 - Normalise, cluster, and name the clusters from marker genes
from sklearn.decomposition import TruncatedSVD
from sklearn.cluster import KMeans

gpos = {g: i for i, g in enumerate(genes)}
cpm = X.multiply(1e4 / np.maximum(cells.n_counts.to_numpy()[:, None], 1)).tocsr()
Ln = cpm.copy(); Ln.data = np.log1p(Ln.data)                       # log normalised expression per cell
v = np.asarray(Ln.power(2).mean(axis=0)).ravel() - np.asarray(Ln.mean(axis=0)).ravel() ** 2
hv = np.argsort(v)[-2000:]
Z = Ln[:, hv].toarray()
Z = (Z - Z.mean(0)) / np.where(Z.std(0) > 0, Z.std(0), 1)
svd = TruncatedSVD(n_components=30, random_state=SEED).fit_transform(np.clip(Z, -10, 10))
svd = svd / np.maximum(np.linalg.norm(svd, axis=1, keepdims=True), 1e-9)
cells["cluster"] = KMeans(n_clusters=N_CLUSTERS, n_init=5, random_state=SEED).fit_predict(svd).astype(str)


def score(gs):
    """Mean normalised expression of a gene set, per cell."""
    idx = [gpos[g] for g in gs if g in gpos]
    return np.asarray(Ln[:, idx].mean(axis=1)).ravel() if idx else np.zeros(len(cells)), len(idx)


MARKERS = {
    "follicle cell (somatic epithelium)": "eya br Fas3 CadN Fas2 Sox14 mirr pnt",
    "germline": "vas nos bam orb osk Rbp9 corolla",
    "escort / inner germarial sheath": "ptc hh tj Wnt4 bnl",
    "muscle sheath": "Mhc Mlc1 up zip Act57B",
    "stalk and polar cells": "upd1 LamC cas Fas3",
    "mitotic (cycling)": "stg CycB polo aurB Det pav PCNA CycA Incenp",
    "early follicle (stages 1-6)": "CadN Fas2 ct SPARC",
    "mid follicle (stages 7-10, endocycling)": "br mirr mid Yp1 Yp2 Yp3",
    "late follicle (stages 9+, endocycling)": "Sox14 Vha16-1 Pvf1 dpp cv-2",
    "vitelline membrane / chorion (amplifying)": "Vm26Aa Vm26Ab Vm32E Vm34Ca Cp36 Cp38 Cp15 Cp16 dec-1 psd",
}
tab = {}
for name, gs in MARKERS.items():
    s, n = score(gs.split())
    tab[f"{name} ({n})"] = pd.Series(s).groupby(cells.cluster.to_numpy()).mean()
tab = pd.DataFrame(tab)
zt = (tab - tab.mean()) / tab.std().replace(0, 1)
zt["cells"] = cells.cluster.value_counts()
display(zt.round(2).sort_values(zt.columns[0], ascending=False))
print("Columns are z-scores across clusters: a high value means that cluster expresses those markers most.")
