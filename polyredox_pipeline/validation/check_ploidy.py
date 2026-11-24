import sys, numpy as np, pandas as pd
from sklearn.metrics import roc_auc_score, average_precision_score
from polyredox import ploidy
S = sys.argv[1]
kappa = float(sys.argv[2]) if len(sys.argv) > 2 else 8.0
f = pd.read_csv(S + "/fedecs_simulated.tsv.gz", sep="\t", index_col=0)
c = pd.read_csv(S + "/cells.tsv", sep="\t"); t = pd.read_csv(S + "/truth.tsv", sep="\t")
post, par = ploidy.call_ploidy(f, c.cell_class.values, kappa=kappa)
d = post.join(t.set_index("cell_id")).assign(cls=c.cell_class.values)
print("not called", d.P_polyploid.isna().sum())
d = d.dropna(subset=["P_polyploid"]); d["poly"] = d.true_copies > 2
for cl, x in d.groupby("cls"):
    print(cl, "n", len(x), "true rate", round(x.poly.mean(), 3), "(doublets", round(x.is_doublet.mean(), 3), ") mean P", round(x.P_polyploid.mean(), 3),
          "AUC", round(roc_auc_score(x.poly, x.P_polyploid), 3), "AP", round(average_precision_score(x.poly, x.P_polyploid), 3))
d["depth_bin"] = pd.qcut(d.n_fragments, 4)
print(d.groupby("depth_bin", observed=True).apply(lambda x: pd.Series({"AUC": roc_auc_score(x.poly, x.P_polyploid), "true": x.poly.mean(), "meanP": x.P_polyploid.mean()})).round(3))
print(pd.crosstab(d.ploidy_call, d.true_copies))
for k, p in par.items(): print(k, {kk: np.round(v, 4) for kk, v in p.items()})
