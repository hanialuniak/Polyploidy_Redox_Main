#@title Cell 13 - How many of those "shared" genes would appear if the germline split were meaningless?
# Cell 12 removed the ambient-mRNA artefacts and still left ~700 genes called shared with the germline. That
# number has never been calibrated, and it has to be, because the two splits are not independent in the way
# the argument assumes: both separate "more mitotic" from "less mitotic, more developed" cells, the germline
# split was made per cell on an M-phase score inside a single cluster, and any gene that tracks cell size,
# library depth or maturity will move on both sides for reasons that have nothing to do with endocycling.
#
# So: shuffle the germline labels, keeping the group sizes and the library structure, recompute the germline
# fold changes, and see how many genes get called shared anyway. That gives the number expected by chance, and
# per gene it gives an empirical p-value for the germline change - which is the thing Cell 11 could not offer.
N_PERM = 200        #@param {type:"integer"}

_gidx = np.flatnonzero(((cells.gl_role.isin([GMIT, GENDO])) & low_g).to_numpy())
_glib = cells.library.to_numpy()[_gidx]
_grole = cells.gl_role.to_numpy()[_gidx]
Xg = X[_gidx][:, [gpos[g] for g in Lg.index]]
keep_lib = sorted({l for l in set(_glib)
                   if ((_glib == l) & (_grole == GMIT)).sum() >= 30
                   and ((_glib == l) & (_grole == GENDO)).sum() >= 30})
print(f"{len(_gidx):,} germline cells across {len(keep_lib)} librar{'y' if len(keep_lib) == 1 else 'ies'} "
      f"with both groups. Shuffling labels within each library {N_PERM} times.")


def gl_fc(role_vec):
    """The same pseudobulk and median-of-ratios normalisation as Cell 10, on a given set of labels."""
    cols, lab2 = [], []
    for l in keep_lib:
        for r in (GMIT, GENDO):
            m = (_glib == l) & (role_vec == r)
            cols.append(np.asarray(Xg[np.flatnonzero(m)].sum(axis=0)).ravel()); lab2.append(r)
    P = np.array(cols).T.astype(float)
    with np.errstate(divide="ignore", invalid="ignore"):
        lp = np.log(np.where(P > 0, P, np.nan))
        ref = np.nanmean(lp, axis=1)
        good = np.isfinite(ref) & np.isfinite(lp).all(axis=1)
        size = np.exp(np.median(lp[good] - ref[good, None], axis=0))
    L = np.log2(P / size / P.sum(axis=0).mean() * 1e6 + 1)
    lab2 = np.array(lab2)
    return L[:, lab2 == GENDO].mean(axis=1) - L[:, lab2 == GMIT].mean(axis=1)


obs = pd.Series(gl_fc(_grole), index=Lg.index)
rng = np.random.default_rng(7)
null = np.empty((N_PERM, len(Lg.index)), dtype=np.float32)
for i in range(N_PERM):
    perm = _grole.copy()
    for l in keep_lib:
        m = _glib == l
        perm[m] = rng.permutation(perm[m])
    null[i] = gl_fc(perm)
    if (i + 1) % 50 == 0:
        print(f"   {i + 1}/{N_PERM}")

# Per-gene p for the germline change. A gene's own 200 shuffles can only reach p = 1/201, which after
# correcting across thousands of genes can never clear any sensible FDR - so the null is pooled across genes
# instead, which is the standard permutation approach and gives the resolution BH needs. Pooling is only
# legitimate between genes whose noise is comparable, and noise scales with expression, so the pool is built
# within expression deciles - the same matching logic used for the random gene sets everywhere else here.
lvl = Lg.mean(axis=1)
bin_id = pd.qcut(lvl.rank(method="first"), 10, labels=False)
aobs = np.abs(obs.to_numpy())
p_emp = np.ones(len(obs))
for b in range(10):
    m = np.flatnonzero((bin_id == b).to_numpy())
    pool = np.sort(np.abs(null[:, m]).ravel())
    p_emp[m] = (len(pool) - np.searchsorted(pool, aobs[m], side="left") + 1) / (len(pool) + 1)
p_emp = pd.Series(p_emp, index=Lg.index)
o = p_emp.sort_values()
q_emp = pd.Series(np.minimum.accumulate((o.to_numpy() * len(o) / np.arange(1, len(o) + 1))[::-1])[::-1],
                  index=o.index).reindex(Lg.index).clip(upper=1.0)
print(f"Per-gene null pooled within expression deciles: {N_PERM * (bin_id == 0).sum():,} null values per bin.")

de3 = de2.join(pd.DataFrame({"germline_p": p_emp, "germline_q": q_emp}))
de3["shared_calibrated"] = (de3.lineage.astype(str).str.startswith("BOTH")) & \
                           (de3.bleed_through_risk == "ok") & (de3.germline_q < 0.10)

# How many shared calls does a meaningless germline split produce? Rebuild the call on each permutation.
f_moves = (de3.log2FC.abs() >= 0.3) & ((de3.follicle_pct > 65) | (de3.follicle_pct < 35))
f_moves = f_moves.reindex(Lg.index).fillna(False).to_numpy()
f_sign = np.sign(de3.log2FC.reindex(Lg.index).fillna(0).to_numpy())
ok_bleed = (de3.bleed_through_risk.reindex(Lg.index) == "ok").to_numpy()


def n_shared(v):
    gp = pd.Series(v).rank(pct=True).to_numpy() * 100
    g_moves = (np.abs(v) >= 0.3) & ((gp > 65) | (gp < 35))
    return int((f_moves & g_moves & ok_bleed & (np.sign(v) == f_sign)).sum())


obs_n = n_shared(obs.to_numpy())
null_n = np.array([n_shared(null[i]) for i in range(N_PERM)])
print(f"\nShared calls observed: {obs_n:,}")
print(f"Shared calls from shuffled germline labels: median {int(np.median(null_n)):,}, "
      f"95th percentile {int(np.percentile(null_n, 95)):,}, range {null_n.min():,}-{null_n.max():,}")
print(f"Empirical p for the count: {((null_n >= obs_n).sum() + 1) / (N_PERM + 1):.3f}")
if np.median(null_n) > 0.5 * obs_n:
    print("\n*** A meaningless split produces most of the shared list. The COUNT is not evidence of shared "
          "regulation, and no gene should be believed on the strength of appearing in it. Use the per-gene "
          "q-values below instead, which is what they are for. ***")
else:
    print("\nThe observed count is well above what shuffling produces, so there is shared signal in the list "
          "as a whole. Individual genes still need their own q-value.")

print(f"\nShared, ambient-checked and calibrated at q < 0.10: {int(de3.shared_calibrated.sum()):,} genes")
best = de3[de3.shared_calibrated]
best = best.reindex(best.log2FC.abs().sort_values(ascending=False).index)
show3 = [c for c in ["log2FC", "germline_log2FC", "germline_q", "gl_detect_max", "expr_rank_pct", "in_set"]
         if c in best.columns]
display(best.head(50)[show3].round(3))

WATCH2 = [g for g in "Rrp1 mre11 geminin Jafrac1 Jafrac2 Ero1L Gclc Trx-2 Trxr-1 Cat Hsc70-3 Sec61beta "
                     "HmgD Df31 D1 Lam ball Bub3 Det".split() if g in de3.index]
print("\nThe genes being followed, now with an empirical p-value on the germline side:")
display(de3.loc[WATCH2, [c for c in ["log2FC", "germline_log2FC", "germline_p", "germline_q", "lineage"]
                         if c in de3.columns]].round(3))
print("""
What this cell can and cannot do. It calibrates the GERMLINE side, which is the weak one: 300-odd cells split
by a score rather than by clustering. It does not calibrate the follicle side, which does not need it - those
effects are large and agree across all three libraries. A gene that survives here has a follicle change, a
germline change in the same direction that shuffling does not reproduce, and enough germline expression that
ambient mRNA is not the explanation. That is the strongest statement these data support about any single gene.
""")
de3.sort_values("log2FC").to_csv(os.path.join(GDIR, "route_e_unbiased_de_calibrated.csv"))
print("Saved to", os.path.join(GDIR, "route_e_unbiased_de_calibrated.csv"))
