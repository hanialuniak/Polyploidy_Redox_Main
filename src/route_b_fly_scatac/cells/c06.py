#@title Cell 6 - Ploidy of every cell, from overlapping fragments
GRID = np.exp(np.linspace(np.log(5e-4), np.log(0.7), 90)); LOG_GRID = np.log(GRID)


def _log_categories(p, kappa, artefact, k):
    t = p - 1; a, b = GRID * kappa, (1 - GRID) * kappa
    out = np.zeros((len(GRID), k))
    for s in range(t + 1):
        out[:, min(s, k - 1)] += np.exp(gammaln(t + 1) - gammaln(s + 1) - gammaln(t - s + 1)
                                        + betaln(s + a, t - s + b) - betaln(a, b))
    out = np.maximum(out + artefact[None, :], 1e-300)
    return np.log(out / out.sum(axis=1, keepdims=True))


def _basis(x, knots):
    return np.column_stack([np.ones_like(x), x] + [np.maximum(x - kn, 0) for kn in knots])


def _e_step(C, n, weights, th, kappa):
    lm = np.zeros((len(C), len(PLOIDIES))); ms = np.zeros_like(lm); vs = np.zeros_like(lm)
    for j, p in enumerate(PLOIDIES):
        lp = _log_categories(p, kappa, th["art"], C.shape[1])
        mean = _basis(np.log(2 * n / p), th["knots"]) @ th["coef"]
        joint = C @ lp.T - 0.5 * ((LOG_GRID[None, :] - mean[:, None]) / th["tau"]) ** 2
        norm = logsumexp(joint, axis=1, keepdims=True); post = np.exp(joint - norm)
        lm[:, j] = norm[:, 0]; ms[:, j] = post @ LOG_GRID; vs[:, j] = post @ LOG_GRID**2 - ms[:, j] ** 2
    lj = lm + np.log(weights)
    return np.exp(lj - logsumexp(lj, axis=1, keepdims=True)), ms, vs


def fit_ploidy(C, kappa=8.0, n_iter=40):
    C = np.asarray(C, float); n = C.sum(axis=1); k = C.shape[1]
    rate = np.maximum((C[:, 1] + 0.5) / (C[:, 0] + 1), 1e-4)
    knots = np.quantile(np.log(n), [0.2, 0.4, 0.6, 0.8])
    th = {"coef": np.linalg.lstsq(_basis(np.log(n), knots), np.log(rate), rcond=None)[0], "knots": knots,
          "tau": 0.5, "art": np.r_[0, 0, np.full(k - 2, 1e-4)]}
    w = np.r_[0.9, np.full(len(PLOIDIES) - 1, 0.1 / (len(PLOIDIES) - 1))]
    for _ in range(n_iter):
        resp, ms, vs = _e_step(C, n, w, th, kappa)
        new_w = np.maximum(resp.mean(0), 1e-5); new_w /= new_w.sum()
        X = np.concatenate([np.log(2 * n / p) for p in PLOIDIES]); Y = ms.T.ravel()
        W = resp.T.ravel(); V = vs.T.ravel(); B = _basis(X, knots); r = np.sqrt(W)
        coef = np.linalg.lstsq(B * r[:, None], Y * r, rcond=None)[0]
        tau = max(float(np.sqrt(np.average((Y - B @ coef) ** 2 + V, weights=W))), 0.05)
        d = resp[:, 0]; art = np.r_[0, 0, np.maximum((d[:, None] * C[:, 2:]).sum(0) / max((d * n).sum(), 1), 1e-7)]
        done = np.max(np.abs(new_w - w)) < 1e-4
        w, th = new_w, {"coef": coef, "knots": knots, "tau": tau, "art": art}
        if done:
            break
    resp, ms, _ = _e_step(C, n, w, th, kappa)
    return resp, (resp * ms).sum(1)


for p in PLOIDIES:
    cells[f"P_{p}C"] = np.nan
for cl, part in cells.groupby("cluster"):
    resp, log_s = fit_ploidy(overlaps.loc[part.cell_id].to_numpy())
    for j, p in enumerate(PLOIDIES):
        cells.loc[part.index, f"P_{p}C"] = resp[:, j]
    cells.loc[part.index, "per_copy_log_s"] = log_s
cells["P_polyploid"] = 1 - cells["P_2C"]
cells["expected_copies"] = sum(cells[f"P_{p}C"] * p for p in PLOIDIES)
cells["fragments_per_copy"] = cells.n_fragments * 2 / cells.expected_copies
best = cells[[f"P_{p}C" for p in PLOIDIES]].to_numpy()
cells["ploidy_call"] = np.where(best.max(1) >= 0.8, np.array([f"{p}C" for p in PLOIDIES])[best.argmax(1)], "uncertain")
display(pd.crosstab(cells.cell_type, cells.ploidy_call))
