"""Probabilistic ploidy calls from fragment overlap counts.

Given a fragment end on one genome copy, each of the other p-1 copies covers that
position with probability s, so x-1 ~ BetaBinomial(p-1, s, kappa). Three things are
learned for each cell class by EM:
  * the expected s as a function of fragments per genome copy,
  * the per-fragment rate of depth >= 3 mapping artefacts in diploid nuclei,
  * the proportions of each ploidy.

This differs from scPloidy's unconstrained categorical mixture. Polyploid components
are tied to the binomial form so that rare polyploid neurons cannot be displaced by a
component that merely splits diploid nuclei by depth.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.special import betaln, gammaln, logsumexp

GRID = np.exp(np.linspace(np.log(5e-4), np.log(0.7), 90))
LOG_GRID = np.log(GRID)


def _log_categories(ploidy, kappa, artefact, k):
    """Return log P(x category | ploidy, s) with shape len(GRID) x k."""
    trials = ploidy - 1
    a, b = GRID * kappa, (1 - GRID) * kappa
    probabilities = np.zeros((len(GRID), k))
    for successes in range(trials + 1):
        log_pmf = (
            gammaln(trials + 1) - gammaln(successes + 1) - gammaln(trials - successes + 1)
            + betaln(successes + a, trials - successes + b) - betaln(a, b)
        )
        probabilities[:, min(successes, k - 1)] += np.exp(log_pmf)
    probabilities = probabilities + artefact[None, :]
    probabilities = np.maximum(probabilities, 1e-300)
    return np.log(probabilities / probabilities.sum(axis=1, keepdims=True))


def _basis(log_depth, knots):
    """Piecewise-linear basis in log fragments per diploid genome."""
    columns = [np.ones_like(log_depth), log_depth]
    columns += [np.maximum(log_depth - k, 0) for k in knots]
    return np.column_stack(columns)


def _log_prior(n_fragments, ploidy, coef, knots, tau):
    mean = _basis(np.log(2 * n_fragments / ploidy), knots) @ coef
    return -0.5 * ((LOG_GRID[None, :] - mean[:, None]) / tau) ** 2 - np.log(tau)


def _e_step(counts, n_fragments, ploidies, weights, theta, kappa, batch=20000):
    k = counts.shape[1]
    n = counts.shape[0]
    log_marginal = np.zeros((n, len(ploidies)))
    mean_log_s = np.zeros_like(log_marginal)
    var_log_s = np.zeros_like(log_marginal)
    for j, ploidy in enumerate(ploidies):
        log_p = _log_categories(ploidy, kappa, theta["artefact"], k)
        for left in range(0, n, batch):
            right = min(left + batch, n)
            joint = counts[left:right] @ log_p.T + _log_prior(
                n_fragments[left:right], ploidy, theta["coef"], theta["knots"], theta["tau"]
            )
            norm = logsumexp(joint, axis=1, keepdims=True)
            log_marginal[left:right, j] = norm[:, 0]
            posterior = np.exp(joint - norm)
            mean_log_s[left:right, j] = posterior @ LOG_GRID
            var_log_s[left:right, j] = posterior @ LOG_GRID**2 - mean_log_s[left:right, j] ** 2
    log_joint = log_marginal + np.log(weights)
    resp = np.exp(log_joint - logsumexp(log_joint, axis=1, keepdims=True))
    return resp, mean_log_s, var_log_s


def fit_group(counts, n_fragments, ploidies=(2, 4, 8), kappa=8.0, n_iter=40, tol=1e-4, seed=0, max_fit=60000):
    """Fit one cell class; return (posterior, theta, posterior mean log s per nucleus).

    s is the per-copy overlap probability, a ploidy-independent measure of library quality.
    """
    counts = np.asarray(counts, dtype=float)
    n_fragments = np.asarray(n_fragments, dtype=float)
    k = counts.shape[1]
    ploidies = list(ploidies)
    rng = np.random.default_rng(seed)
    subset = (
        rng.choice(len(counts), max_fit, replace=False) if len(counts) > max_fit else np.arange(len(counts))
    )
    fit_counts, fit_n = counts[subset], n_fragments[subset]

    overlap_rate = np.maximum((fit_counts[:, 1] + 0.5) / (fit_counts[:, 0] + 1.0), 1e-4)
    knots = np.quantile(np.log(fit_n), [0.2, 0.4, 0.6, 0.8])
    coef = np.linalg.lstsq(_basis(np.log(fit_n), knots), np.log(overlap_rate), rcond=None)[0]
    weights = np.array([0.97] + [0.03 / (len(ploidies) - 1)] * (len(ploidies) - 1))
    artefact = np.r_[0, 0, np.full(k - 2, 1e-4)]
    theta = {"coef": coef, "knots": knots, "tau": 0.5, "artefact": artefact}

    for _ in range(n_iter):
        resp, mean_log_s, var_log_s = _e_step(fit_counts, fit_n, ploidies, weights, theta, kappa)
        new_weights = np.maximum(resp.mean(axis=0), 1e-5)
        new_weights /= new_weights.sum()
        design, target, weight, variance = [], [], [], []
        for j, ploidy in enumerate(ploidies):
            design.append(np.log(2 * fit_n / ploidy))
            target.append(mean_log_s[:, j])
            weight.append(resp[:, j])
            variance.append(var_log_s[:, j])
        design, target = np.concatenate(design), np.concatenate(target)
        weight, variance = np.concatenate(weight), np.concatenate(variance)
        basis = _basis(design, knots)
        root = np.sqrt(weight)
        coef = np.linalg.lstsq(basis * root[:, None], target * root, rcond=None)[0]
        residual = target - basis @ coef
        tau = max(float(np.sqrt(np.average(residual**2 + variance, weights=weight))), 0.05)
        diploid = resp[:, 0]
        denominator = max(float(np.sum(diploid * fit_n)), 1.0)
        new_artefact = np.r_[0, 0, np.maximum(np.sum(diploid[:, None] * fit_counts[:, 2:], axis=0) / denominator, 1e-7)]
        change = float(np.max(np.abs(new_weights - weights)))
        weights = new_weights
        theta = {"coef": coef, "knots": knots, "tau": tau, "artefact": new_artefact}
        if change < tol:
            break
    resp, mean_log_s, _ = _e_step(counts, n_fragments, ploidies, weights, theta, kappa)
    per_copy_log_s = np.sum(resp * mean_log_s, axis=1)
    theta = {"coef": theta["coef"].tolist(), "knots": theta["knots"].tolist(), "tau": theta["tau"],
             "artefact": theta["artefact"].tolist(), "weights": weights.tolist()}
    return resp, theta, per_copy_log_s


def call_ploidy(counts_df, group, ploidies=(2, 4, 8), kappa=8.0, min_fragments=1000, max_fragments=100000,
                min_group_cells=200, threshold=0.8, seed=0):
    """Return per-nucleus ploidy posteriors and per-group parameters.

    counts_df has f1..fk columns and an index of cell ids; group aligns with it.
    Nuclei outside the fragment range get missing posteriors.
    """
    columns = [c for c in counts_df.columns if c.startswith("f") and c[1:].isdigit()]
    counts = counts_df[columns].to_numpy(dtype=float)
    n_fragments = counts.sum(axis=1)
    usable = (n_fragments >= min_fragments) & (n_fragments <= max_fragments)
    labels = pd.Series(np.asarray(group, dtype=object), index=counts_df.index).astype(object)
    sizes = labels[usable].value_counts()
    small = sizes[sizes < min_group_cells].index
    labels = labels.where(~labels.isin(small), "pooled_small_groups")

    names = [f"P_{p}C" for p in ploidies]
    output = pd.DataFrame(np.nan, index=counts_df.index, columns=names, dtype=float)
    output.insert(0, "n_fragments", n_fragments.astype(int))
    output["per_copy_log_s"] = np.nan
    parameters = {}
    for name in pd.unique(labels[usable]):
        mask = usable & (labels.to_numpy() == name)
        resp, theta, log_s = fit_group(counts[mask], n_fragments[mask], ploidies, kappa, seed=seed)
        output.loc[output.index[mask], names] = resp
        output.loc[output.index[mask], "per_copy_log_s"] = log_s
        parameters[str(name)] = theta
    output["P_polyploid"] = 1 - output[names[0]]
    best = output[names].fillna(-1).to_numpy()
    calls = np.array([f"{ploidies[i]}C" for i in best.argmax(axis=1)], dtype=object)
    calls[best.max(axis=1) < threshold] = "uncertain"
    calls[~usable] = "not_called"
    output["ploidy_call"] = calls
    return output, parameters
