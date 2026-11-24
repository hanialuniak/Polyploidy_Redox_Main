"""Within-donor association tests with donor-level inference.

The exposure varies within donor, so donor x cell type x region strata are absorbed as
fixed effects. Age, sex, PMI and donor-level pathology are removed by construction.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats


def absorb(matrix, codes, n_strata):
    matrix = np.array(matrix, dtype=float)
    single = matrix.ndim == 1
    if single:
        matrix = matrix[:, None]
    counts = np.bincount(codes, minlength=n_strata).astype(float)
    for j in range(matrix.shape[1]):
        means = np.bincount(codes, weights=matrix[:, j], minlength=n_strata) / np.maximum(counts, 1)
        matrix[:, j] -= means[codes]
    return matrix[:, 0] if single else matrix


def _orthonormal(Z):
    if Z is None or Z.shape[1] == 0:
        return None, 0
    U, s, _ = np.linalg.svd(Z, full_matrices=False)
    rank = int(np.sum(s > s.max() * 1e-10))
    return U[:, :rank], rank


def within_donor_test(Y, x, strata, clusters, Z=None, n_flips=9999, seed=0):
    """Test each column of Y against x after absorbing strata and adjusting for Z."""
    Y = np.asarray(Y, dtype=float)
    Y = Y[:, None] if Y.ndim == 1 else Y
    strata_codes, strata_names = pd.factorize(pd.Series(np.asarray(strata, dtype=object)))
    cluster_codes, cluster_names = pd.factorize(pd.Series(np.asarray(clusters, dtype=object)))
    n_strata, n_clusters = len(strata_names), len(cluster_names)
    n = len(x)

    x_tilde = absorb(x, strata_codes, n_strata)
    Y_tilde = absorb(Y, strata_codes, n_strata)
    Q, rank = _orthonormal(absorb(Z, strata_codes, n_strata) if Z is not None else None)
    if Q is not None:
        x_tilde = x_tilde - Q @ (Q.T @ x_tilde)
        Y_tilde = Y_tilde - Q @ (Q.T @ Y_tilde)
    sxx = float(x_tilde @ x_tilde)
    if sxx <= 0:
        raise ValueError("The exposure has no within-stratum variation.")
    beta = x_tilde @ Y_tilde / sxx
    residual = Y_tilde - x_tilde[:, None] * beta
    robust_scores = _cluster_sum(x_tilde[:, None] * residual, cluster_codes, n_clusters)
    null_scores = _cluster_sum(x_tilde[:, None] * Y_tilde, cluster_codes, n_clusters)
    correction = n_clusters / max(n_clusters - 1, 1) * (n - 1) / max(n - n_strata - rank - 1, 1)
    se = np.sqrt(np.sum(robust_scores**2, axis=0) * correction) / sxx
    critical = stats.t.ppf(0.975, max(n_clusters - 1, 1))

    rng = np.random.default_rng(seed)
    flips = rng.choice([-1.0, 1.0], size=(n_flips, n_clusters))
    observed = null_scores.sum(axis=0)
    flipped = flips @ null_scores
    p_flip = (1 + np.sum(np.abs(flipped) >= np.abs(observed) - 1e-12, axis=0)) / (n_flips + 1)
    leverage = _cluster_sum((x_tilde**2)[:, None], cluster_codes, n_clusters).ravel()
    effective = float(leverage.sum() ** 2 / np.sum(leverage**2))
    with np.errstate(divide="ignore", invalid="ignore"):
        t_value = np.where(se > 0, beta / se, 0.0)
    return pd.DataFrame({
        "beta": beta, "se_cluster": se, "ci_low": beta - critical * se, "ci_high": beta + critical * se,
        "t": t_value, "p_signflip": p_flip, "n_nuclei": n, "n_donors": n_clusters,
        "effective_donors": effective,
    })


def _cluster_sum(matrix, codes, n_clusters):
    out = np.zeros((n_clusters, matrix.shape[1]))
    for j in range(matrix.shape[1]):
        out[:, j] = np.bincount(codes, weights=matrix[:, j], minlength=n_clusters)
    return out


def benjamini_hochberg(p_values):
    p = np.asarray(p_values, dtype=float)
    out = np.full(p.shape, np.nan)
    valid = np.flatnonzero(~np.isnan(p))
    if valid.size == 0:
        return out
    order = valid[np.argsort(p[valid])]
    adjusted = p[order] * valid.size / np.arange(1, valid.size + 1)
    out[order] = np.minimum(np.minimum.accumulate(adjusted[::-1])[::-1], 1.0)
    return out


def technical_covariates(frame, columns, degree=3):
    """Return standardised polynomial terms for log-scaled technical columns."""
    terms = []
    for column in columns:
        value = np.log1p(frame[column].to_numpy(float))
        value = (value - value.mean()) / max(value.std(), 1e-9)
        terms.extend(value**power for power in range(1, degree + 1))
    return np.column_stack(terms)
