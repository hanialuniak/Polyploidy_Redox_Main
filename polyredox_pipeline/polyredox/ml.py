"""Does redox state add predictive information about ploidy beyond technical factors?

Posterior ploidy is used as a soft label: each nucleus appears once as positive with
weight P(polyploid) and once as negative with weight 1 - P(polyploid). Splits are by
donor. Redox's gain is compared with the gain from expression-matched random programmes.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.metrics import average_precision_score
from sklearn.model_selection import GroupKFold


def _model(seed):
    return HistGradientBoostingClassifier(
        max_iter=120, learning_rate=0.06, max_leaf_nodes=15, min_samples_leaf=40,
        l2_regularization=1.0, early_stopping=False, random_state=seed,
    )


def _soft(X, p, w):
    return np.vstack([X, X]), np.r_[np.ones(len(p)), np.zeros(len(p))], np.r_[p * w, (1 - p) * w]


def _expected_ap(p, w, score):
    return average_precision_score(
        np.r_[np.ones(len(p)), np.zeros(len(p))], np.r_[score, score], sample_weight=np.r_[p * w, (1 - p) * w]
    )


def _cross_val_scores(X, p, w, folds, seed):
    prediction = np.full(len(p), np.nan)
    for train, test in folds:
        X_train, y_train, w_train = _soft(X[train], p[train], w[train])
        prediction[test] = _model(seed).fit(X_train, y_train, sample_weight=w_train).predict_proba(X[test])[:, 1]
    return prediction


def subsample(p, max_cells, seed):
    """Keep likely polyploids and downsample others, returning indices and weights."""
    rng = np.random.default_rng(seed)
    weights = np.ones(len(p))
    if len(p) <= max_cells:
        return np.arange(len(p)), weights
    positives = np.flatnonzero(p >= 0.5)
    negatives = np.flatnonzero(p < 0.5)
    n_keep = max(max_cells - positives.size, 1000)
    n_keep = min(n_keep, negatives.size)
    kept = rng.choice(negatives, n_keep, replace=False)
    weights[kept] = negatives.size / n_keep
    index = np.r_[positives, kept]
    return index, weights[index]


def evaluate(technical, redox, null_features, p, donors, strata_for_shuffle, feature_names, cohorts=None,
             n_splits=5, max_cells=200000, seed=0):
    """Return summary rows, permutation importance and cross-cohort rows."""
    index, w = subsample(p, max_cells, seed)
    technical, redox, p, donors = technical[index], redox[index], p[index], np.asarray(donors, dtype=object)[index]
    null_features = [n[index] for n in null_features]
    shuffle_groups = np.asarray(strata_for_shuffle, dtype=object)[index]
    n_splits = min(n_splits, len(pd.unique(donors)))
    folds = list(GroupKFold(n_splits=n_splits).split(technical, groups=donors))

    base_prediction = _cross_val_scores(technical, p, w, folds, seed)
    full = np.column_stack([technical, redox])
    full_prediction = _cross_val_scores(full, p, w, folds, seed)
    base_ap = _expected_ap(p, w, base_prediction)
    full_ap = _expected_ap(p, w, full_prediction)
    gains = [
        _expected_ap(p, w, _cross_val_scores(np.column_stack([technical, n]), p, w, folds, seed)) - base_ap
        for n in null_features
    ]
    gains = np.array(gains, dtype=float)
    gain = full_ap - base_ap
    summary = {
        "prevalence": float(np.sum(p * w) / np.sum(w)), "ap_technical": base_ap, "ap_technical_plus_redox": full_ap,
        "ap_gain": gain, "null_gain_mean": float(gains.mean()) if gains.size else np.nan,
        "null_gain_q95": float(np.quantile(gains, 0.95)) if gains.size else np.nan,
        "p_empirical": float((1 + np.sum(gains >= gain)) / (gains.size + 1)) if gains.size else np.nan,
        "n_null_draws": int(gains.size), "n_nuclei_used": int(len(p)),
    }

    rng = np.random.default_rng(seed)
    rows = []
    start = technical.shape[1]
    for j, name in enumerate(feature_names):
        drops = []
        for train, test in folds:
            X_train, y_train, w_train = _soft(full[train], p[train], w[train])
            model = _model(seed).fit(X_train, y_train, sample_weight=w_train)
            intact = _expected_ap(p[test], w[test], model.predict_proba(full[test])[:, 1])
            shuffled = full[test].copy()
            column = shuffled[:, start + j].copy()
            for group in pd.unique(shuffle_groups[test]):
                members = np.flatnonzero(shuffle_groups[test] == group)
                column[members] = rng.permutation(column[members])
            shuffled[:, start + j] = column
            drops.append(intact - _expected_ap(p[test], w[test], model.predict_proba(shuffled)[:, 1]))
        rows.append({"programme": name, "ap_drop_mean": float(np.mean(drops)), "ap_drop_sd": float(np.std(drops))})
    importance = pd.DataFrame(rows).sort_values("ap_drop_mean", ascending=False)

    cross = []
    if cohorts is not None:
        cohorts = np.asarray(cohorts, dtype=object)[index]
        if len(pd.unique(cohorts)) >= 2:
            for cohort in pd.unique(cohorts):
                test, train = np.flatnonzero(cohorts == cohort), np.flatnonzero(cohorts != cohort)
                outcome = {}
                for label, X in (("technical", technical), ("plus_redox", full)):
                    X_train, y_train, w_train = _soft(X[train], p[train], w[train])
                    fitted = _model(seed).fit(X_train, y_train, sample_weight=w_train)
                    outcome[label] = _expected_ap(p[test], w[test], fitted.predict_proba(X[test])[:, 1])
                cross.append({"held_out_cohort": cohort, "ap_technical": outcome["technical"],
                              "ap_technical_plus_redox": outcome["plus_redox"],
                              "ap_gain": outcome["plus_redox"] - outcome["technical"]})
    return summary, importance, pd.DataFrame(cross)
