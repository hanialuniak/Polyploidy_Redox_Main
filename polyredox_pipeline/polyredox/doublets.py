"""Compare inferred polyploidy with the rate expected from same-class doublets.

Two diploid nuclei in one droplet are indistinguishable from one tetraploid nucleus by
fragment overlap. This module estimates how much of the signal doublets could explain.
"""
from __future__ import annotations

import numpy as np
import pandas as pd


def expected_rates(cells, rate_per_1000=0.008):
    """Add 'expected_doublet' for each nucleus; needs library, cell_class, n_recovered."""
    frame = cells.copy()
    multiplet = np.minimum(frame["n_recovered"].to_numpy(float) / 1000 * rate_per_1000, 0.5)
    share = frame.groupby(["library", "cell_class"], observed=True)["cell_class"].transform("size") / frame.groupby(
        "library", observed=True
    )["cell_class"].transform("size")
    share = share.to_numpy(float)
    frame["expected_doublet"] = multiplet * share / (1 - multiplet + multiplet * share)
    return frame


def excess_table(cells, classes=None, n_boot=1000, seed=0):
    """Summarise observed polyploidy and doublet-independent excess by class."""
    rng = np.random.default_rng(seed)
    frame = cells.dropna(subset=["P_polyploid"])
    if classes is not None:
        frame = frame[frame["cell_class"].isin(classes)]
    rows = []
    for name, group in frame.groupby("cell_class", observed=True):
        by_donor = group.groupby("donor", observed=True).agg(
            observed=("P_polyploid", "sum"), expected=("expected_doublet", "sum"), n=("P_polyploid", "size")
        )
        values = by_donor.to_numpy(float)
        draws = rng.integers(0, len(values), (n_boot, len(values)))
        sums = values[draws].sum(axis=1)
        excess_draws = (sums[:, 0] - sums[:, 1]) / sums[:, 2]
        observed, expected, total = values.sum(axis=0)
        by_library = group.groupby("library", observed=True).agg(
            observed=("P_polyploid", "mean"), expected=("expected_doublet", "mean"), n=("P_polyploid", "size")
        )
        intercept = slope = np.nan
        if len(by_library) >= 5 and by_library["expected"].std() > 0:
            design = np.column_stack([np.ones(len(by_library)), by_library["expected"]])
            root = np.sqrt(by_library["n"].to_numpy(float))
            intercept, slope = np.linalg.lstsq(design * root[:, None], by_library["observed"].to_numpy() * root, rcond=None)[0]
        rate, expected_rate = observed / total, expected / total
        rows.append({
            "cell_class": name, "n_nuclei": int(total), "n_donors": len(values),
            "observed_polyploid_rate": rate, "expected_doublet_rate": expected_rate,
            "excess_rate": rate - expected_rate,
            "excess_ci_low": np.quantile(excess_draws, 0.025),
            "excess_ci_high": np.quantile(excess_draws, 0.975),
            "estimated_purity": max(0.0, 1 - expected_rate / rate) if rate > 0 else np.nan,
            "library_intercept": intercept, "library_slope": slope,
        })
    return pd.DataFrame(rows)
