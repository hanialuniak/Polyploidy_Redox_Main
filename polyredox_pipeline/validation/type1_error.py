"""False-positive rate of the within-donor test when polyploidy has no effect.

Donors differ in size, polyploid rate and noise level, which is what makes cluster
inference fragile. Run: PYTHONPATH=. python3 validation/type1_error.py
"""
import numpy as np
from polyredox import stats

rng = np.random.default_rng(7)
for n_donors in (8, 12, 24):
    for rate in (0.02, 0.08):
        rejections_flip, rejections_t, reps = 0, 0, 1000
        for _ in range(reps):
            sizes = rng.integers(50, 400, n_donors)
            donor = np.repeat(np.arange(n_donors), sizes)
            donor_rate = rng.beta(2, 2 / rate - 2, n_donors)
            truth = rng.random(donor.size) < donor_rate[donor]
            p = np.clip(np.where(truth, rng.beta(8, 2, donor.size), rng.beta(1, 30, donor.size)), 0, 1)
            noise = rng.lognormal(0, 0.6, n_donors)[donor]
            depth = rng.normal(size=donor.size) + 0.5 * truth
            y = rng.normal(size=donor.size) * noise + 0.3 * depth + rng.normal(0, 1, n_donors)[donor]
            result = stats.within_donor_test(y, p, donor, donor, np.column_stack([depth]), n_flips=999,
                                             seed=int(rng.integers(1e9)))
            rejections_flip += result["p_signflip"].iloc[0] < 0.05
            from scipy import stats as st
            rejections_t += 2 * st.t.sf(abs(result["t"].iloc[0]), n_donors - 1) < 0.05
        print(f"donors={n_donors:>2} polyploid~{rate:.2f}  sign-flip: {rejections_flip / reps:.3f}  cluster-t: {rejections_t / reps:.3f}")
