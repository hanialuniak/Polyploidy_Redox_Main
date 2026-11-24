import gzip

import numpy as np
import pandas as pd
import scipy.sparse as sp
import pytest
from scipy.stats import rankdata

from polyredox import doublets, fragments, ploidy, scoring, stats


def brute_force(codes, starts, ends, n, k=6):
    out = np.zeros((n, k), dtype=int)
    unique = sorted(set(zip(codes, starts, ends)))
    for c, s, e in unique:
        depth = sum(1 for c2, s2, e2 in unique if c2 == c and s2 <= s < e2)
        out[c, min(depth, k) - 1] += 1
    return out


def test_overlap_counts_match_brute_force():
    rng = np.random.default_rng(0)
    codes = rng.integers(0, 5, 400)
    starts = rng.integers(0, 600, 400)
    ends = starts + rng.integers(1, 120, 400)
    assert np.array_equal(fragments.count_overlaps(codes, starts, ends, 5), brute_force(codes, starts, ends, 5))


def test_abutting_fragments_do_not_overlap_and_duplicates_collapse():
    counts = fragments.count_overlaps([0, 0, 0], [10, 20, 20], [20, 30, 30], 1)
    assert counts[0, 0] == 2 and counts[0, 1:].sum() == 0


def test_file_parser_matches_arrays_and_excludes(tmp_path):
    path = tmp_path / "f.tsv.gz"
    with gzip.open(path, "wt") as handle:
        handle.write("# comment\n")
        handle.write("chr1\t10\t50\tA\t1\nchr1\t20\t60\tA\t2\nchr1\t30\t70\tZ\t1\n")
        handle.write("chr2\t10\t50\tB\t1\nchrX\t10\t50\tA\t1\nchrX\t11\t50\tA\t1\n")
    counts = fragments.library_overlaps(str(path), np.array(["A", "B"], dtype=object), chunksize=2)
    assert counts[0, :2].tolist() == [1, 1] and counts[1, 0] == 1
    exclude = {"chr1": (np.array([20]), np.array([21]))}
    counts = fragments.library_overlaps(str(path), np.array(["A", "B"], dtype=object), exclude=exclude)
    assert counts[0, :2].tolist() == [1, 0]


def test_unsorted_chromosomes_are_rejected(tmp_path):
    path = tmp_path / "f.tsv"
    path.write_text("chr1\t1\t5\tA\t1\nchr2\t1\t5\tA\t1\nchr1\t9\t15\tA\t1\n")
    with pytest.raises(ValueError):
        fragments.library_overlaps(str(path), np.array(["A"], dtype=object))


def test_scores_match_definition():
    rng = np.random.default_rng(1)
    X = sp.csr_matrix(rng.poisson(0.7, (30, 200)))
    M = 40
    ranks = scoring.rank_matrix(X, M)
    members = [3, 17, 44, 90, 120]
    score = scoring.score_sets(ranks, [members], M)[:, 0]
    dense = X.toarray().astype(float)
    for i in range(30):
        row = dense[i]
        r = rankdata(-row, method="average")
        expressed = row > 0
        capped = np.where(expressed & (r <= M), r, M + 1)[members]
        expected = 1 - (capped.sum() - 5 * 6 / 2) / (5 * M)
        assert score[i] == pytest.approx(expected)


def test_within_donor_beta_matches_explicit_regression():
    rng = np.random.default_rng(2)
    n = 400
    strata = rng.integers(0, 12, n)
    x = rng.random(n)
    z = rng.normal(size=n)
    y = 0.7 * x + 0.3 * z + strata * 0.4 + rng.normal(size=n)
    dummies = pd.get_dummies(strata).to_numpy(float)
    design = np.column_stack([x, z, dummies])
    expected = np.linalg.lstsq(design, y, rcond=None)[0][0]
    result = stats.within_donor_test(y, x, strata, strata, z[:, None], n_flips=99)
    assert result["beta"].iloc[0] == pytest.approx(expected)


def test_benjamini_hochberg_known_values():
    adjusted = stats.benjamini_hochberg([0.01, 0.04, 0.03, np.nan])
    assert adjusted[:3] == pytest.approx([0.03, 0.04, 0.04]) and np.isnan(adjusted[3])


def test_expected_doublet_formula():
    cells = pd.DataFrame({"library": ["a"] * 4, "cell_class": ["x", "x", "x", "y"], "n_recovered": [10000] * 4})
    result = doublets.expected_rates(cells)
    D, p = 0.08, 0.75
    assert result["expected_doublet"].iloc[0] == pytest.approx(D * p / (1 - D + D * p))


def test_ploidy_model_probabilities_are_valid():
    artefact = np.r_[0, 0, 1e-4, 1e-4, 1e-4, 1e-4]
    for p in (2, 4, 8):
        assert np.exp(ploidy._log_categories(p, 8.0, artefact, 6)).sum(axis=1) == pytest.approx(1.0)
