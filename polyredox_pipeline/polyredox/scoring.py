"""Rank-based (UCell-style) gene-set scoring and expression-matched random gene sets."""
from __future__ import annotations

import numpy as np
import scipy.sparse as sp
from scipy.stats import rankdata

EXCLUDE_SET = "CELLCYCLE_EXCLUDE"


def read_gmt(paths):
    sets = {}
    for path in paths if isinstance(paths, (list, tuple)) else [paths]:
        with open(path) as handle:
            for line in handle:
                fields = line.rstrip("\n").split("\t")
                if len(fields) >= 3:
                    sets[fields[0]] = list(dict.fromkeys(g for g in fields[2:] if g))
    return sets


def rank_matrix(X, max_rank=1500):
    """Return a sparse cells x genes matrix of max_rank + 1 - rank for top-ranked genes."""
    X = sp.csr_matrix(X)
    rows, cols, values = [], [], []
    for i in range(X.shape[0]):
        data = X.data[X.indptr[i]:X.indptr[i + 1]]
        genes = X.indices[X.indptr[i]:X.indptr[i + 1]]
        if data.size == 0:
            continue
        ranks = rankdata(-data, method="average")
        keep = ranks <= max_rank
        rows.append(np.full(keep.sum(), i))
        cols.append(genes[keep])
        values.append(max_rank + 1 - ranks[keep])
    if not rows:
        return sp.csr_matrix(X.shape)
    return sp.csr_matrix((np.concatenate(values), (np.concatenate(rows), np.concatenate(cols))), shape=X.shape)


def indicator(gene_index_sets, n_genes):
    rows = np.concatenate([np.asarray(g, dtype=int) for g in gene_index_sets])
    cols = np.concatenate([np.full(len(g), j) for j, g in enumerate(gene_index_sets)])
    return sp.csc_matrix((np.ones(rows.size), (rows, cols)), shape=(n_genes, len(gene_index_sets)))


def score_sets(ranks, gene_index_sets, max_rank=1500):
    """Return a cells x sets array of UCell-style scores."""
    gene_matrix = indicator(gene_index_sets, ranks.shape[1])
    sizes = np.array([len(g) for g in gene_index_sets], dtype=float)
    credit = np.asarray((ranks @ gene_matrix).todense(), dtype=float)
    rank_sum = sizes * (max_rank + 1) - credit
    return 1 - (rank_sum - sizes * (sizes + 1) / 2) / (sizes * max_rank)


def resolve_sets(sets, genes, min_genes=5):
    """Map names to indices, remove cell-cycle genes, and report what was kept."""
    position = {g: i for i, g in enumerate(genes)}
    excluded = {position[g] for g in sets.get(EXCLUDE_SET, []) if g in position}
    resolved, skipped = {}, {}
    for name, members in sets.items():
        if name == EXCLUDE_SET:
            continue
        found = sorted({position[g] for g in members if g in position} - excluded)
        if len(found) >= min_genes:
            resolved[name] = found
        else:
            skipped[name] = len(found)
    return resolved, skipped, excluded


class MatchedSampler:
    """Draw random gene sets matching a programme's expression-bin composition."""

    def __init__(self, X, resolved, excluded, n_bins=20, seed=0):
        X = sp.csr_matrix(X)
        library = np.maximum(np.asarray(X.sum(axis=1)).ravel(), 1)
        normalised = sp.diags(1e4 / library) @ X
        mean = np.asarray(normalised.mean(axis=0)).ravel()
        self.rng = np.random.default_rng(seed)
        detected = mean > 0
        order = np.argsort(np.argsort(mean[detected]))
        self.bins = np.full(X.shape[1], -1)
        self.bins[detected] = (order * n_bins // detected.sum()).astype(int)
        programme_genes = set().union(*resolved.values()) if resolved else set()
        blocked = programme_genes | set(excluded)
        self.pool = {
            b: np.array([g for g in np.flatnonzero(self.bins == b) if g not in blocked], dtype=int)
            for b in range(n_bins)
        }
        self.fallback = np.array([g for g in np.flatnonzero(detected) if g not in blocked], dtype=int)

    def draw(self, members):
        picked = set()
        for gene in members:
            pool = self.pool.get(int(self.bins[gene]), np.array([], dtype=int))
            if pool.size <= 2:
                pool = self.fallback
            for candidate in self.rng.permutation(pool):
                if int(candidate) not in picked:
                    picked.add(int(candidate))
                    break
        return sorted(picked)
