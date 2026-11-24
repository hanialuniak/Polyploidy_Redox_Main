"""Run the full analysis from a YAML configuration."""
from __future__ import annotations

import json
import logging
import os
import time
from importlib import resources

import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp
import yaml
from joblib import Parallel, delayed

from . import doublets, fragments, ml, ploidy, report, scoring, stats

log = logging.getLogger("polyredox")
REQUIRED = ["cell_id", "library", "barcode", "donor", "cell_class", "cell_type"]
OPTIONAL = ["region", "cohort", "pathology", "doublet_score"]
DEFAULTS = {
    "analysis_classes": None,
    "k": 6,
    "autosome_regex": fragments.HUMAN_AUTOSOMES,
    "exclude_bed": None,
    "fedecs_table": None,
    "n_jobs": 1,
    "ploidies": [2, 4, 8],
    "kappa": 8.0,
    "min_fragments": 1000,
    "max_fragments": 100000,
    "doublet_rate_per_1000": 0.008,
    "doublet_score_max": None,
    "max_rank": 1500,
    "rna_downsample_quantile": 0.1,
    "extra_covariates": [],
    "min_genes_in_set": 5,
    "n_flips": 9999,
    "n_matched_sets": 999,
    "n_ml_null_draws": 49,
    "max_ml_cells": 100000,
    "seed": 0,
}


def load_config(path):
    with open(path) as handle:
        config = yaml.safe_load(handle)
    base = os.path.dirname(os.path.abspath(path))
    merged = {**DEFAULTS, **config}
    for key in ("outdir", "cells", "libraries", "exclude_bed", "fedecs_table"):
        if merged.get(key) and not os.path.isabs(merged[key]):
            merged[key] = os.path.join(base, merged[key])
    rna = dict(merged["rna"])
    if not os.path.isabs(rna["path"]):
        rna["path"] = os.path.join(base, rna["path"])
    merged["rna"] = rna
    merged["genesets"] = [
        p if os.path.isabs(p) else os.path.join(base, p) for p in merged.get("genesets", [])
    ] or [str(resources.files("polyredox") / "data" / "starter_genesets.gmt")]
    return merged


def read_cells(config):
    columns = config.get("columns") or {}
    cells = pd.read_csv(config["cells"], sep="\t", low_memory=False)
    cells = cells.rename(columns={source: name for name, source in columns.items() if source != name})
    missing = [c for c in REQUIRED if c not in cells.columns]
    if missing:
        raise ValueError(f"cells table lacks required columns: {missing}")
    if cells["cell_id"].duplicated().any():
        raise ValueError("cell_id values must be unique")
    for name in ("region", "cohort"):
        if name not in cells.columns:
            cells[name] = "all"
    for name in ("cell_id", "library", "barcode", "donor", "cell_class", "cell_type", "region", "cohort"):
        cells[name] = cells[name].astype(str)
    return cells


def read_libraries(config, cells):
    if config.get("libraries"):
        libraries = pd.read_csv(config["libraries"], sep="\t", dtype={"library": str})
    else:
        libraries = pd.DataFrame({"library": cells["library"].unique()})
    counted = cells.groupby("library").size().rename("n_in_table").reset_index()
    libraries = libraries.merge(counted, on="library", how="right")
    if "n_recovered" not in libraries.columns:
        libraries["n_recovered"] = np.nan
    estimated = libraries["n_recovered"].isna()
    libraries.loc[estimated, "n_recovered"] = libraries.loc[estimated, "n_in_table"]
    libraries["n_recovered_estimated"] = estimated
    return libraries


def _one_library(row, cells, config, exclude):
    subset = cells.loc[cells["library"] == row["library"], ["cell_id", "barcode"]]
    counts = fragments.library_overlaps(
        row["fragments"], subset["barcode"].to_numpy(dtype=object), config["autosome_regex"], exclude, config["k"]
    )
    frame = pd.DataFrame(counts, columns=[f"f{i}" for i in range(1, config["k"] + 1)], index=subset["cell_id"])
    return frame


def overlap_table(config, cells, libraries):
    if config.get("fedecs_table"):
        table = pd.read_csv(config["fedecs_table"], sep="\t", index_col=0)
        table.index = table.index.astype(str)
        return table.reindex(cells["cell_id"]).fillna(0).astype(int)
    if "fragments" not in libraries.columns:
        raise ValueError("Set fedecs_table or give a fragments column in the libraries table.")
    missing_paths = [p for p in libraries["fragments"] if not os.path.exists(str(p))]
    if missing_paths:
        raise FileNotFoundError(f"Missing fragment files, first: {missing_paths[0]}")
    exclude = fragments.load_exclude_regions(config.get("exclude_bed"))
    frames = Parallel(n_jobs=config["n_jobs"])(
        delayed(_one_library)(row, cells, config, exclude) for _, row in libraries.iterrows()
    )
    return pd.concat(frames).reindex(cells["cell_id"])


def load_rna(rna_config, cell_ids):
    """Return (cells x genes CSR raw counts, gene names, aligned cell ids)."""
    if rna_config.get("format", "mtx") == "mtx":
        directory = rna_config["path"]
        def locate(stem):
            for name in (stem, stem + ".gz"):
                if os.path.exists(os.path.join(directory, name)):
                    return os.path.join(directory, name)
            raise FileNotFoundError(os.path.join(directory, stem))
        matrix = sp.csr_matrix(scipy.io.mmread(locate("matrix.mtx")).T)
        barcodes = pd.read_csv(locate("barcodes.tsv"), sep="\t", header=None, dtype=str)[0].to_numpy(dtype=object)
        features = pd.read_csv(locate("features.tsv"), sep="\t", header=None, dtype=str)
        genes = features[1 if features.shape[1] > 1 else 0].to_numpy(dtype=object)
        if features.shape[1] > 2:
            # Cell Ranger ARC matrices also hold ATAC peaks; keep gene expression only.
            expression = (features[2] == "Gene Expression").to_numpy()
            matrix, genes = matrix[:, expression], genes[expression]
    elif rna_config["format"] == "h5ad":
        import anndata

        data = anndata.read_h5ad(rna_config["path"], backed="r")
        keep = data.obs_names.isin(cell_ids)
        data = data[keep].to_memory()
        layer = rna_config.get("layer")
        matrix = sp.csr_matrix(data.layers[layer] if layer else data.X)
        barcodes = data.obs_names.to_numpy(dtype=object)
        symbol_column = rna_config.get("gene_symbol_column")
        genes = (data.var[symbol_column] if symbol_column else data.var_names).to_numpy(dtype=object)
    else:
        raise ValueError(f"Unknown rna format: {rna_config['format']}")
    position = pd.Index(barcodes).get_indexer(cell_ids)
    found = position >= 0
    aligned = matrix[position[found]]
    if not np.allclose(aligned.data, np.round(aligned.data)):
        raise ValueError("RNA matrix is not integer counts; set rna.layer to the raw UMI layer.")
    return aligned, genes, np.asarray(cell_ids, dtype=object)[found]


def thin_counts(matrix, target, seed):
    """Binomially thin each nucleus to exactly about target UMIs; drop nuclei below it.

    Rank-based scores depend on library size, and polyploid nuclei have more RNA.
    Thinning removes that technical path without conditioning on library size.
    """
    rng = np.random.default_rng(seed)
    matrix = sp.csr_matrix(matrix, dtype=np.int64, copy=True)
    totals = np.asarray(matrix.sum(axis=1)).ravel()
    keep = totals >= target
    probability = np.repeat(np.where(totals > 0, target / np.maximum(totals, 1), 0.0), np.diff(matrix.indptr))
    matrix.data = rng.binomial(matrix.data, np.minimum(probability, 1.0))
    matrix.eliminate_zeros()
    return matrix[keep], keep


def zscore_by(values, groups):
    values = np.array(values, dtype=float)
    groups = np.asarray(groups, dtype=object)
    for group in pd.unique(groups):
        mask = groups == group
        mean = values[mask].mean(axis=0)
        sd = values[mask].std(axis=0)
        values[mask] = (values[mask] - mean) / np.where(sd > 0, sd, 1.0)
    return values


def run_association(scores, names, x, strata, donors, Z, n_flips, seed, scope):
    table = stats.within_donor_test(scores, x, strata, donors, Z, n_flips, seed)
    table.insert(0, "programme", names)
    table.insert(0, "scope", scope)
    is_redox = np.array([n.startswith("REDOX_") for n in names])
    table["p_adj"] = np.nan
    table.loc[is_redox, "p_adj"] = stats.benjamini_hochberg(table.loc[is_redox, "p_signflip"].to_numpy())
    table.loc[~is_redox, "p_adj"] = stats.benjamini_hochberg(table.loc[~is_redox, "p_signflip"].to_numpy())
    return table


def specificity_table(config, ranks, resolved, sampler, names, scores, x, frame, Z, strata, primary):
    n_draws = config["n_matched_sets"]
    rows = []
    stress = [i for i, n in enumerate(names) if n.startswith("CONTROL_")]
    donors = frame["donor"].to_numpy(dtype=object)
    classes = frame["cell_class"].to_numpy(dtype=object)
    primary = primary.set_index("programme")
    for index, name in enumerate(names):
        if not name.startswith("REDOX_"):
            continue
        observed_t = float(primary.loc[name, "t"])
        null_t = []
        for start in range(0, n_draws, 50):
            batch = [sampler.draw(resolved[name]) for _ in range(min(50, n_draws - start))]
            matched = zscore_by(scoring.score_sets(ranks, batch, config["max_rank"]), classes)
            null_t.append(stats.within_donor_test(matched, x, strata, donors, Z, n_flips=1, seed=0)["t"].to_numpy())
        null_t = np.concatenate(null_t) if null_t else np.array([])
        p_matched = float((1 + np.sum(np.abs(null_t) >= abs(observed_t))) / (null_t.size + 1))
        worst_p, worst_control = np.nan, None
        if stress:
            contrasts = scores[:, [index]] - scores[:, stress]
            result = stats.within_donor_test(contrasts, x, strata, donors, Z, config["n_flips"], config["seed"] + 1)
            same_direction = np.sign(result["beta"].to_numpy()) == np.sign(primary.loc[name, "beta"])
            p_values = np.where(same_direction, result["p_signflip"].to_numpy(), 1.0)
            position = int(np.argmax(p_values))
            worst_p, worst_control = float(p_values[position]), names[stress[position]]
        rows.append({
            "programme": name, "t_observed": observed_t, "t_matched_mean": float(null_t.mean()) if null_t.size else np.nan,
            "t_matched_sd": float(null_t.std()) if null_t.size else np.nan, "p_matched": p_matched,
            "p_vs_worst_stress_control": worst_p, "worst_stress_control": worst_control,
        })
    return pd.DataFrame(rows)


def shared_genes(resolved, names, threshold=0.2):
    """For each REDOX_ set, list other REDOX_ sets that share at least threshold of its genes."""
    out = {}
    for name in names:
        if not name.startswith("REDOX_"):
            continue
        own = set(resolved[name])
        partners = [other for other in names if other != name and other.startswith("REDOX_")
                    and len(own & set(resolved[other])) / len(own) >= threshold]
        out[name] = ",".join(partners)
    return out


def depth_diagnostic(cells, bins=10):
    """Polyploid rate by fragment-depth decile; check this against expectations."""
    frame = cells.dropna(subset=["P_polyploid"]).copy()
    frame["depth_decile"] = frame.groupby("cell_class", observed=True)["n_fragments"].transform(
        lambda v: pd.qcut(v.rank(method="first"), bins, labels=False) + 1
    )
    return frame.groupby(["cell_class", "depth_decile"], observed=True).agg(
        n=("P_polyploid", "size"), median_fragments=("n_fragments", "median"),
        mean_P_polyploid=("P_polyploid", "mean"), expected_doublet=("expected_doublet", "mean"),
    ).reset_index()


def depth_sensitivity(scores, names, x, strata, donors, Z, per_copy, n_flips, seed):
    """Refit redox programmes under two other depth adjustments.

    No single depth covariate is guaranteed unbiased, because depth is partly caused by
    ploidy. A result is kept only if it holds without adjustment and in the
    interquartile range of per-copy depth.
    """
    redox = [i for i, n in enumerate(names) if n.startswith("REDOX_")]
    low, high = np.quantile(per_copy, [0.25, 0.75])
    middle = (per_copy >= low) & (per_copy <= high)
    unadjusted = stats.within_donor_test(scores[:, redox], x, strata, donors, None, n_flips, seed)
    trimmed = stats.within_donor_test(scores[middle][:, redox], x[middle], strata[middle], donors[middle],
                                      Z[middle], n_flips, seed)
    return pd.DataFrame({
        "programme": [names[i] for i in redox],
        "beta_unadjusted": unadjusted["beta"].to_numpy(), "p_unadjusted": unadjusted["p_signflip"].to_numpy(),
        "beta_depth_iqr": trimmed["beta"].to_numpy(), "p_depth_iqr": trimmed["p_signflip"].to_numpy(),
    })


def classify(row):
    if not (row["p_adj"] < 0.05):
        return "no_association"
    if not (row["p_matched"] < 0.05):
        return "not_above_matched_gene_sets"
    if pd.isna(row["p_vs_worst_stress_control"]) or row["p_vs_worst_stress_control"] >= 0.05:
        return "not_specific_to_redox"
    for beta, p in ((row["beta_unadjusted"], row["p_unadjusted"]), (row["beta_depth_iqr"], row["p_depth_iqr"])):
        if np.sign(beta) != np.sign(row["beta"]) or p >= 0.05:
            return "depth_sensitive"
    return "supported_association"


def run(config):
    outdir = config["outdir"]
    os.makedirs(outdir, exist_ok=True)
    seed = config["seed"]
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s", datefmt="%H:%M:%S")
    started = time.time()
    step = lambda message: log.info("[%5.0fs] %s", time.time() - started, message)
    cells = read_cells(config)
    libraries = read_libraries(config, cells)
    cells = cells.merge(libraries[["library", "n_recovered"]], on="library", how="left")

    step(f"{len(cells):,} nuclei in {len(libraries)} libraries; counting fragment overlaps")
    overlaps = overlap_table(config, cells, libraries)
    step("Calling ploidy")
    overlaps.to_csv(os.path.join(outdir, "fedecs.tsv.gz"), sep="\t")
    posterior, parameters = ploidy.call_ploidy(
        overlaps, cells["cell_class"].to_numpy(dtype=object), tuple(config["ploidies"]), config["kappa"],
        config["min_fragments"], config["max_fragments"], seed=seed,
    )
    with open(os.path.join(outdir, "ploidy_params.json"), "w") as handle:
        json.dump(parameters, handle, indent=1)
    cells = cells.merge(posterior.reset_index(names="cell_id"), on="cell_id", how="left")
    cells = doublets.expected_rates(cells, config["doublet_rate_per_1000"])
    cells.drop(columns=["barcode"]).to_csv(os.path.join(outdir, "ploidy_calls.tsv.gz"), sep="\t", index=False)

    classes = config["analysis_classes"]
    excess = doublets.excess_table(cells, seed=seed)
    excess.to_csv(os.path.join(outdir, "doublet_excess.tsv"), sep="\t", index=False)
    depth_diagnostic(cells).to_csv(os.path.join(outdir, "ploidy_by_depth.tsv"), sep="\t", index=False)

    frame = cells.dropna(subset=["P_polyploid"])
    if classes:
        frame = frame[frame["cell_class"].isin(classes)]
    limit = config.get("doublet_score_max")
    if limit is not None and "doublet_score" in frame.columns:
        frame = frame[frame["doublet_score"] <= limit]
    if frame.empty:
        raise ValueError("No nuclei remain after filtering.")

    step(f"Loading RNA for {len(frame):,} nuclei")
    matrix, genes, kept = load_rna(config["rna"], frame["cell_id"].to_numpy(dtype=object))
    frame = frame.set_index("cell_id").loc[kept].reset_index()
    frame["rna_umis"] = np.asarray(matrix.sum(axis=1)).ravel()
    frame["rna_genes"] = np.asarray((matrix > 0).sum(axis=1)).ravel()
    target = int(np.quantile(frame["rna_umis"], config["rna_downsample_quantile"]))
    matrix, keep = thin_counts(matrix, target, seed)
    step(f"Thinned RNA to {target:,} UMIs per nucleus; dropped {int((~keep).sum()):,} nuclei below it")
    frame = frame.loc[keep].reset_index(drop=True)
    resolved, skipped, excluded = scoring.resolve_sets(
        scoring.read_gmt(config["genesets"]), genes, config["min_genes_in_set"]
    )
    if not any(n.startswith("REDOX_") for n in resolved):
        raise ValueError("No REDOX_ gene set has enough genes in this dataset; check gene symbols.")
    names = list(resolved)
    step("Scoring gene programmes")
    ranks = scoring.rank_matrix(matrix, config["max_rank"])
    raw_scores = scoring.score_sets(ranks, [resolved[n] for n in names], config["max_rank"])
    scores = zscore_by(raw_scores, frame["cell_class"].to_numpy(dtype=object))
    pd.DataFrame(scores, columns=names).assign(cell_id=frame["cell_id"].to_numpy()).to_csv(
        os.path.join(outdir, "programme_scores.tsv.gz"), sep="\t", index=False
    )

    # Total fragments and UMIs are partly caused by ploidy. Adjusting for them conditions on a
    # collider: at equal total fragments, a tetraploid has half the per-copy quality. Use
    # per-copy ATAC depth plus user-supplied ratio metrics such as TSS enrichment or FRiP.
    copies = sum(frame[f"P_{p}C"].to_numpy(float) * p for p in config["ploidies"])
    frame["expected_copies"] = copies
    frame["fragments_per_copy"] = frame["n_fragments"] * 2 / copies
    technical = ["fragments_per_copy"]
    Z = stats.technical_covariates(frame, technical)
    log_s = frame["per_copy_log_s"].to_numpy(float)
    log_s = (log_s - log_s.mean()) / max(log_s.std(), 1e-9)
    Z = np.column_stack([Z, log_s, log_s**2])
    extra = [c for c in config["extra_covariates"] if c in frame.columns]
    if extra:
        values = frame[extra].to_numpy(float)
        Z = np.column_stack([Z, (values - values.mean(axis=0)) / np.maximum(values.std(axis=0), 1e-9)])
    x = frame["P_polyploid"].to_numpy(float)
    donors = frame["donor"].to_numpy(dtype=object)
    strata = (frame["donor"] + "|" + frame["region"] + "|" + frame["cell_type"]).to_numpy(dtype=object)
    step("Within-donor association tests")
    primary = run_association(scores, names, x, strata, donors, Z, config["n_flips"], seed, "all_analysed")
    tables = [primary]
    for name in pd.unique(frame["cell_class"]):
        mask = (frame["cell_class"] == name).to_numpy()
        if mask.sum() > 50 and frame.loc[mask, "donor"].nunique() >= 4:
            tables.append(run_association(scores[mask], names, x[mask], strata[mask], donors[mask], Z[mask],
                                          config["n_flips"], seed, str(name)))
    association = pd.concat(tables, ignore_index=True)

    observed = frame["P_polyploid"].sum()
    expected = frame["expected_doublet"].sum()
    purity = max(0.0, 1 - expected / observed) if observed > 0 else np.nan
    association["beta_purity_adjusted"] = association["beta"] / purity if purity and purity > 0 else np.nan

    step(f"Comparing with {config['n_matched_sets']} expression-matched gene sets per programme")
    sampler = scoring.MatchedSampler(matrix, resolved, excluded, seed=seed)
    specificity = specificity_table(config, ranks, resolved, sampler, names, scores, x, frame, Z, strata, primary)
    sensitivity = depth_sensitivity(scores, names, x, strata, donors, Z, frame["fragments_per_copy"].to_numpy(float),
                                    config["n_flips"], seed)
    verdict = primary.merge(specificity, on="programme").merge(sensitivity, on="programme")
    verdict["verdict"] = verdict.apply(classify, axis=1)
    verdict["shares_20pct_genes_with"] = verdict["programme"].map(shared_genes(resolved, names))
    association.to_csv(os.path.join(outdir, "association.tsv"), sep="\t", index=False)
    verdict.to_csv(os.path.join(outdir, "redox_verdicts.tsv"), sep="\t", index=False)

    interaction = pd.DataFrame()
    if "pathology" in frame.columns and frame["pathology"].notna().all() and frame["pathology"].std() > 0:
        centred = frame["pathology"].to_numpy(float) - frame["pathology"].mean()
        interaction = stats.within_donor_test(
            scores, x * centred, strata, donors, np.column_stack([Z, x]), config["n_flips"], seed
        )
        interaction.insert(0, "programme", names)
        interaction["p_adj"] = stats.benjamini_hochberg(interaction["p_signflip"].to_numpy())
        interaction.to_csv(os.path.join(outdir, "interaction_pathology.tsv"), sep="\t", index=False)

    step(f"Machine learning with {config['n_ml_null_draws']} matched null draws")
    redox_index = [i for i, n in enumerate(names) if n.startswith("REDOX_")]
    redox_names = [names[i] for i in redox_index]
    classes_array = frame["cell_class"].to_numpy(dtype=object)
    log_technical = np.column_stack([np.log1p(frame[technical].to_numpy(float)), frame["per_copy_log_s"].to_numpy(float)])
    if extra:
        log_technical = np.column_stack([log_technical, frame[extra].to_numpy(float)])
    cell_types = pd.get_dummies(frame["cell_type"]).to_numpy(float)
    ml_technical = np.column_stack([log_technical, cell_types])
    null_features = []
    for _ in range(config["n_ml_null_draws"]):
        draw = [sampler.draw(resolved[n]) for n in redox_names]
        null_features.append(zscore_by(scoring.score_sets(ranks, draw, config["max_rank"]), classes_array))
    summary, importance, cross = ml.evaluate(
        ml_technical, scores[:, redox_index], null_features, x, donors, strata, redox_names,
        frame["cohort"].to_numpy(dtype=object), max_cells=config["max_ml_cells"], seed=seed,
    )
    with open(os.path.join(outdir, "ml_summary.json"), "w") as handle:
        json.dump(summary, handle, indent=1)
    importance.to_csv(os.path.join(outdir, "ml_importance.tsv"), sep="\t", index=False)
    cross.to_csv(os.path.join(outdir, "ml_cross_cohort.tsv"), sep="\t", index=False)

    step("Writing report")
    report.write(outdir, config, frame, excess, verdict, interaction, summary, importance, cross, skipped, purity, libraries)
    return {"verdict": verdict, "excess": excess, "ml": summary, "cross": cross, "interaction": interaction,
            "frame": frame, "association": association}
