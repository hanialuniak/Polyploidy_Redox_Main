"""Simulated Multiome data where the right answer is known.

Scenarios
  effect   polyploid neurons up-regulate NRF2 and glutathione programmes only
  null     no ploidy effect, but redox scores track sequencing depth and donor age
  general  polyploid neurons up-regulate redox and all other stress programmes

Every scenario includes same-class doublets, mapping artefacts at depth >= 3, varied
peak accessibility, and donor age that raises both polyploidy and redox.
"""
from __future__ import annotations

import gzip
import os
from importlib import resources

import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp
import yaml

from .fragments import count_overlaps
from .scoring import EXCLUDE_SET, read_gmt

SCENARIOS = ("effect", "null", "general")
CLASSES = {"Exc": (0.55, ["Exc_L23", "Exc_L5"]), "Inh": (0.30, ["Inh_Pvalb", "Inh_Sst"]), "Oligo": (0.15, ["Oligo"])}
PEAK_SPACING, PEAK_WIDTH, N_PEAKS = 5000, 400, 15000


BACKGROUND_BINS = 600000


def _nucleus_fragments(rng, cumulative, copies, per_copy, artefact_clusters, in_peaks):
    """in_peaks is the nucleus's fraction of fragments in peaks; the rest fall in background."""
    starts, ends = [], []
    for _ in range(copies):
        peaked = rng.random(per_copy) < in_peaks
        peak = np.searchsorted(cumulative, rng.random(per_copy))
        start = peak * PEAK_SPACING + rng.integers(0, PEAK_WIDTH, per_copy)
        background = (N_PEAKS + 400 + rng.integers(0, BACKGROUND_BINS, per_copy)) * PEAK_SPACING
        start = np.where(peaked, start, background + rng.integers(0, 4000, per_copy))
        end = start + np.clip(60 + rng.exponential(110, per_copy).astype(np.int64), 50, 900)
        order = np.argsort(start)
        start, end = start[order], end[order]
        keep = np.r_[True, start[1:] >= np.maximum.accumulate(end)[:-1]]
        starts.append(start[keep])
        ends.append(end[keep])
    for _ in range(rng.poisson(artefact_clusters * per_copy * copies / 8000)):
        site = (N_PEAKS + rng.integers(0, 300)) * PEAK_SPACING
        depth = rng.integers(3, 6)
        start = site + rng.integers(0, 40, depth)
        starts.append(start)
        ends.append(start + rng.integers(120, 200, depth))
    return np.concatenate(starts), np.concatenate(ends)


def simulate(outdir, scenario="effect", n_donors=24, nuclei_per_library=250, seed=1, per_copy_fragments=4000,
             write_fragments=False, polyploid_scale=1.0, n_background_genes=2500):
    if scenario not in SCENARIOS:
        raise ValueError(scenario)
    os.makedirs(outdir, exist_ok=True)
    rng = np.random.default_rng(seed)
    programme_sets = read_gmt(str(resources.files("polyredox") / "data" / "starter_genesets.gmt"))
    named = sorted(set().union(*programme_sets.values()))
    genes = named + [f"BG{i:05d}" for i in range(n_background_genes)]
    gene_position = {g: i for i, g in enumerate(genes)}
    base_mean = rng.lognormal(0, 1.3, len(genes))
    class_factor = {c: rng.lognormal(0, 0.35, len(genes)) for c in CLASSES}
    peak_weight = rng.gamma(0.5, 1.0, N_PEAKS)
    class_cumulative = {}
    for name in CLASSES:
        weight = peak_weight * rng.gamma(2.0, 0.5, N_PEAKS)
        class_cumulative[name] = np.cumsum(weight) / weight.sum()

    records, library_rows, count_rows, rna_blocks, fragment_tables = [], [], [], [], []
    class_names = list(CLASSES)
    class_share = np.array([CLASSES[c][0] for c in class_names])
    redox = [n for n in programme_sets if n.startswith("REDOX_")]
    stress = [n for n in programme_sets if n.startswith("CONTROL_")]
    for d in range(n_donors):
        donor, library = f"D{d:03d}", f"L{d:03d}"
        age = rng.uniform(65, 95)
        age_z = (age - 80) / 9
        pathology = float(np.clip(0.5 + 0.2 * age_z + rng.normal(0, 0.2), 0, 1))
        n_recovered = int(rng.uniform(3000, 10000))
        multiplet = n_recovered / 1000 * 0.008
        library_rows.append({"library": library, "n_recovered": n_recovered})
        donor_slope = rng.normal(1.0, 0.4)
        for i in range(nuclei_per_library):
            class_index = rng.choice(len(class_names), p=class_share)
            name = class_names[class_index]
            cell_type = rng.choice(CLASSES[name][1])
            rate = 0.01 if name == "Oligo" else (0.05 + 0.03 * age_z + 0.03 * pathology) * polyploid_scale
            homotypic = multiplet * class_share[class_index] / (1 - multiplet + multiplet * class_share[class_index])
            is_doublet = rng.random() < homotypic
            if is_doublet:
                copies = 4
            elif rng.random() < np.clip(rate, 0, 0.6):
                copies = 8 if rng.random() < 0.05 else 4
            else:
                copies = 2
            depth_z = rng.normal(0, 1)
            per_copy = int(per_copy_fragments * np.exp(0.5 * depth_z))
            cell_id = f"{library}_{i:05d}"
            in_peaks = float(np.clip(rng.beta(6, 7), 0.1, 0.9))
            starts, ends = _nucleus_fragments(rng, class_cumulative[name], copies, per_copy, 6, in_peaks)
            count_rows.append(count_overlaps(np.zeros(starts.size, dtype=np.int64), starts, ends, 1)[0])
            if write_fragments:
                fragment_tables.append(pd.DataFrame({
                    "chrom": [f"chr{1 + c}" for c in (starts // PEAK_SPACING) % 22], "start": starts, "end": ends,
                    "barcode": cell_id, "count": 1, "library": library,
                }))

            true_polyploid = copies > 2 and not is_doublet
            effect = np.zeros(len(genes))
            for programme in redox + stress:
                members = [gene_position[g] for g in programme_sets[programme]]
                activity = rng.normal(0, 0.15)
                if programme in redox:
                    activity += 0.3 * age_z
                    if scenario == "null":
                        activity += 0.4 * depth_z
                if true_polyploid and name != "Oligo":
                    response = 0.5 * donor_slope
                    if scenario == "effect" and programme in ("REDOX_NRF2_TARGETS", "REDOX_GLUTATHIONE"):
                        activity += response
                    if scenario == "general":
                        activity += response
                effect[members] += activity
            library_size = 5000 * np.exp(0.3 * depth_z + rng.normal(0, 0.35)) * (copies / 2) ** 0.8
            mean = base_mean * class_factor[name] * np.exp(effect)
            rna_blocks.append(sp.csr_matrix(rng.poisson(library_size * mean / mean.sum())[None, :]))
            records.append({
                "cell_id": cell_id, "library": library, "barcode": cell_id, "donor": donor, "cell_class": name,
                "cell_type": str(cell_type), "region": "MTG", "cohort": "A" if d % 2 == 0 else "B",
                "pathology": pathology, "age": age, "true_copies": copies, "is_doublet": is_doublet,
                "depth_z": depth_z,
            })

    cells = pd.DataFrame(records)
    truth = ["age", "true_copies", "is_doublet", "depth_z"]
    cells[truth + ["cell_id"]].to_csv(os.path.join(outdir, "truth.tsv"), sep="\t", index=False)
    cells.drop(columns=truth).to_csv(os.path.join(outdir, "cells.tsv"), sep="\t", index=False)
    libraries = pd.DataFrame(library_rows)

    columns = [f"f{i}" for i in range(1, 7)]
    pd.DataFrame(np.vstack(count_rows), columns=columns, index=cells["cell_id"]).to_csv(
        os.path.join(outdir, "fedecs_simulated.tsv.gz"), sep="\t"
    )
    if write_fragments:
        all_fragments = pd.concat(fragment_tables, ignore_index=True)
        directory = os.path.join(outdir, "fragments")
        os.makedirs(directory, exist_ok=True)
        paths = {}
        for library, part in all_fragments.groupby("library"):
            part = part.assign(number=part["chrom"].str[3:].astype(int)).sort_values(["number", "start"])
            path = os.path.join(directory, f"{library}.tsv.gz")
            part[["chrom", "start", "end", "barcode", "count"]].to_csv(path, sep="\t", header=False, index=False)
            paths[library] = os.path.abspath(path)
        libraries["fragments"] = libraries["library"].map(paths)
    libraries.to_csv(os.path.join(outdir, "libraries.tsv"), sep="\t", index=False)

    directory = os.path.join(outdir, "rna")
    os.makedirs(directory, exist_ok=True)
    scipy.io.mmwrite(os.path.join(directory, "matrix.mtx"), sp.vstack(rna_blocks).T.tocoo(), field="integer")
    pd.Series(cells["cell_id"]).to_csv(os.path.join(directory, "barcodes.tsv"), header=False, index=False)
    pd.Series(genes).to_csv(os.path.join(directory, "features.tsv"), header=False, index=False)

    config = {
        "outdir": "results", "cells": "cells.tsv", "libraries": "libraries.tsv",
        "rna": {"format": "mtx", "path": "rna"}, "analysis_classes": ["Exc", "Inh"],
        "n_matched_sets": 199, "n_ml_null_draws": 19, "n_flips": 9999,
    }
    if not write_fragments:
        config["fedecs_table"] = "fedecs_simulated.tsv.gz"
    with open(os.path.join(outdir, "config.yaml"), "w") as handle:
        yaml.safe_dump(config, handle)
    return config
