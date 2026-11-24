"""Simulate a small dataset with fragment files and run the whole pipeline on it."""
import numpy as np
import pandas as pd

from polyredox import pipeline, simulate


def test_fragment_files_reproduce_in_memory_counts_and_pipeline_runs(tmp_path):
    simulate.simulate(str(tmp_path), "effect", n_donors=6, nuclei_per_library=60, seed=3,
                      per_copy_fragments=1500, write_fragments=True, n_background_genes=600)
    config = pipeline.load_config(str(tmp_path / "config.yaml"))
    config.update(n_matched_sets=9, n_ml_null_draws=1, n_flips=99, min_fragments=200, n_jobs=2)
    result = pipeline.run(config)

    from_files = pd.read_csv(tmp_path / "results" / "fedecs.tsv.gz", sep="\t", index_col=0)
    in_memory = pd.read_csv(tmp_path / "fedecs_simulated.tsv.gz", sep="\t", index_col=0)
    assert np.array_equal(from_files.loc[in_memory.index].to_numpy(), in_memory.to_numpy())
    for name in ("report.md", "redox_verdicts.tsv", "ploidy_calls.tsv.gz", "doublet_excess.tsv",
                 "ploidy_by_depth.tsv", "ml_summary.json"):
        assert (tmp_path / "results" / name).exists()
    assert set(result["verdict"]["verdict"]) <= {
        "supported_association", "no_association", "not_above_matched_gene_sets",
        "not_specific_to_redox", "depth_sensitive"}
