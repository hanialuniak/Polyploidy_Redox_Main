# Route C, mouse hepatocytes — output tables

From `notebooks/RouteC_liver_redox_enzymes.ipynb` run on Supplementary Table S1 of
Lanz et al. 2024, *Nature Structural & Molecular Biology* (PRIDE PXD052786).

| File | Contents |
| --- | --- |
| `redox_enzymes_per_protein.csv` | Each redox enzyme's log2 change in all five comparisons, with its percentile among all proteins |
| `redox_sets.csv` | Per enzyme set: median change, the matched random-set median, p and FDR |

Comparisons in the columns:

- `ploidy_same_dna_conc` — small 4n vs small 2n: twice the DNA, twice the size, same DNA per volume
- `ploidy_same_dna_conc_large` — the same contrast in large cells, an independent repeat
- `ploidy_same_size` — small 4n vs large 2n: same size, twice the DNA per volume
- `size_2n`, `size_4n` — twice the size at constant ploidy, so the genome is diluted
