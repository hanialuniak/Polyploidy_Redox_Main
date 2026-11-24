# Validation on simulated data

Everything here uses simulated data (`polyredox simulate`), so the right answer is known. None of
it shows the method works on real nuclei; that needs sorted 2C/4C nuclei.

## Ploidy calls (Exc neurons, 24 donors, about 6,000 nuclei)

| Depth (fragments per genome copy) | AUC | Mean P(polyploid) vs true rate, lowest depth quartile |
|---|---|---|
| ~4,000 | 0.999 | 0.012 vs 0.007 |
| ~1,500 | 0.967 | 0.060 vs 0.005 (overstated) |

The test data include varied peak accessibility, background fragments, mapping artefacts and
same-class doublets. In an earlier version, the artefact rate was modelled per nucleus instead of
per fragment, and low-depth calls were badly miscalibrated. The model therefore assumes artefacts
scale with fragment count.

## False-positive rate of the within-donor test (1,000 null datasets each, nominal 0.05)

| Donors | Polyploid rate | Sign-flip | Cluster t |
|---|---|---|---|
| 8 | 2% | 0.047 | 0.061 |
| 8 | 8% | 0.054 | 0.055 |
| 12 | 2% | 0.054 | 0.047 |
| 12 | 8% | 0.035 | 0.045 |
| 24 | 2% | 0.047 | 0.054 |
| 24 | 8% | 0.068 | 0.065 |

Monte Carlo error is about ±0.014.

## End-to-end scenarios (`validation/run_scenarios.sh`)

| Scenario | Truth | Pipeline verdicts |
|---|---|---|
| effect | NRF2 and glutathione programmes up in polyploid neurons | Both `supported_association`. Thioredoxin was also supported; it shares a third of its genes with the NRF2 set, which the report flags. Other programmes were not specific or showed no association. |
| null | No ploidy effect; redox tracks per-copy library quality and donor age | All seven programmes significant on the primary test (β ≈ 0.12–0.19 SD), but none passed the specificity checks |
| general | Redox and all other stress programmes up | All `not_specific_to_redox` |
| all | Doublets about a third of the polyploid signal in excitatory neurons | Doublet-expected rate recovered; excess CI above 0 for neurons |

## Failures found and fixed

1. **Adjusting for raw depth created a false association (β ≈ −1.1 SD) under the null.** Ploidy
   raises total fragments, so at equal fragments a tetraploid has lower per-copy quality. The
   pipeline now adjusts for per-copy depth, and thins RNA to equal UMIs instead of regressing on
   UMIs.
2. **Residual bias remains.** Per-copy depth is measured imperfectly, partly because identical
   fragments from different copies collapse in dense peaks. This left β ≈ 0.15 SD under the null,
   which is why verdicts also require specificity and stability across depth adjustments.
3. **The machine-learning gain is not protective on its own.** Under the null, adding redox
   scores raised AP by 0.008. That was the smallest possible empirical p (0.05) with 19 draws,
   because random programmes lowered AP. Read the ML result only alongside the verdicts.
