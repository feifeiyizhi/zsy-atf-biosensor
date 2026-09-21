# mmc2 Phase 2 ruggedness baseline

Implemented and ran:

```text
work/mmc2_ruggedness.py
```

Method follows the paper's reported definition:

```text
symmetric kNN graph, k = 34
normalized Dirichlet energy = yᵀLy / (N yᵀy)
95% node subsampling, 1000 replicates, seed=42
```

Outputs:

```text
mmc2_ruggedness_nodes.csv
mmc2_ruggedness_summary.csv
mmc2_ruggedness.json
```

## Baseline values

| Dataset/length stratum | N | full energy | subsample mean | SD |
|---|---:|---:|---:|---:|
| DMS / 59 aa | 1122 | 0.0108520 | 0.0108317 | 0.0004850 |
| LGF / 63 aa | 935 | 0.0050404 | 0.0050403 | 0.0001352 |
| LGF / 61 aa | 153 | 0.0215382 | 0.0215240 | 0.0004596 |
| LGF / 62 aa | 45 | 0.0525770 | 0.0525255 | 0.0014420 |

## Important limitation

The source LGF sequences contain lengths 49, 58, 60, 61, 62, and 63 aa. The paper's graph construction is based on one-hot representations of an aligned/truncated DBD sequence space. The current baseline therefore computes only within equal-length strata and is **not yet the exact paper reproduction for LGF**.

The DMS set is a clean 59-aa single-mutant space and is directly comparable to a fixed-length Hamming graph. Exact LGF reproduction requires importing the source alignment or reconstructing the paper's 60-column aligned DBD representation. Until then, LGF stratum values must be treated as diagnostic baselines, not the paper's reported ruggedness number.
