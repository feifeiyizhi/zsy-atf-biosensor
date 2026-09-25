# Current phase

04 — Oracle Baseline Validation (`▶ IN PROGRESS`), with mutation-sensitive structure and high-quality data foundation work running in parallel.

# Data location and public/private boundary

Real data, PDB collections, full landscapes, and model weights are under `/volume/schen04/ray/0724/`. GitHub is a de-identified public mirror. Only code, schemas, checksummed/reviewed summaries, and small audit artifacts belong in the repository; raw CSVs, PDB collections, full landscapes, credentials, and multi-GB weights stay on the shared volume.

# Structure/data foundation status

- Full-length LacI ESMFold assets exist on the shared volume: 1094 unique sequences, 27 currently visible PDB files in the checked directory, and reviewed full-length analysis artifacts.
- Existing reviewed result: mutation-local structural response is present, but grouped held-out activity gain from full-length ESMFold features is negative/uncertain (`Delta R2 = -0.021`, CI crosses zero). This is a negative predictive result, not a failed engineering run.
- The stronger current LacI feature is WT-template contact perturbation; ESMFold per-variant deltas remain a separate feature family.
- GB1 measured/imputed landscape assets exist on the shared volume and must remain separately labeled.

# Pipeline

✅ 01 Scope → ✅ 02 Benchmark → ✅ 03 Graph → ▶ 04 Oracle → ○ 05 Cross-assay → ○ 06 Partial observation → ○ 07 Learned planner → ○ 08 338lib → ○ 09 Ablations/paper

# Phase exit criteria

- [x] Random, greedy, beam, k-step lookahead, and oracle implemented.
- [x] Identical graph, starts, horizon, and seed used in smoke comparison.
- [x] Four additional suitable assays evaluated with 3 seeds each.
- [x] 600 total tasks and 123 `VALLEY_REQUIRED` tasks are recorded.
- [x] Paired task-level bootstrap confidence intervals computed across 450 tasks.
- [ ] Failure cases inspected in detail.
- [x] `work/results/experiment_registry.csv` populated.

# Latest experiment scale

- assays: 4 additional assays (A4, D7PM05, F7YBW8/Aakre, GCN4)
- proteins: 4
- tasks: 600
- VALLEY_REQUIRED tasks: 123
- seeds: 3 per assay
- horizon: 3
- beam width: 4

# Latest key result

Lookahead mean regret was lower than greedy on A4 (0.041 vs 0.400), D7PM05 (1665 vs 2406), and F7YBW8 (0.016 vs 0.065), but worse on GCN4 (0.109 vs 0.081). The result is heterogeneous rather than universally positive.

# Evidence level

PRELIMINARY. Paired bootstrap intervals now support assay-specific claims: lookahead improves regret over greedy on A4 and F7YBW8; D7PM05 is uncertain overall but favorable on valley tasks; GCN4 is negative overall but favorable on valley tasks. This is not yet a universal planning claim.

# Current blockers

None. Bootstrap aggregation is complete; the gate remains open for failure inspection and a defensible ruggedness/planning analysis.

# Next 3 actions

1. Compute per-assay repeated-seed confidence intervals and paired task-level comparisons.
2. Inspect GCN4 failures and verify whether the valley definition or lookahead objective explains the negative result.
3. Inspect GCN4 failures, then compute/join ruggedness descriptors to `experiment_registry.csv`.
