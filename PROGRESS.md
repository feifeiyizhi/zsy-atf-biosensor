# Current phase

04 — Oracle Baseline Validation (`▶ IN PROGRESS`), with data/structure foundation work running in parallel

# Pipeline

✅ 01 Scope → ✅ 02 Benchmark → ✅ 03 Graph → ▶ 04 Oracle → ○ 05 Cross-assay → ○ 06 Partial observation → ○ 07 Learned planner → ○ 08 338lib → ○ 09 Ablations/paper

# Phase exit criteria

- [x] Random, greedy, beam, k-step lookahead, and oracle implemented.
- [x] Identical graph, starts, horizon, and seed used in smoke comparison.
- [x] Four additional suitable assays evaluated with 3 seeds each.
- [x] 600 total tasks and 123 `VALLEY_REQUIRED` tasks are recorded.
- [ ] Confidence intervals / formal repeated-seed aggregation.
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

PRELIMINARY. Cross-assay signal exists, but one assay shows lookahead underperforming greedy and formal confidence intervals/failure analysis are not complete.

# Current blockers

None. The current gate is statistical aggregation and failure analysis before Phase 04 can be marked done.

# Next 3 actions

1. Compute per-assay repeated-seed confidence intervals and paired task-level comparisons.
2. Inspect GCN4 failures and verify whether the valley definition or lookahead objective explains the negative result.
3. Compute ruggedness descriptors and join planning advantage to `experiment_registry.csv`.
