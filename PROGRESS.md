# Current phase

06 — Closed-pool Partial Observation (`▶ IN PROGRESS`) after verified Phase 04 and Phase 05 exit gates.

# Data location and public/private boundary

Real data and full graph tables remain under `/volume/schen04/ray/0724/`. The Git repository carries code and compact reproducible summaries; ignored raw CSV graph/task tables remain on the shared volume.

# Pipeline

✅ 01 Scope → ✅ 02 Benchmark → ✅ 03 Graph → ✅ 04 Oracle Closeout → ✅ 05 Landscape Characterization → ▶ 06 Partial Observation → ○ 07 Learned Planner → ○ 08 338lib → ○ 09 Paper

# Phase 04 verified result

- 4 assays, 4 proteins, 600 paired task records, 123 `VALLEY_REQUIRED`, 3 seeds/assay, horizon 3, beam width 4.
- Unequal forced-continuation semantics were corrected; all path/edge/horizon/seed/regret/valley audits pass.
- Corrected exact lookahead equals the bounded oracle. The former GCN4 terminal loss was an implementation artifact; old outputs remain archived.
- Raw and normalized effects plus unique-start clustered confidence intervals are frozen.

# Phase 05 verified result

- Assay diagnostics: 4/4; task diagnostics: 600/600.
- WELL_CONNECTED: A4, F7YBW8.
- MODERATELY_CONNECTED: D7PM05.
- SPARSE_NEGATIVE_CONTROL: GCN4 (95.5% isolated; LCC 57/2638).
- Primary Phase 06 assays: A4, D7PM05, F7YBW8. GCN4 remains as a smoke/control case.
- Assay-level relationships are explicitly `EXPLORATORY — N=4 ASSAYS`; no ruggedness-generalization claim is made.

# Scientific interpretation

**SUPPORTED:** Planning advantage exists on some observed protein fitness landscapes under full information.

**NOT ESTABLISHED:** Planning is generally superior across protein landscapes or recoverable when fitness is hidden. Phase 06 now tests the latter under closed-pool partial observation.

# Current gate

The strict Phase 06 method gate has passed; the paired scientific benchmark is not yet complete.

Primary protocol:

- `LINEAGE_WALK`: each query must be an unqueried Hamming-1 neighbor of the current genotype; the queried node becomes current.
- `ADAPTIVE_GREEDY` uses the shared Ridge surrogate with H=1.
- `ADAPTIVE_LOOKAHEAD` uses the same surrogate with H=3, predicted terminal fitness, first-action execution, reveal, retrain, and replan.
- `BUDGETED_WALK_ORACLE` supplies the primary budget-feasible regret reference.
- Strict valley success is computed only from the actual executed lineage.
- Query-zero top-threshold success is retained as `initially_solved`, not counted as policy success.
- Pre-query prediction diagnostics are secondary to optimization metrics.

Eleven deterministic invariant tests pass. A 32-record real-data validation across all four assays passed path, budget, oracle, retraining, and pairing audits.

The earlier 40-record smoke remains frozen as engineering-only `ACTIVE_FRONTIER_SEARCH`. It does **not** establish planning advantage: F7YBW8 surrogate policies collapsed to the same poor outcome and GCN4 saturated rapidly.

Frozen paired benchmark design: 154 unique starts across all four assays, up to 20 tasks per assay/stratum (F7YBW8 retains all 14 available valley starts), 3 policy seeds, 4 methods, and budgets 10/20/50/100, for 7,392 expected records.

Phase 06 remains **IN PROGRESS**. Phase 07 has not started.

# Next 3 actions

1. Run and audit the 7,392-record paired `LINEAGE_WALK` benchmark.
2. Produce paired bootstrap intervals for ADAPTIVE_LOOKAHEAD versus ADAPTIVE_GREEDY in ALL / VALLEY_REQUIRED / NO_VALLEY_REQUIRED strata.
3. Report assay-specific positive and negative cases, surrogate-error diagnostics, and the Phase 07 gate decision without forcing a positive result.

# Blockers

None.
