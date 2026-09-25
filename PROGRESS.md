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

Phase 06 mechanics and methodology-aligned smoke are verified; the full repeated-task/repeated-seed benchmark has not run.

Verified smoke scope:

- 2 deliberately selected valley-required tasks: F7YBW8 primary assay plus GCN4 sparse negative/control.
- 1 seed, 5 paired methods, budgets 10/20/50/100, 40 records total.
- Shared deterministic calibration and shared mutation-feature Ridge for all deployable policies.
- `STATIC_SURROGATE_LOOKAHEAD` is static graph search over the current surrogate; it does not model future observations or retraining.
- `CLAIRVOYANT_REFERENCE` uses hidden pool fitness but is not claimed to be the finite-budget optimum.
- Target attainment and strict observed downhill target-path evidence are separate metrics.
- Top-1%/top-5% success is calibration-adjusted; all 40 smoke records had already reached top-5% during calibration, so policy top-5% success is not informative and is reported as NA.
- Surrogate R²/Spearman are explicitly in-sample observed-training-set diagnostics, not hidden-candidate generalization estimates.
- Raw fitness AUC is retained, with normalized calibration-improvement AUC added for cross-assay scale handling.
- Python syntax checks and 9/9 deterministic tests pass; paired groups, clipping, leakage boundary, idempotent append, and all 24 prior registry rows were verified.

The smoke validates benchmark mechanics only. It does **not** establish a partial-observation planning advantage.

# Next 3 actions

1. Freeze this corrected Phase 06 smoke as its own coherent commit.
2. Run a non-cherry-picked, paired multi-task/multi-seed benchmark across all suitable assays while retaining GCN4 as a control.
3. Add task/seed-clustered paired uncertainty summaries and failure-case analysis before drawing Phase 06 scientific conclusions.

# Blockers

None.
