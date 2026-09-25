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

Implement one shared mutation-feature Ridge surrogate, paired calibration and budgets, hidden-fitness-safe policies, registry append, deterministic checks, and a multi-assay smoke test.

# Next 3 actions

1. Implement shared mutation-feature Ridge and hidden-fitness-safe query policies.
2. Add paired iterative query metrics and budgets 10/20/50/100 with clipping.
3. Run deterministic tests and a smoke benchmark including a primary assay and GCN4 control.

# Blockers

None.
