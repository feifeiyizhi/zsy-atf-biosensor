# Current phase

06 — Closed-pool Partial Observation (`✅ DONE`). Phase 07 is `HOLD / NOT STARTED` after the scientific gate.

# Data location and public/private boundary

Real data and full graph tables remain under `/volume/schen04/ray/0724/`. The Git repository carries code and compact reproducible summaries; ignored raw CSV graph/task tables remain on the shared volume.

# Pipeline

✅ 01 Scope → ✅ 02 Benchmark → ✅ 03 Graph → ✅ 04 Oracle Closeout → ✅ 05 Landscape Characterization → ✅ 06 Partial Observation → ○ 07 Learned Planner (`HOLD`) → ○ 08 338lib → ○ 09 Paper

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

**SUPPORTED:** Under strict hidden-fitness `LINEAGE_WALK`, H=3 adaptive lookahead helps on some assays and budgets, notably A4 and F7YBW8 at budget 100.

**NOT ESTABLISHED:** H=3 planning is generally superior, is specifically better on valley-required tasks, or improves strict valley crossing. D7PM05 is a significant harm case and GCN4 is largely a sparse saturation/control case.

# Phase 06 closeout

The primary paired benchmark is complete:

- 154 unique starts: 74 `VALLEY_REQUIRED`, 80 `NO_VALLEY_REQUIRED`.
- 3 policy seeds, 4 methods, requested budgets 10/20/50/100.
- 7,392 records and 1,848 exact H=1/H=3 pairs.
- 10,000 task-clustered bootstrap iterations after averaging policy seeds within unique tasks.
- 0 audit errors across path, edge, no-repeat, budget, oracle, retraining, horizon, first-action, valley, calibration-pairing, and registry checks.
- The prior 64 registry rows remain byte-identical; 7,392 strict-lineage rows were appended.

Pooled H=3 minus H=1 normalized-regret deltas were +0.0013, -0.0030, -0.0158, and -0.0216 at budgets 10, 20, 50, and 100. Every pooled 95% paired bootstrap interval includes zero. At budget 100, A4 and F7YBW8 favor H=3, D7PM05 significantly favors H=1, and GCN4 is inconclusive.

The earlier 40-record smoke remains frozen as engineering-only `ACTIVE_FRONTIER_SEARCH` and is not used for valley-crossing claims. Full paired raw CSV/JSON remain on the shared volume; compact audit, analysis, and report artifacts are versioned.

Phase 06 is **DONE**. Phase 07 is **HOLD / NOT STARTED** because pooled and valley-specific H=3 benefits are not established and D7PM05 shows significant harm.

# Next 3 actions

1. Test a fixed dead-end-aware objective against the frozen terminal-fitness H=3 policy without introducing a learned planner.
2. Run controlled shared-calibration/surrogate-quality ablations that preserve current legal neighbors.
3. Reconsider the Phase 07 gate only after paired Phase 06 ablations preserve both positive and negative assay results.

# Blockers

None.
