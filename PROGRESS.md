# Current phase

05 — Cross-assay Landscape Characterization (`▶ IN PROGRESS`) after a verified Phase 04 closeout.

# Data location and public/private boundary

Real data and full graph tables remain under `/volume/schen04/ray/0724/`. The Git repository carries code and compact reproducible summaries; ignored raw CSV graph/task tables remain on the shared volume.

# Pipeline

✅ 01 Scope → ✅ 02 Benchmark → ✅ 03 Graph → ✅ 04 Oracle Closeout → ▶ 05 Landscape Characterization → ○ 06 Partial Observation → ○ 07 Learned Planner → ○ 08 338lib → ○ 09 Paper

# Phase 04 closeout

- 4 assays, 4 proteins, 600 paired task records, 3 seeds/assay, horizon 3, beam width 4.
- Recomputed `VALLEY_REQUIRED`: 123/600. The pasted value 132 is not reproduced.
- Audit found and corrected unequal stopping semantics: greedy/oracle could stop early while beam/lookahead were forced to continue.
- Corrected exact lookahead equals the bounded oracle on all 600 records; old results are retained under `work/results/phase04_precloseout/`.
- All path, edge, horizon, seed sampling, regret-sign, and valley-label audits pass.
- Clustered paired confidence intervals resample unique start states within assay because GCN4 repeats connected starts across seeds.
- Raw regret is retained, with normalized regret and assay-IQR standardized planning gain added.

# Scientific interpretation

**SUPPORTED:** Planning advantage exists on some observed protein fitness landscapes under full information.

**NOT ESTABLISHED:** Planning is generally superior across protein landscapes. The four-assay result is preliminary; GCN4 remains a sparse graph/control case, and oracle performance does not imply recovery under hidden fitness.

# Current gate

Phase 05 must generate consistent assay/task diagnostics, explicit graph-quality classes, exploratory-only N=4 associations, and a defensible Phase 06 assay set.

# Next 3 actions

1. Generate assay-level graph/fitness/task diagnostic metrics with explicit NA reasons.
2. Generate task-level difficulty and valley descriptors from corrected Phase 04 records.
3. Classify assay suitability, close Phase 05, then immediately implement Phase 06 closed-pool partial observation.

# Blockers

None.
