# Phase 06 LINEAGE_WALK Report

## PROJECT PIPELINE

✅ 01 Scope
✅ 02 Benchmark
✅ 03 Graph
✅ 04 Oracle
✅ 05 Landscape Characterization
✅ 06 Partial Observation
○ 07 Learned Planner — **HOLD / NOT STARTED**
○ 08 338lib
○ 09 Paper

## CURRENT PHASE & GATE

Phase 06 is complete under the predefined exit criteria. The Phase 07 scientific gate is **not passed**.

The paired benchmark does not establish a general partial-observation advantage for H=3 adaptive lookahead over H=1 adaptive greedy. It instead establishes strong assay dependence, including a clear negative case where lookahead hurts.

## SESSION DELTA

- Preserved the earlier engineering smoke as secondary `ACTIVE_FRONTIER_SEARCH` evidence.
- Replaced the primary benchmark with strict `LINEAGE_WALK`.
- Implemented receding-horizon `ADAPTIVE_GREEDY` (H=1) and `ADAPTIVE_LOOKAHEAD` (H=3) with the same Ridge surrogate.
- Added `BUDGETED_WALK_ORACLE`, strict executed-path valley evidence, query-zero calibration labels, and pre-query diagnostics.
- Passed 11 deterministic invariant tests and a 32-record four-assay validation.
- Completed and audited 7,392 benchmark records with zero audit errors.
- Produced task-clustered paired bootstrap intervals from 1,848 exact H=1/H=3 pairs.

## METHOD DEFINITION

### Primary protocol

`LINEAGE_WALK` requires every query to be an unqueried Hamming-1 neighbor of the current genotype. The selected node becomes current. Teleporting across previously observed branches is forbidden.

### Paired policies

- `ADAPTIVE_GREEDY`: H=1.
- `ADAPTIVE_LOOKAHEAD`: H=3.

Both use the same binary mutation-token Kernel Ridge surrogate, start-only calibration, features, query budget, task, and policy seed. At every step the surrogate is fitted to observed measurements, the policy plans, executes one action, reveals the true measured fitness, retrains, and replans.

The lookahead objective is exactly predicted terminal fitness. No information-gain reward or hidden-fitness heuristic is used.

### Evaluation references

- `BUDGETED_WALK_ORACLE`: highest true measured fitness reachable within the requested number of legal steps, with a valid shortest path. This supplies primary regret.
- `clairvoyant_global_ceiling`: descriptive only; it is not the primary oracle.

### Strict valley crossing

A success requires a Phase 04 `VALLEY_REQUIRED` task, an actual downhill transition in the executed lineage, and later recovery on that same lineage to the task target fitness. Query-zero success is excluded.

## EXPERIMENT SCALE

- 4 assays, including GCN4 as the sparse negative/control case.
- 154 unique starts.
- 74 `VALLEY_REQUIRED` and 80 `NO_VALLEY_REQUIRED` starts selected before policy execution.
- 3 policy seeds.
- 4 methods.
- Requested budgets 10, 20, 50, and 100; actual budgets retained.
- 7,392 total records.
- 1,848 exact H=1/H=3 paired records.
- 10,000 task-clustered bootstrap iterations after averaging policy seeds within each unique task.

F7YBW8 had only 14 available unique valley starts; all 14 were retained rather than duplicated or discarded.

## ENGINEERING PROGRESS

- 11/11 deterministic method-gate tests passed.
- Full audit: 7,392/7,392 records, 154/154 tasks, 1,848/1,848 groups.
- All executed paths passed start, no-repeat, edge, and budget checks.
- All adaptive records passed retraining-count, horizon, and first-action checks.
- All oracle paths and primary normalized regrets were reproduced.
- All strict valley labels were recomputed from executed fitness paths.
- Prior 64 registry rows remained byte-identical; 7,392 Phase 06 rows were appended.
- 1,056 records were initially solved and retained with explicit labels.
- 5,149 records ended before the requested budget because a legal unqueried continuation was unavailable or the evaluation oracle target was reached by a shorter path.

## SCIENTIFIC PROGRESS

### Primary pooled result

Delta is `ADAPTIVE_LOOKAHEAD - ADAPTIVE_GREEDY`; negative normalized-regret delta favors lookahead.

| Budget | Stratum | Task clusters | Regret delta | 95% paired bootstrap CI |
|---:|---|---:|---:|---:|
| 10 | ALL | 151 | +0.0013 | [-0.0318, +0.0328] |
| 20 | ALL | 151 | -0.0030 | [-0.0390, +0.0324] |
| 50 | ALL | 151 | -0.0158 | [-0.0532, +0.0203] |
| 100 | ALL | 151 | -0.0216 | [-0.0607, +0.0158] |
| 10 | VALLEY_REQUIRED | 74 | +0.0105 | [-0.0353, +0.0522] |
| 20 | VALLEY_REQUIRED | 74 | +0.0272 | [-0.0254, +0.0761] |
| 50 | VALLEY_REQUIRED | 74 | +0.0034 | [-0.0536, +0.0583] |
| 100 | VALLEY_REQUIRED | 74 | -0.0004 | [-0.0603, +0.0577] |
| 10 | NO_VALLEY_REQUIRED | 77 | -0.0074 | [-0.0570, +0.0378] |
| 20 | NO_VALLEY_REQUIRED | 77 | -0.0320 | [-0.0830, +0.0155] |
| 50 | NO_VALLEY_REQUIRED | 77 | -0.0344 | [-0.0837, +0.0122] |
| 100 | NO_VALLEY_REQUIRED | 77 | -0.0420 | [-0.0920, +0.0035] |

Every pooled confidence interval includes zero. The normalized best-gain AUC and strict-valley-success comparisons also fail to show a stable pooled H=3 advantage.

## PAIRED RESULT TABLE

Budget-100 assay-level normalized-regret comparison:

| Assay | Stratum | Task clusters | H=3 − H=1 | 95% paired bootstrap CI | Interpretation |
|---|---|---:|---:|---:|---|
| A4 | ALL | 40 | -0.0826 | [-0.1193, -0.0493] | H=3 better |
| A4 | VALLEY | 20 | -0.0726 | [-0.1358, -0.0177] | H=3 better |
| A4 | NO_VALLEY | 20 | -0.0926 | [-0.1326, -0.0545] | H=3 better |
| D7PM05 | ALL | 40 | +0.1220 | [+0.0415, +0.2014] | H=3 worse |
| D7PM05 | VALLEY | 20 | +0.1957 | [+0.0976, +0.3015] | H=3 worse |
| D7PM05 | NO_VALLEY | 20 | +0.0483 | [-0.0738, +0.1619] | inconclusive |
| F7YBW8 | ALL | 34 | -0.0733 | [-0.1450, -0.0192] | H=3 better |
| F7YBW8 | VALLEY | 14 | -0.1434 | [-0.3043, -0.0190] | H=3 better, smaller stratum |
| F7YBW8 | NO_VALLEY | 20 | -0.0242 | [-0.0432, -0.0053] | H=3 better |
| GCN4 | ALL | 37 | -0.0634 | [-0.1606, +0.0169] | inconclusive / many ties |
| GCN4 | VALLEY | 20 | -0.0241 | [-0.1462, +0.0606] | inconclusive |
| GCN4 | NO_VALLEY | 17 | -0.1095 | [-0.2819, +0.0095] | inconclusive |

D7PM05 is a required negative case: at budget 100, H=3 averaged only 3.25 legal queries versus 12.48 for H=1 and had worse normalized regret. The current terminal-fitness surrogate objective can drive H=3 into short dead ends.

At budget 100, F7YBW8 was the only assay where both adaptive policies consistently used all 100 queries. A4 and GCN4 commonly terminated much earlier; therefore requested-budget comparisons must always retain actual-budget diagnostics.

## INTERPRETATION

### Supported

- H=3 adaptive planning can outperform H=1 on some observed landscapes and budgets, notably A4 and F7YBW8 at larger budgets.
- Planning can also hurt: D7PM05 shows a clear adverse H=3 effect, strongest in valley-required tasks.
- The effect is assay-dependent rather than universal.
- GCN4 remains a sparse saturation/control case with many ties.
- A budgeted walk oracle leaves substantial headroom in several assays, so information/model quality and path-selection behavior remain important bottlenecks.

### Not established

- H=3 does not show a pooled advantage across all four assays.
- The advantage is not concentrated reliably in `VALLEY_REQUIRED` tasks.
- H=3 does not improve strict valley-crossing success overall.
- The current diagnostics do not establish that surrogate error causally explains failures.
- There is no basis for a universal planning-superiority or ruggedness law.

### Phase 07 gate

**HOLD. Do not begin Phase 07 yet.**

The correct next step is not a larger neural model. The current simple surrogate plus terminal-fitness planning objective produces both positive and negative effects and can terminate early in dead ends. Before learned planning, Phase 06 follow-up should isolate whether the bottleneck is prediction quality, terminal-only planning, or lack of dead-end awareness.

## LIMITATIONS

- Four assays are insufficient for broad cross-landscape laws.
- F7YBW8 has only 14 unique valley tasks.
- Many requested budgets are clipped by strict lineage termination; actual budgets differ across policies as a genuine outcome of their selected trajectories.
- The start-only surrogate is intentionally simple and initially data-poor.
- Pre-query RMSE is assay-scale dependent; its exploratory correlations are not causal evidence.
- Bootstrap uncertainty is task-clustered within this four-assay benchmark and does not represent a population of proteins.
- `BUDGETED_WALK_ORACLE` is evaluation-only and uses true hidden fitness.

## NEXT 3 ACTIONS

1. Add a dead-end-aware but still fixed/simple planning objective and compare it against the frozen terminal-fitness H=3 policy without introducing a learned planner.
2. Add controlled surrogate-quality ablations (for example, additional shared initial observations that do not consume current legal neighbors) to separate information failure from planning failure.
3. Re-run only the necessary paired Phase 06 ablations, preserving A4/F7YBW8 positive cases, D7PM05 harm, and GCN4 saturation before reconsidering Phase 07.

## BLOCKERS

No engineering blocker. The blocker is scientific: pooled H=3 benefit and valley-specific benefit are not established, and D7PM05 shows significant harm.
