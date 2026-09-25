# Phase 06 Method Gate

## Status

**ENGINEERING/METHODOLOGICAL GATE: PASSED**

**SCIENTIFIC GATE: NOT YET PASSED**

Phase 06 remains in progress. Phase 07 must not begin until the paired multi-task/multi-seed benchmark and paired uncertainty analysis are complete.

## Frozen engineering smoke

The previously committed closed-pool smoke is retained as an engineering-only `ACTIVE_FRONTIER_SEARCH` result. It permits querying from the frontier of any observed node and therefore is not the primary trajectory benchmark and is not used for valley-crossing claims.

Its existing static planning policy is `STATIC_SURROGATE_LOOKAHEAD`. Its hidden-fitness reference is interpreted as a `CLAIRVOYANT_GLOBAL_REFERENCE` / descriptive global ceiling, not a budget-feasible walk oracle.

The smoke demonstrates execution only. It does not demonstrate a partial-observation planning advantage. In particular:

- F7YBW8: surrogate-guided methods collapsed to the same poor best value in the smoke.
- GCN4: methods saturated rapidly.
- The smoke contained only two tasks and one seed.
- Negative results are retained in the committed smoke artifacts.

## Primary protocol: LINEAGE_WALK

The primary benchmark now uses a strict sequential lineage:

1. Start at `s0`.
2. Select only an unqueried Hamming-1 neighbor of the current genotype.
3. Commit the query before revealing measured fitness.
4. Move current state to the queried genotype.
5. Refit the shared Ridge surrogate after every reveal.
6. Continue until budget exhaustion or no legal move remains.

Every consecutive pair in `executed_path` must be a graph edge. Teleporting to another observed branch is forbidden.

Formal implementation: `work/run_lineage_walk.py`.

## Paired methods

- `RANDOM_WALK`
- `ADAPTIVE_GREEDY`: receding horizon `H=1`
- `ADAPTIVE_LOOKAHEAD`: receding horizon `H=3`
- `BUDGETED_WALK_ORACLE`: evaluation-only true-fitness oracle

`ADAPTIVE_GREEDY` and `ADAPTIVE_LOOKAHEAD` use the same binary mutation-token Kernel Ridge model, feature representation, calibration, task, seed, and query budget. Their intended methodological difference is planning horizon.

The lookahead objective is exactly **predicted terminal fitness**. Beam search limits computation; only the first planned action is executed. The true fitness is then revealed, appended, and the surrogate is retrained before replanning.

Formal benchmark calibration is start-only. This preserves all unqueried one-step actions and aligns the deployable and oracle action spaces.

## Oracle

For each requested budget `B`, `BUDGETED_WALK_ORACLE` finds the highest measured fitness reachable from the start within at most `B` legal graph steps and stores a valid shortest path to that node.

The global component maximum is retained only as `clairvoyant_global_ceiling`; it is not used for primary regret.

## Strict valley evidence

A strict valley crossing requires:

- a Phase 04 `VALLEY_REQUIRED` task;
- an actual downhill edge in the executed lineage;
- later recovery on the same lineage to the task target fitness;
- the target not already satisfied at query zero.

Stored evidence includes executed node/fitness paths, downhill count, maximum drop, minimum path fitness, first downhill query, recovery query, and strict success.

## Calibration and diagnostics

Records store initial best fitness, top-1/top-5 status, and `initially_solved`. Thresholds reached during calibration are not counted as policy query success.

Before each reveal, deployable surrogate policies record predicted and true fitness, prediction error, action rank, legal candidate count, planned path, and training observation count. Run-level diagnostics include pre-query RMSE, Spearman where supported, and surrogate top-choice ranking accuracy. These are secondary diagnostics.

## Verified invariants

Eleven deterministic invariant tests pass:

1. every executed transition is a graph edge;
2. no branch teleporting;
3. deployable planners do not accept hidden fitness;
4. Ridge training set grows after every reveal;
5. greedy uses `H=1`;
6. lookahead uses `H=3` and executes only its first action;
7. budgeted oracle respects the walk budget;
8. strict valley success requires executed downhill-and-recovery evidence;
9. query-zero threshold success is labelled initially solved;
10. paired policies share calibration;
11. calibration does not consume current legal neighbors.

A real-data validation across all four assays produced 32 records (4 tasks × 1 seed × 4 methods × 2 budgets). Path, budget, oracle, retraining, and pairing audits all passed. This validation is not a scientific result.

## Frozen paired benchmark design

- Assays: all four existing assays, including GCN4.
- Task sampling: deterministic unique-start stratification across the Phase 04 task pool.
- Target: up to 20 tasks per assay per stratum.
- Available sample: 154 unique starts.
  - A4: 20 valley + 20 non-valley.
  - D7PM05: 20 valley + 20 non-valley.
  - F7YBW8: all 14 available valley + 20 non-valley.
  - GCN4: 20 valley + 20 non-valley.
- Policy seeds: 1, 2, 3.
- Budgets: 10, 20, 50, 100, with actual budget recorded.
- Expected records: 7,392.

Primary paired comparison: `ADAPTIVE_LOOKAHEAD - ADAPTIVE_GREEDY`, reported for ALL, VALLEY_REQUIRED, and NO_VALLEY_REQUIRED with paired bootstrap confidence intervals.
