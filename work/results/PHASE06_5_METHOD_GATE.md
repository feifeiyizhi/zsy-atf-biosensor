# Phase 06.5 Diagnostic Method Gate

## Status

**PASS — ready for the supervised frozen-task run.**

Phase 06 remains **DONE**. Phase 07 remains **HOLD / NOT STARTED**.

## Frozen benchmark reuse

The diagnostic runner reads `phase06_lineage_paired.csv` and requires exactly the existing:

- 154 assay/task IDs and starts;
- start-only calibration sets;
- seeds 1/2/3;
- budgets 10/20/50/100;
- strict no-repeat Hamming-1 `LINEAGE_WALK` graph;
- frozen adaptive H1/H3 records.

It does not sample new tasks or append to the experiment registry.

## True-fitness horizon control

`TRUE_FITNESS_POLICY_DIAGNOSTIC` is evaluation-only. For H1–H5 it:

1. enumerates legal no-repeat continuations under the current observed set;
2. maximizes true terminal fitness with the frozen deterministic tie break;
3. executes only the first action;
4. advances, observes, and replans.

This is the Phase 06 receding-horizon protocol with true fitness substituted for Ridge predictions. It is not the Phase 04 oracle and is not deployable.

Exact search uses branch-and-bound with an admissible relaxed terminal-fitness bound. Every decision records search expansions and exactness. Any capped trajectory is marked non-exact and cannot support an oracle-opportunity conclusion.

## Surrogate horizon sensitivity

H1–H5 use the same Phase 06 mutation-token Ridge, alpha 1.0, calibration, beam width 16, terminal-fitness objective, reveal/retrain/replan loop, and first-action execution.

The four-assay validation reproduced all frozen H1/H3 paths and fitness values across 96 task/budget/seed comparisons.

## Actionable valleys

The original Phase 04 target node is recovered from the frozen `task_id` and audited against start and target fitness. This avoids treating a different equal-fitness node as the original target.

`PHASE06_ACTIONABLE_VALLEY` is budget-specific and requires:

- original `VALLEY_REQUIRED` provenance;
- the original target requires at least one true downhill edge;
- the target is reachable by a simple strict-lineage path;
- the shortest target path fits the evaluated query budget.

The original `VALLEY_REQUIRED` value is never overwritten.

## Local decisions and downhill evidence

Before each action the runner records legal-neighbor count, true/predicted best, ranks, top-1/top-3 quality, local action regret, and local Spearman when defined.

For adaptive H3 downhill actions it additionally records predicted/true downhill, predicted terminal gain, exact true conditional gain within H after the selected first action, actual recovery, queries to recovery, and independent legal recovery-path existence.

## Validation

- Python syntax: PASS.
- Deterministic tests: 12/12 PASS.
- Four assays × one frozen task × H1–H5 × two policy types × four budgets: 160 trajectory records.
- Independent audit: PASS, zero errors.
- TRUE H1–H5 exactness: 160/160 rows exact.
- Frozen H1/H3 reproduction: 96/96 task/budget/seed comparisons.
- Exhaustive-versus-optimized pilot: 160/160 trajectory rows identical in path, budget, fitness, regret, and exactness.

## Full-run parameters

- Tasks: frozen 154.
- Horizons: 1, 2, 3, 4, 5.
- Policies: true-fitness diagnostic and same-Ridge surrogate diagnostic.
- Budgets: 10, 20, 50, 100.
- TRUE search cap: 5,000,000 expansions per decision; any cap is a visible non-exact result and a hard scientific blocker for that task.
- GCN4 simple-path count cap: 5,000,000, with exact/capped status stored separately.
- Expected trajectory snapshots: 6,160.
- Expected frozen H1/H3 reproduction comparisons: 3,696.
