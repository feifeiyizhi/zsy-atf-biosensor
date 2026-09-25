# MASTER_PLAN

## Research question
Can multi-step planning outperform greedy mutation selection for protein fitness/activity optimization on real experimental fitness landscapes?

## Primary endpoint
Given the same experimental/query budget, does multi-step planning find higher-fitness protein variants than greedy search?

## Canonical pipeline and exit gates

### 01 SCOPE FREEZE — ✅ DONE
Exit criteria:
- [x] Protein fitness/activity and sequential mutation optimization are the primary scope.
- [x] Ligand-specific response and bile-acid prediction are secondary/archive topics.
- [x] Valley-required definition is explicit and non-contradictory.

### 02 FITNESS BENCHMARK / DATASET ELIGIBILITY — ✅ DONE
Exit criteria:
- [x] `benchmark_manifest.csv` exists.
- [x] Eligibility logic is tested on real ProteinGym data.
- [x] Eligible assays are identified.
- [x] Rejected assays retain explicit reasons.

### 03 LANDSCAPE GRAPH CONSTRUCTION — ✅ DONE
Exit criteria:
- [x] Unified mutation-set indexed graph builder works across eligible assays.
- [x] Edges represent observed Hamming-1 mutation steps.
- [x] Nodes, edges, connected components, mutation depths, and diagnostics are generated.
- [x] SPG1/Wu and GFP smoke graphs are verified.

### 04 ORACLE BASELINE VALIDATION — ✅ DONE
Exit criteria:
- [x] Random, greedy, beam, k-step lookahead, and oracle are implemented.
- [x] All methods use identical task constraints, including the option to stop before the maximum horizon.
- [x] Adequately sized task sets across four assays are evaluated (600 paired records).
- [x] `VALLEY_REQUIRED` tasks are represented and independently audited (123/600).
- [x] Repeated seeds and clustered paired confidence intervals are available.
- [x] Failure cases and the pre-closeout forced-continuation defect are documented.
- [x] Corrected results are recorded in `experiment_registry.csv` without overwriting the original rows.
Frozen result: `work/results/PHASE04_ORACLE_REPORT.md` and `phase04_oracle_summary.{csv,json}`.

### 05 CROSS-ASSAY / RUGGEDNESS VALIDATION — ✅ DONE
Exit criteria:
- [x] Multiple proteins/assays are evaluated with the same API.
- [x] Planning advantage is joined to consistent landscape and task descriptors.
- [x] Sparse/negative-control assays are retained rather than cherry-picked.
- [x] Benchmark-quality classes and Phase 06 suitability are explicit.
- [x] All assay-level relationships are labeled exploratory (N=4).
Frozen result: `work/results/PHASE05_LANDSCAPE_REPORT.md`, `assay_landscape_diagnostics.csv`, and `task_difficulty.csv`.

### 06 PARTIAL-OBSERVATION BENCHMARK — ▶ IN PROGRESS
Exit criteria:
- [ ] Fixed query budgets are evaluated.
- [ ] All methods receive identical initial observations and budgets.
- [ ] Optimization metrics, not only prediction metrics, are reported.
- [ ] Repeated seeds and failure cases are recorded.

### 07 LEARNED GRAPHWALKS PLANNER — ○ NOT STARTED
Exit criteria:
- [ ] Simple models are compared before a learned planner.
- [ ] Target/model leakage is prevented.
- [ ] Held-out protein/landscape evaluation is complete.
- [ ] Planner is compared with all fixed baselines.

### 08 338LIB EXTERNAL CASE STUDY — ○ NOT STARTED
Exit criteria:
- [ ] Existing 338lib results are reproduced from code.
- [ ] Same evaluation API is applied where data permit.
- [ ] Missing round-0/lib0 measurements are not invented.

### 09 ABLATIONS / PAPER FIGURES / FINAL STORY — ○ NOT STARTED
Exit criteria:
- [ ] Ablations and null/permutation checks are complete.
- [ ] Experiment registry is the source for figures/tables.
- [ ] Positive and negative datasets are reported.
- [ ] Final claims match evidence strength.

## Stable implementation entry points
- `work/collect_mutation_response_data.py`
- `work/build_fitness_benchmark.py`
- `work/build_fitness_graph.py`
- `work/run_oracle_benchmark.py`

## Scientific definitions
- A graph edge is an observed Hamming-1 mutation between measured genotypes, not an embedding neighbor.
- `VALLEY_REQUIRED` means every non-repeating path from the start to a superior reachable target within the horizon contains at least one strictly downhill edge.
- Monotonic paths, valley-crossing paths, and globally better trajectories are reported separately.
- Fitness/activity is the primary target; ligand-specific response and bile-acid prediction are archived secondary work.
