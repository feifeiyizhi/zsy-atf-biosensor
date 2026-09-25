# MASTER_PLAN

## Research question
Can multi-step planning outperform greedy mutation selection for protein fitness/activity optimization on real experimental fitness landscapes?

## Primary endpoint
Given the same experimental/query budget, does multi-step planning find higher-fitness protein variants than greedy search?

## Pipeline status
- [x] ProteinGym metadata integration
- [ ] assay eligibility audit
- [ ] benchmark manifest
- [ ] graph construction
- [ ] landscape diagnostics
- [ ] baseline algorithms
- [ ] oracle GraphWalks
- [ ] partial-observation GraphWalks
- [ ] cross-protein evaluation
- [ ] 338lib external case study
- [ ] ablations
- [ ] paper-ready figures/tables

## Current stage
Phase A assay audit completed on a 10,000-row real ProteinGym smoke manifest. The compiler found two assays and classified both as `PREDICTION_ONLY`: no observed multi-mutants, no Hamming-1 nested edges, and no depth-2 paths in this prefix. The next run must use the full archive before concluding that no ProteinGym assay is eligible.

## Completed outputs
- `MASTER_PLAN.md`: frozen operational research plan and valley definition.
- `work/build_fitness_benchmark.py`: assay-level compiler with explicit GraphWalk eligibility and diagnostics.
- `work/collect_mutation_response_data.py`: provenance-preserving row collector with explicit missingness and protein-level splits.
- `work/results/evolution_graph.json`: prior 338lib GraphWalks diagnostic (secondary case study; not yet on the common evaluation API).
- `work/results/mmc2_landscape/`: existing LacI/GalR fitness benchmark outputs.

## Blockers
- ProteinGym substitution data are mostly single-mutant assays; GraphWalks eligibility must be measured assay by assay.
- Processed ProteinGym tables do not provide response uncertainty/replicate counts; missing values remain explicit.

## Next automatic action
Run `python3 work/build_fitness_benchmark.py` against the cached ProteinGym archive, audit eligibility, then build graphs and run oracle greedy/lookahead comparisons only for eligible assays. Do not promote prediction-only assays into the GraphWalk benchmark.

## Scientific definitions
- A graph edge is an observed Hamming-1 mutation between measured genotypes, not an embedding neighbor.
- `VALLEY_REQUIRED` means every path from the start to a superior reachable target contains at least one strictly downhill edge; an existing downhill neighbor alone is insufficient.
- Monotonic paths, valley-crossing paths, and globally better trajectories are reported separately.
- Fitness/activity is the primary target; ligand-specific response and bile-acid prediction are archived secondary work.
