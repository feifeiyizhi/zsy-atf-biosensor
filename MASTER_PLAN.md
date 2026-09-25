# MASTER_PLAN

## Research question
Can multi-step planning outperform greedy mutation selection for protein fitness/activity optimization on real experimental fitness landscapes?

## Primary endpoint
Given the same experimental/query budget, does multi-step planning find higher-fitness protein variants than greedy search?

## Pipeline status
- [x] ProteinGym metadata integration
- [x] assay eligibility audit
- [x] benchmark manifest
- [x] graph construction
- [x] landscape diagnostics
- [x] baseline algorithms
- [x] oracle GraphWalks
- [ ] partial-observation GraphWalks
- [ ] cross-protein evaluation
- [ ] 338lib external case study
- [ ] ablations
- [ ] paper-ready figures/tables

## Current stage
Phase C oracle task validation has started. On 50 deterministic SPG1/Wu tasks at horizon 4, the oracle found 12 `VALLEY_REQUIRED` tasks; greedy had positive regret on all 12. Across all tasks, oracle mean terminal fitness was 5.047 versus greedy 3.619 (mean regret 1.429). This is an initial smoke result, not a final claim; larger balanced task samples and other assays are required.

## Completed outputs
- `MASTER_PLAN.md`: frozen operational research plan and valley definition.
- `work/build_fitness_graph.py`: mutation-set indexed observed Hamming-1 graph builder.
- `work/run_oracle_benchmark.py`: fixed-horizon oracle greedy comparison and valley task audit.
- `work/results/benchmarks/SPG1_STRSG_Wu_2016/`: 149,360-node graph and first 50-task oracle run.
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
