# Phase 04 Oracle Report

**Evidence level: PRELIMINARY.** This is a paired, full-information oracle benchmark on four observed assays, not a general protein-design result.

## Audit and correction

The closeout audit found that the pre-closeout implementation allowed greedy/oracle to stop early but forced beam/lookahead to continue to the horizon. The registered 4-assay × 3-seed × 50-start design was rerun after giving every optimizing policy the same up-to-horizon stopping action. Pre-closeout files are retained under `work/results/phase04_precloseout/`.

The reproducible label count is **123 VALLEY_REQUIRED records out of 600**. The earlier on-disk and progress count was 123; the pasted value 132 is not reproduced and was not substituted.

Regret remains `oracle terminal fitness - policy terminal fitness`; negative `lookahead - greedy` regret delta favors lookahead. Raw magnitudes must not be compared across assays.

Cross-assay metrics are (1) normalized regret, `regret / (oracle - start)`, only when the oracle improvement exceeds the assay-specific epsilon, and (2) paired terminal gain divided by assay fitness IQR. Epsilon is `max(1e-12, 1e-9 × max(|IQR|, 1))`. Confidence intervals resample unique start states, so repeated GCN4 starts across seeds are clustered.

## Current result table (ALL TASKS)

| Assay | Tasks | Valley | Raw regret delta | 95% clustered CI | IQR-standardized gain | Lookahead better | Ties |
|---|---:|---:|---:|---:|---:|---:|---:|
| A4_HUMAN_Seuma_2022 | 150 | 22 | -0.400154 | [-0.509571, -0.296377] | 0.134892 | 39.3% | 60.7% |
| D7PM05_CLYGR_Somermeyer_2022 | 150 | 36 | -2406.12 | [-3479.81, -1479.56] | 0.0956355 | 49.3% | 50.7% |
| F7YBW8_MESOW_Aakre_2015 | 150 | 14 | -0.0646283 | [-0.10311, -0.0324926] | 5.15532 | 32.7% | 67.3% |
| GCN4_YEAST_Staller_2018 | 150 | 51 | -0.0810632 | [-0.120974, -0.034803] | 0.237046 | 45.3% | 54.7% |

## 1. Does a planning advantage exist at all?

**SUPPORTED:** Planning advantage exists on some observed protein fitness landscapes. Exact bounded lookahead reaches the bounded oracle under the corrected shared constraints; its paired advantage over greedy is positive on tasks where greedy is locally trapped.

## 2. Is it universal?

**NOT ESTABLISHED:** Planning is generally superior across protein landscapes. Only four assays were evaluated, and many tasks are ties. Full-information success also does not establish query-limited success.

## 3. On which assays does it appear?

The clustered paired interval favors lookahead on: A4_HUMAN_Seuma_2022, D7PM05_CLYGR_Somermeyer_2022, F7YBW8_MESOW_Aakre_2015, GCN4_YEAST_Staller_2018. Group-specific VALLEY_REQUIRED and NO_VALLEY_REQUIRED rows are in the CSV/JSON summary.

## 4. What are the important failure cases?

GCN4 remains an important sparse graph/control case: most nodes are isolated and sampled connected starts repeat across seeds. Its previously reported terminal loss was caused by unequal forced-continuation semantics and does not survive the corrected audit. The sparse observability limitation remains real and is carried into Phase 05 classification rather than hidden.

## 5. What is NOT yet established?

This does not establish universal superiority, a ruggedness law, open-world protein design, transfer to unseen proteins, or recovery of oracle value when fitness is hidden. Those questions require Phases 05–07.

## Reproducibility status

All 4 assay audits passed with zero path, graph, horizon, seed-sampling, regret, or valley-label errors.
