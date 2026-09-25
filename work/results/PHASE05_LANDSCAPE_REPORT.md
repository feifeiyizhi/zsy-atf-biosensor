# Phase 05 Landscape Characterization

**Status: VERIFIED.** Purpose: explain benchmark observability and task structure, not infer a ruggedness law from four assays.

## Assay classification

| Assay | Class | LCC fraction | Isolated fraction | Greedy failure | Valley tasks | Phase 06 primary |
|---|---|---:|---:|---:|---:|---|
| A4_HUMAN_Seuma_2022 | WELL_CONNECTED | 0.963 | 0.037 | 0.393 | 0.147 | yes |
| D7PM05_CLYGR_Somermeyer_2022 | MODERATELY_CONNECTED | 0.491 | 0.509 | 0.493 | 0.240 | yes |
| F7YBW8_MESOW_Aakre_2015 | WELL_CONNECTED | 1.000 | 0.000 | 0.327 | 0.093 | yes |
| GCN4_YEAST_Staller_2018 | SPARSE_NEGATIVE_CONTROL | 0.0216 | 0.955 | 0.453 | 0.340 | no; retain as control |

Definitions are data-driven and frozen in `phase05_diagnostics_metadata.json`. WELL_CONNECTED requires LCC fraction ≥0.90 and isolated fraction ≤0.10; MODERATELY_CONNECTED requires LCC fraction ≥0.25 and isolated fraction ≤0.60; SPARSE_NEGATIVE_CONTROL has edges but LCC fraction <0.10 or isolated fraction ≥0.90. Graphs without usable edges or 20 eligible starts are UNSUITABLE_FOR_LOCAL_WALK.

## Diagnostic interpretation

A4 and F7YBW8 support broad local walks. D7PM05 supports a large connected benchmark component but roughly half its observed nodes are isolated, so component membership must be explicit. GCN4 remains visible: 95.5% of nodes are isolated and its largest component contains only 57/2638 nodes. It is scientifically useful as an observability/control case but is not a primary Phase 06 landscape.

The task table contains all 600 corrected Phase 04 records, including start fitness, bounded oracle target, valley geometry, raw and normalized regret, component size, and start degree. Metrics unsupported by a task are NA with a reason; in particular valley width is NA for non-valley tasks.

## Exploratory relationships

**EXPLORATORY — N=4 ASSAYS.** Descriptive Spearman values are stored in `phase05_exploratory_associations.json`. No p-values or generalization claims are reported. Extreme ranks (including correlations of ±1) are expected with four points and must not be interpreted as evidence that ruggedness predicts planning performance. Task rows are clustered within assays and are not independent protein landscapes.

## Phase 05 exit gate

- Phase 04 results frozen: yes.
- Cross-assay normalized metrics: yes.
- Assay diagnostics: yes.
- Task difficulty diagnostics: yes, 600/600 records.
- Positive and sparse/control assays retained: yes.
- Four-assay interpretation conservative: yes.
- Unsupported ruggedness-generalization claim avoided: yes.
- Phase 06 primary assays explicit: A4, D7PM05, F7YBW8; GCN4 retained for smoke/control analysis.
