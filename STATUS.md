# STATUS

> Last updated: 2026-09-20  
> The project has credible selection and structural signals. Ligand-specific mutation effects still require paired phenotype validation. GraphWalks is an independent methods track.

## Where we are now

| ESTABLISHED | NOT ESTABLISHED | NEXT |
|---|---|---|
| Q164R is enriched across all six bile-acid selection pools | Cross-bile-acid prediction ability | Paired genotype × ligand phenotype |
| The 17 designed sites are in the LBD | Target-specific mutation prediction | Round-0 / lib0 baseline |
| Pocket ↔ HTH is a long-range structural relationship | Pool frequency is not mutation effect | Prospective validation |
| GraphWalks shows a real but weak valley signal | The recommender has not passed prospective validation | Strong-valley benchmark |

## Project interpretation

The strongest current statement is:

> **Q164R is the strongest shared mutation across six bile-acid selections and is located at the ligand-pocket wall, with structural evidence supporting a role in ligand recognition.**

This is an enrichment and structural observation. It is **not yet a universal response-switch claim**: causal ligand-specific function requires matched genotype × ligand × response measurements.

## Research tracks

### A. Biosensor biology

```text
bile-acid selection
        ↓
shared / ligand-associated mutations
        ↓
structural interpretation
        ↓
paired genotype × ligand phenotype
        ↓
prospective mutation validation
```

**Current bottleneck:** paired phenotype.

### B. GraphWalks method

```text
338lib landscape
        ↓
directed mutation graph
        ↓
fitness valleys
        ↓
greedy vs lookahead
        ↓
verifiable planning benchmark
```

**Current next step:** build the strong-valley benchmark.

## Evidence

See [`docs/EVIDENCE.md`](docs/EVIDENCE.md) for the source files and reproducible outputs behind this status.

## Explicit non-claims

- Cross-bile-acid or cross-protein generalization is not established.
- Target-specific mutation prediction is not established.
- Selection-pool frequency is not a paired activity measurement.
- v4/v5 rankings are exploratory candidate prioritization, not prospective validation.
- The mmc2 LacI/GalR data are a separate landscape-method benchmark, not evidence of transfer to the TetR biosensor.
