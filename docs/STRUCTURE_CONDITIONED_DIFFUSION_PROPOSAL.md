# Structure-conditioned diffusion for protein fitness landscapes

> Proposal status: research design, not a validated model. The current repository contains biosensor selection data, an AncLacI/mmc2 benchmark, and exploratory structure/landscape analyses; it does not yet contain a trained diffusion model or prospective mutation validation.

## 1. Central hypothesis

Protein mutation experiments are sparse observations of a high-dimensional fitness landscape. A mutation is not an isolated categorical event: its effect depends on sequence context, local three-dimensional environment, conformational state, assay condition, and epistatic interactions with other substitutions.

The project will test whether a **structure-conditioned generative model** can reconstruct a target protein's local fitness/response landscape from sparse mutation measurements, then use posterior samples of that landscape to select mutations for experimental validation.

The central claim is deliberately stronger and more specific than “use diffusion to predict fitness”:

> **Structure-conditioned diffusion can generate uncertainty-aware local fitness landscapes from sparse sequence–structure–activity observations, enabling landscape-level mutation discovery.**

The generated landscape is a posterior simulator/prioritization tool, not a replacement for measurement.

## 2. Scope and terminology

Use two labels separately:

- **Fitness landscape:** a measured or modeled assay score, such as DNA-binding/repression enrichment, growth fitness, stability, or enzyme activity.
- **Response landscape:** a ligand- or condition-specific phenotype, such as biosensor output under a defined bile acid.

The mmc2 LacI/GalR data provide a DNA-binding/repression fitness proxy. They are useful for method development but are not evidence of bile-acid response or transfer to the TetR biosensor. Selection-pool frequency is not itself a paired mutation-response measurement.

For a target protein with wild-type sequence `s_0`, mutation set `M`, structure/context `c`, and assay condition `a`, define a local landscape:

```text
F(M | s_0, structure, context, assay) -> response
```

Because the full combinatorial space is usually unmeasured, the model should represent a distribution:

```text
p(F_unknown | observed mutations, sequence, structure, assay)
```

rather than a single deterministic score.

## 3. System architecture

```text
Public DMS / MAVE / mutational-scanning data
          + target-protein measurements
          + sequence and structure records
                              |
                              v
                 mutation–structure–activity table
                              |
              sequence encoder + structure encoder
                              |
                              v
             masked conditional landscape diffusion
                              |
              posterior landscape samples (multiple)
                              |
        +---------------------+----------------------+
        v                     v                      v
   peaks / basins       uncertainty regions      epistasis
        +---------------------+----------------------+
                              v
       uncertainty-aware candidate and experiment selection
                              |
                              v
                         new measurements
                              |
                              +----> active-learning update
```

Diffusion is not used as a generic structure generator. Its output is a **fitness/response field over a constrained local mutation space**.

## 4. Aim 1 — Build the mutation–structure–activity resource

### 4.1 Public pretraining data

Collect public mutation-response datasets with explicit metadata, prioritizing:

- ProteinGym DMS assays;
- MaveDB MAVE datasets;
- deep mutational scanning and multiplexed functional assays;
- selected multi-mutant combinatorial landscapes such as GB1 and LacI/GalR;
- stability, binding, activity, growth, and fluorescence assays only when their measurement semantics are retained.

The collector must preserve `protein_id`, wild-type sequence, mutant notation, mutation count, assay, raw/normalized score, uncertainty or replicate information, organism, and source accession. Do not pool scores across assays as if they were one universal activity scale.

The initial implementation entry point is:

```text
work/collect_mutation_response_data.py
```

Raw downloads remain outside the public repository. The repository should contain the schema, checksums, source manifest, and reviewed derived summaries—not an uncontrolled mirror of third-party datasets.

### 4.2 Target-project data

For the TetR bile-acid biosensor, the target table should eventually contain:

```text
variant sequence
mutation list
ligand / condition
selection round or assay batch
response measurement
replicate uncertainty
control / WT normalization
```

The critical missing experiment is a paired genotype × ligand × response table. Until it exists, the project can model selection fitness and structural hypotheses, but should not call them ligand-specific response effects.

### 4.3 Structural records

For each wild type and selected variants, combine experimental structures where available with predicted structures from ESMFold/AlphaFold-class tools. The first structural representation should be mutation-sensitive but computationally bounded:

- residue-level coordinates and local frames;
- distance/contact-map changes;
- secondary structure;
- solvent accessibility;
- pocket geometry and residue environment;
- hydrogen-bond/contact changes;
- interface and DNA-contact features;
- confidence scores and model disagreement;
- optional local flexibility proxies.

Do not generate expensive full structures for every public variant in the first pass. Start with wild type, measured target variants, and high-value local neighborhoods.

## 5. Aim 2 — Learn a mutation-sensitive structural representation

Represent a mutation `i: a -> b` using both sequence and structure:

```text
z_mut = [z_sequence(WT), z_structure(WT),
         z_position(i), z_residue(a), z_residue(b),
         delta_local_structure, assay_condition]
```

Candidate encoders:

- ESM2 or a comparable protein language model for sequence context;
- residue-graph or geometric encoder for local 3D neighborhoods;
- explicit delta features from WT versus mutant prediction for prioritized mutations;
- ligand/assay encoder for condition-specific response;
- mutation mask indicating editable positions and allowed amino acids.

The representation must be evaluated by ablation:

1. sequence only;
2. structure only;
3. sequence + structure;
4. sequence + structure + assay condition;
5. sequence + structure + measured sparse observations.

If structure features do not improve held-out landscape reconstruction, they should not be retained merely because they are biologically plausible.

## 6. Aim 3 — Train a masked conditional landscape diffusion model

### 6.1 Landscape object

For a defined protein, assay, and editable mutation set, construct a vector or graph field:

```text
Y = {fitness(m_1), fitness(m_2), ..., fitness(m_n)}
```

For combinatorial data, add pair or higher-order entries and a mutation graph whose edges connect nearby genotypes. The 2D plots in the paper are projections of this high-dimensional object, not the landscape itself.

### 6.2 Masked training task

During training, randomly hide measured mutation values and condition the model on the visible subset:

```text
observed: A->V +0.8, A->F -1.4, A->G -0.1
masked:   A->L ?, A->S ?, ...
```

The diffusion model learns to denoise a noisy landscape while receiving:

- wild-type sequence;
- wild-type structure representation;
- mutation identities and editable-position mask;
- assay/ligand condition;
- observed fitness values and observation mask;
- measurement uncertainty and replicate count;
- optional homolog/evolutionary context.

The output is a set of plausible completed landscapes. Their mean gives a posterior estimate; their spread gives model uncertainty; their cross-position covariance exposes possible epistasis.

### 6.3 Two-stage training

A single model trained on all proteins can learn protein-family or assay artifacts instead of mutation-to-function rules. Use two stages:

**Stage I — general prior**

Train on many proteins and assays after within-assay normalization, with protein and assay identifiers retained. The goal is to learn generic local mutation/structure/function regularities and a calibrated uncertainty prior.

**Stage II — target adaptation**

Adapt or condition the prior using sparse measurements from the target protein. Fine-tuning must be evaluated with protein-level and position-level holdouts, not only random rows.

Diffusion is justified only if it adds value over a simpler conditional imputer or probabilistic regressor. The proposed model must therefore be compared with masked regression, Gaussian process/Bayesian optimization, sequence-only PLM models, and sequence+structure regressors.

## 7. Mutation discovery on the generated landscape

Do not rank only by posterior mean. Use an uncertainty-aware acquisition score, for example:

```text
score(m) = mean[F(m)]
           + beta * std[F(m)]
           + gamma * diversity(m)
           - penalties(m)
```

Penalties include:

- excessive mutation count;
- structural instability or confidence loss;
- distance outside the training distribution;
- synthesis or assay infeasibility;
- predicted expression problems;
- duplicate or near-duplicate candidates.

Select a balanced panel:

1. high posterior mean — exploitation;
2. high uncertainty with plausible upside — exploration;
3. high predicted positive epistasis — interaction testing;
4. random or conservative controls — calibration.

A candidate is a hypothesis for experiment, not a validated mutation.

## 8. Evaluation plan

The central evaluation is landscape reconstruction, not only single-mutation correlation.

### Task A — masked fitness recovery

Hide measured mutations at several observation fractions (for example 10%, 25%, 50%, 75%) and evaluate:

- Spearman/Pearson correlation;
- RMSE or negative log likelihood;
- calibration of prediction intervals;
- top-k recovery of measured high-value variants.

### Task B — generalization splits

Use:

- random variant split;
- position holdout;
- amino-acid substitution holdout;
- sequence-cluster holdout;
- protein holdout;
- LGF → DMS and DMS → LGF transfer only as clearly labeled stress tests.

### Task C — landscape-level metrics

Measure:

- peak and basin recovery;
- rank correlation of local neighborhoods;
- valley/ridge topology preservation;
- pairwise epistasis sign recovery;
- uncertainty calibration near observed versus extrapolative regions.

### Task D — prospective validation

Freeze the model, withhold a genuinely new set of mutations, and test whether selected candidates outperform:

- random selection;
- PLM-only ranking;
- structure-only ranking;
- ordinary supervised regression;
- Bayesian optimization baseline.

The most important prospective metrics are top-k experimental hit rate, enrichment over random, and calibration of uncertainty—not only global regression accuracy.

## 9. Experimental active-learning loop

```text
initial public prior + target sparse measurements
                         |
                         v
             generate posterior landscapes
                         |
                         v
          choose exploitation/exploration/control panel
                         |
                         v
                 paired wet-lab measurements
                         |
                         v
              update target-specific conditioning
                         |
                         +---- repeat
```

For the biosensor, each round should preserve ligand identity and include WT, library baseline, technical replicates, and a defined response readout. This is required to distinguish a mutation effect from pool composition or selection frequency.

## 10. Risks and falsifiable checkpoints

- **Data heterogeneity:** if assay normalization destroys biological meaning, model each assay family separately.
- **Structure noise:** if predicted structural deltas are unstable, use confidence-weighted features and treat structure as a prior, not a label.
- **Sparse combinatorics:** begin with single mutants and selected double-mutant neighborhoods before claiming broad landscape reconstruction.
- **Diffusion overkill:** if masked regression or a Gaussian process matches diffusion on held-out landscapes, diffusion is not justified for that task.
- **Family leakage:** split by protein or sequence cluster before reporting generalization.
- **Transfer failure:** mmc2/LacI success does not establish transfer to TetR or bile-acid response.
- **No paired phenotype:** without genotype × ligand × response measurements, report selection fitness and structural prioritization only.

Go/no-go checkpoints:

1. A clean mutation-response manifest with assay provenance and uncertainty.
2. Structure features improve at least one held-out split without worsening calibration.
3. Masked diffusion beats strong probabilistic baselines on landscape metrics.
4. Candidate ranking improves prospective top-k hit rate over random and PLM baselines.
5. Repeated active-learning rounds show information gain rather than only exploitation.

## 11. Specific aims and expected contribution

**Aim 1 — Resource:** build a provenance-preserving mutation–structure–activity dataset and target-protein sparse-observation format.

**Aim 2 — Model:** train and validate a masked, structure-conditioned diffusion model that generates posterior samples of local fitness landscapes.

**Aim 3 — Discovery:** use uncertainty-aware peak, basin, and epistasis search to nominate mutations and test them in a closed active-learning loop.

The expected contribution is not a generic activity predictor. It is a reproducible framework for **generative reconstruction of protein fitness landscapes from sparse structure–activity observations**, with explicit uncertainty, epistasis analysis, and prospective mutation discovery.

## 12. Current implementation boundary

Already available in this repository:

- mmc2/LacI-GalR data preparation and ruggedness scripts;
- TetR selection and structure analyses;
- exploratory sequence/structure baselines;
- `work/collect_mutation_response_data.py` as the public-data collection entry point.

Not yet completed:

- current public-data download integration and full provenance manifest;
- mutation-sensitive structural feature table;
- masked landscape diffusion implementation;
- held-out landscape benchmark;
- prospective paired genotype × ligand × response validation.

These boundaries are part of the proposal: no diffusion result or mutation-discovery claim should be presented before the corresponding data and validation steps are completed.
