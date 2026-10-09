# Model Card: Multi-Omics Survival Models, Lung Adenocarcinoma

## Model details
- **Models:** Cox proportional hazards models (lifelines, small L2 penalty): clinical only (age, sex, stage); plus
  driver mutations and tumor mutational burden; plus four gene-expression signatures. A stage-stratified model as a
  sensitivity analysis.
- **Developer:** Oshin Miranda. **License:** MIT.

## Intended use
- **Intended:** research and education on clinico-genomic survival analysis and multi-omics integration.
- **Out of scope:** prognosis for individual patients or any treatment decision.

## Data
TCGA Lung Adenocarcinoma, PanCancer Atlas, via the public cBioPortal API: 566 patients; 485 with complete clinical,
mutation, and expression data for model comparison (177 deaths). Median follow-up 21.5 months.

## Performance
| Model | Cross-validated C-index (5 x 5-fold) | Apparent C-index |
|---|---|---|
| Clinical | 0.660 | 0.669 |
| + mutations and TMB | 0.670 | 0.697 |
| + expression signatures | 0.673 | 0.707 |

Molecular layers add little beyond stage for overall prediction; STK11 mutation (HR 1.80) and the proliferation
signature (HR 1.45 per SD) are independently prognostic and robust to stage stratification.

## Representation and fairness
- 86.7% European genetic ancestry, 1.6% East Asian. EGFR-driven tumors, far more common in East Asian patients,
  are underrepresented, so related findings may not generalize.
- Subgroup-specific performance was not estimated: the cohort is too small for reliable subgroup C-indices.

## Limitations
- No treatment or smoking data; both confound survival in lung cancer.
- Short follow-up and moderate size limit power; proportional hazards questionable for sex, KRAS, and stage III.
- Expression standardization was done across the whole cohort (unsupervised) before cross-validation.
- Not externally validated in an independent cohort.
