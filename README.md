# Lung Adenocarcinoma Clinico-Genomic Analysis (TCGA)

How much do tumor mutations and copy-number changes add to clinical factors for predicting survival in lung
adenocarcinoma? Using TCGA PanCancer Atlas data (566 patients) from the public cBioPortal API.

**Status:** in progress. Data download and patient-level feature construction are complete; survival analysis and
multi-omics model comparison are next.

## Quick start
```bash
pip install -r requirements.txt
pip install -e .
pytest
python -m luad.fetch       # download from the cBioPortal API into data/raw/
python -m luad.features    # one row per patient: clinical + genomic features
```
