"""Download TCGA lung adenocarcinoma (PanCancer Atlas) data from the public cBioPortal REST API.

Run:  python -m luad.fetch
Writes to data/raw/:
  clinical_patient.csv   one row per patient (survival, stage, age, sex, smoking, ...)
  clinical_sample.csv    one row per tumor sample (including tumor mutational burden)
  mutations_drivers.csv  somatic mutations in key lung cancer genes
  cna_drivers.csv        amplifications and deep deletions in the same genes

Using the API (instead of a manual download) makes the data step reproducible: anyone can rerun it.
"""
import csv
import json
import time
import urllib.request
from pathlib import Path

API = "https://www.cbioportal.org/api"
STUDY = "luad_tcga_pan_can_atlas_2018"
GENES = ["EGFR", "KRAS", "TP53", "STK11", "KEAP1", "MET", "ALK", "BRAF", "ERBB2",
         "SMARCA4", "PIK3CA", "NF1", "RBM10", "CDKN2A"]
OUT = Path("data/raw")

# Gene expression panel: compact, interpretable signatures rather than all ~20,000 genes (the cohort has fewer
# than 200 deaths, so a few signature scores are far less prone to overfitting than thousands of genes).
SIGNATURES = {
    "proliferation": ["MKI67", "TOP2A", "CCNB1", "CDK1", "BUB1", "AURKA", "PLK1", "MCM2", "FOXM1", "BIRC5"],
    "cytotoxic_immune": ["CD8A", "CD8B", "GZMA", "GZMB", "PRF1", "CXCL9", "CXCL10", "IFNG", "NKG7"],
    "nrf2_targets": ["NQO1", "AKR1C1", "AKR1C3", "GCLM", "TXNRD1", "SRXN1", "G6PD", "ME1"],
    "immune_checkpoint": ["CD274", "PDCD1", "LAG3", "CTLA4"],
}


def call(path, body=None):
    url = f"{API}{path}"
    data = json.dumps(body).encode() if body is not None else None
    req = urllib.request.Request(url, data=data, method="POST" if body is not None else "GET",
                                 headers={"Content-Type": "application/json", "Accept": "application/json",
                                          "User-Agent": "Mozilla/5.0 (lung-cancer-genomics research script)"})
    for attempt in range(3):
        try:
            with urllib.request.urlopen(req, timeout=120) as r:
                return json.load(r)
        except Exception as e:  # retry transient network errors
            if attempt == 2:
                raise RuntimeError(f"Request failed: {url}\n{e}")
            time.sleep(3)


def pivot_clinical(rows, key):
    """The API returns long format (one row per attribute); turn it into one row per patient or sample."""
    table = {}
    for r in rows:
        table.setdefault(r[key], {key: r[key]})[r["clinicalAttributeId"]] = r["value"]
    return list(table.values())


def write_csv(path, rows):
    cols = sorted({k for r in rows for k in r}, key=lambda c: (c not in ("patientId", "sampleId"), c))
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        w.writerows(rows)
    print(f"{path}: {len(rows)} rows, {len(cols)} columns")


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    study = call(f"/studies/{STUDY}")
    print(f"Study: {study['name']} | {study.get('allSampleCount', '?')} samples")

    for kind, key in (("PATIENT", "patientId"), ("SAMPLE", "sampleId")):
        rows = call(f"/studies/{STUDY}/clinical-data?clinicalDataType={kind}&projection=SUMMARY&pageSize=10000000")
        write_csv(OUT / f"clinical_{kind.lower()}.csv", pivot_clinical(rows, key))

    genes = call("/genes/fetch?geneIdType=HUGO_GENE_SYMBOL&projection=SUMMARY", GENES)
    entrez = {g["entrezGeneId"]: g["hugoGeneSymbol"] for g in genes}
    body = {"sampleListId": f"{STUDY}_all", "entrezGeneIds": list(entrez)}

    muts = call(f"/molecular-profiles/{STUDY}_mutations/mutations/fetch?projection=DETAILED", body)
    write_csv(OUT / "mutations_drivers.csv", [{
        "sampleId": m["sampleId"], "patientId": m["patientId"], "gene": entrez.get(m["entrezGeneId"]),
        "proteinChange": m.get("proteinChange"), "mutationType": m.get("mutationType"),
        "variantType": m.get("variantType")} for m in muts])

    fetch_expression()

    cna = call(f"/molecular-profiles/{STUDY}_gistic/discrete-copy-number/fetch"
               f"?discreteCopyNumberEventType=HOMDEL_AND_AMP&projection=SUMMARY", body)
    write_csv(OUT / "cna_drivers.csv", [{
        "sampleId": c["sampleId"], "patientId": c["patientId"], "gene": entrez.get(c["entrezGeneId"]),
        "alteration": c["alteration"]} for c in cna])  # 2 = amplification, -2 = deep deletion


def fetch_expression():
    """RNA-seq expression for the signature genes, from the study's RNA-seq mRNA profile."""
    profiles = call(f"/studies/{STUDY}/molecular-profiles")
    rna = [p["molecularProfileId"] for p in profiles if p.get("molecularAlterationType") == "MRNA_EXPRESSION"
           and "rna_seq" in p["molecularProfileId"] and "Zscores" not in p["molecularProfileId"]]
    if not rna:
        raise RuntimeError("No RNA-seq expression profile found; available: "
                           + ", ".join(p["molecularProfileId"] for p in profiles))
    profile = sorted(rna, key=len)[0]  # the plain expression profile (shortest id), not a derived one
    symbols = sorted({g for genes in SIGNATURES.values() for g in genes})
    genes = call("/genes/fetch?geneIdType=HUGO_GENE_SYMBOL&projection=SUMMARY", symbols)
    entrez = {g["entrezGeneId"]: g["hugoGeneSymbol"] for g in genes}
    missing = set(symbols) - set(entrez.values())
    data = call(f"/molecular-profiles/{profile}/molecular-data/fetch?projection=SUMMARY",
                {"sampleListId": f"{STUDY}_all", "entrezGeneIds": list(entrez)})
    write_csv(OUT / "expression_panel.csv", [{
        "sampleId": d["sampleId"], "patientId": d["patientId"], "gene": entrez.get(d["entrezGeneId"]),
        "value": d["value"]} for d in data])
    print(f"Expression profile: {profile} | genes found {len(entrez)} of {len(symbols)}"
          + (f" | not found: {sorted(missing)}" if missing else ""))


if __name__ == "__main__":
    main()
