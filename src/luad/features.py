"""Build one analysis-ready row per patient from the raw cBioPortal files.

Run:  python -m luad.features
Writes data/processed/patients.csv and, in results/: a cohort summary (Table 1), driver alteration frequencies,
a data-quality report, and a genetic-ancestry representation check.
"""
import re
from pathlib import Path

import pandas as pd

RAW, PROCESSED, RESULTS = Path("data/raw"), Path("data/processed"), Path("results")

# Mutation classes that change the protein. Silent and non-coding changes are excluded.
PROTEIN_ALTERING = {"Missense_Mutation", "Nonsense_Mutation", "Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site",
                    "In_Frame_Del", "In_Frame_Ins", "Translation_Start_Site", "Nonstop_Mutation"}
DRIVER_GENES = ["TP53", "KRAS", "EGFR", "STK11", "KEAP1", "NF1", "BRAF", "SMARCA4", "RBM10", "PIK3CA", "ERBB2",
                "MET", "ALK", "CDKN2A"]
STAGE_GROUPS = {"I": "I", "II": "II", "III": "III", "IV": "IV"}


def stage_group(raw):
    """'STAGE IIIA' -> 'III'. Unknown or missing values stay missing (never guessed)."""
    if not isinstance(raw, str):
        return None
    m = re.match(r"STAGE\s+(IV|III|II|I)", raw.strip().upper())
    return STAGE_GROUPS[m.group(1)] if m else None


def event_flag(status):
    """cBioPortal codes survival status as '1:DECEASED' / '0:LIVING' (or '1:Recurred/Progressed')."""
    if not isinstance(status, str) or ":" not in status:
        return None
    return int(status.split(":")[0])


def egfr_classic_activating(protein_change, mutation_type):
    """Classic EGFR activating mutations: L858R, or an exon 19 in-frame deletion (starting at residues 744-753)."""
    if protein_change == "L858R":
        return True
    m = re.match(r"^[A-Z](\d+)_", str(protein_change))
    return mutation_type == "In_Frame_Del" and m is not None and 744 <= int(m.group(1)) <= 753


def build(raw=RAW):
    pat = pd.read_csv(raw / "clinical_patient.csv")
    smp = pd.read_csv(raw / "clinical_sample.csv")
    mut = pd.read_csv(raw / "mutations_drivers.csv")
    cna = pd.read_csv(raw / "cna_drivers.csv")

    df = pd.DataFrame({
        "patient_id": pat.patientId,
        "age": pd.to_numeric(pat.AGE, errors="coerce"),
        "sex": pat.SEX.str.title(),
        "stage": pat.AJCC_PATHOLOGIC_TUMOR_STAGE.map(stage_group),
        "genetic_ancestry": pat.GENETIC_ANCESTRY_LABEL,
        "os_months": pd.to_numeric(pat.OS_MONTHS, errors="coerce"),
        "os_event": pat.OS_STATUS.map(event_flag),
        "pfs_months": pd.to_numeric(pat.PFS_MONTHS, errors="coerce"),
        "pfs_event": pat.PFS_STATUS.map(event_flag),
    })

    # Tumor-level genomic summaries (one primary tumor sample per patient; sample IDs ending in -01)
    smp = smp[smp.sampleId.str.endswith("-01")].copy()
    smp["patient_id"] = smp.sampleId.str[:12]
    smp = smp.rename(columns={"TMB_NONSYNONYMOUS": "tmb", "FRACTION_GENOME_ALTERED": "fraction_genome_altered",
                              "ANEUPLOIDY_SCORE": "aneuploidy_score"})
    df = df.merge(smp[["patient_id", "tmb", "fraction_genome_altered", "aneuploidy_score"]], on="patient_id", how="left")

    # Patient x gene mutation flags (protein-altering only)
    kept = mut[mut.mutationType.isin(PROTEIN_ALTERING)]
    flags = pd.crosstab(kept.patientId, kept.gene).clip(upper=1)
    for gene in DRIVER_GENES:
        df[f"mut_{gene}"] = df.patient_id.map(flags[gene] if gene in flags else {}).fillna(0).astype(int)
    # Clinically actionable alterations
    g12c = set(kept[(kept.gene == "KRAS") & (kept.proteinChange == "G12C")].patientId)
    egfr = kept[kept.gene == "EGFR"]
    is_act = pd.Series([egfr_classic_activating(p, t) for p, t in zip(egfr.proteinChange, egfr.mutationType)],
                       index=egfr.index, dtype=bool)  # a boolean Series, so an empty selection still works
    egfr_act = set(egfr.loc[is_act, "patientId"])
    df["kras_g12c"] = df.patient_id.isin(g12c).astype(int)
    df["egfr_activating"] = df.patient_id.isin(egfr_act).astype(int)

    # Copy-number flags: 2 = amplification, -2 = deep deletion
    df["amp_MET"] = df.patient_id.isin(cna[(cna.gene == "MET") & (cna.alteration == 2)].patientId).astype(int)
    df["amp_ERBB2"] = df.patient_id.isin(cna[(cna.gene == "ERBB2") & (cna.alteration == 2)].patientId).astype(int)
    df["del_CDKN2A"] = df.patient_id.isin(cna[(cna.gene == "CDKN2A") & (cna.alteration == -2)].patientId).astype(int)

    quality = {
        "patients": len(pat),
        "mutation rows downloaded": len(mut),
        "mutation rows excluded (not protein-altering)": int((~mut.mutationType.isin(PROTEIN_ALTERING)).sum()),
        "missing or unknown stage": int(df.stage.isna().sum()),
        "missing overall survival time or status": int((df.os_months.isna() | df.os_event.isna()).sum()),
        "overall survival time <= 0": int((df.os_months <= 0).sum()),
        "missing tumor mutational burden": int(df.tmb.isna().sum()),
        "missing genetic ancestry label": int(df.genetic_ancestry.isna().sum()),
    }
    return df, quality


def summaries(df):
    n = len(df)
    pct = lambda s: f"{100 * s.mean():.1f}%"
    t1 = [("Patients", f"{n}"), ("Age, median (IQR)", f"{df.age.median():.0f} ({df.age.quantile(.25):.0f}-"
                                                      f"{df.age.quantile(.75):.0f})"),
          ("Female", pct(df.sex == "Female"))]
    t1 += [(f"Stage {s}", pct(df.stage == s)) for s in ["I", "II", "III", "IV"]] + [("Stage unknown", pct(df.stage.isna()))]
    t1 += [("Deaths (overall survival events)", f"{int(df.os_event.sum())} ({pct(df.os_event == 1)})"),
           ("Median follow-up, months", f"{df.os_months.median():.1f}"),
           ("Tumor mutational burden, median (IQR)", f"{df.tmb.median():.1f} ({df.tmb.quantile(.25):.1f}-"
                                                     f"{df.tmb.quantile(.75):.1f})")]
    alt_cols = [c for c in df.columns if c.startswith(("mut_", "amp_", "del_"))] + ["kras_g12c", "egfr_activating"]
    freq = (df[alt_cols].mean().mul(100).round(1).sort_values(ascending=False)
            .rename("percent_of_patients").rename_axis("alteration").reset_index())
    ancestry = (df.genetic_ancestry.fillna("missing").value_counts().rename("patients").rename_axis("genetic_ancestry")
                .reset_index())
    ancestry["percent"] = (100 * ancestry.patients / n).round(1)
    return pd.DataFrame(t1, columns=["characteristic", "value"]), freq, ancestry


def main():
    df, quality = build()
    PROCESSED.mkdir(parents=True, exist_ok=True)
    RESULTS.mkdir(exist_ok=True)
    df.to_csv(PROCESSED / "patients.csv", index=False)
    t1, freq, ancestry = summaries(df)
    t1.to_csv(RESULTS / "table1.csv", index=False)
    freq.to_csv(RESULTS / "alteration_frequencies.csv", index=False)
    ancestry.to_csv(RESULTS / "genetic_ancestry.csv", index=False)
    pd.Series(quality).rename("count").to_csv(RESULTS / "data_quality.csv")
    print("Data quality:\n" + pd.Series(quality).to_string())
    print("\nTable 1:\n" + t1.to_string(index=False))
    print("\nAlteration frequencies (% of patients):\n" + freq.to_string(index=False))
    print("\nGenetic ancestry (representation check):\n" + ancestry.to_string(index=False))


if __name__ == "__main__":
    main()
