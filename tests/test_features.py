import pandas as pd

from luad.features import build, egfr_classic_activating, event_flag, stage_group


def test_stage_grouping():
    assert stage_group("STAGE IIIA") == "III" and stage_group("Stage IB") == "I" and stage_group("STAGE IV") == "IV"
    assert stage_group(None) is None and stage_group("[Not Available]") is None


def test_survival_status_parsing():
    assert event_flag("1:DECEASED") == 1 and event_flag("0:LIVING") == 0 and event_flag(None) is None


def test_egfr_activating_detection():
    assert egfr_classic_activating("L858R", "Missense_Mutation")
    assert egfr_classic_activating("E746_A750del", "In_Frame_Del")
    assert not egfr_classic_activating("T790M", "Missense_Mutation")
    assert not egfr_classic_activating("E746_A750del", "Frame_Shift_Del")


def test_build_flags(tmp_path):
    a, b = "TCGA-AA-0001", "TCGA-AA-0002"          # TCGA patient barcodes are 12 characters
    pd.DataFrame({"patientId": [a, b], "AGE": [60, 70], "SEX": ["Female", "Male"],
                  "AJCC_PATHOLOGIC_TUMOR_STAGE": ["STAGE IB", None], "GENETIC_ANCESTRY_LABEL": ["EUR", None],
                  "OS_MONTHS": [10.0, 5.0], "OS_STATUS": ["1:DECEASED", "0:LIVING"],
                  "PFS_MONTHS": [8.0, 5.0], "PFS_STATUS": ["1:Recurred/Progressed", "0:CENSORED"]}
                 ).to_csv(tmp_path / "clinical_patient.csv", index=False)
    pd.DataFrame({"sampleId": [a + "-01", b + "-01"], "TMB_NONSYNONYMOUS": [5.0, 2.0],
                  "FRACTION_GENOME_ALTERED": [.2, .1], "ANEUPLOIDY_SCORE": [10, 3]}
                 ).to_csv(tmp_path / "clinical_sample.csv", index=False)
    pd.DataFrame({"patientId": [a, a, b], "sampleId": [a + "-01", a + "-01", b + "-01"],
                  "gene": ["KRAS", "TP53", "TP53"], "mutationType": ["Missense_Mutation", "Silent", "Nonsense_Mutation"],
                  "proteinChange": ["G12C", "P72P", "R213*"], "variantType": ["SNP"] * 3}
                 ).to_csv(tmp_path / "mutations_drivers.csv", index=False)
    pd.DataFrame({"patientId": [b], "sampleId": [b + "-01"], "gene": ["MET"], "alteration": [2]}
                 ).to_csv(tmp_path / "cna_drivers.csv", index=False)
    df, quality = build(tmp_path)
    row = df.set_index("patient_id")
    assert row.loc[a, "mut_KRAS"] == 1 and row.loc[a, "kras_g12c"] == 1
    assert row.loc[a, "mut_TP53"] == 0          # silent mutation excluded
    assert row.loc[b, "mut_TP53"] == 1 and row.loc[b, "amp_MET"] == 1
    assert row.loc[a, "stage"] == "I" and pd.isna(row.loc[b, "stage"])
    assert row.loc[a, "tmb"] == 5.0
    assert quality["mutation rows excluded (not protein-altering)"] == 1
