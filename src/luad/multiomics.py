"""Multi-omics: does gene expression add prognostic information beyond clinical factors and mutations?

Run:  python -m luad.multiomics      (after `python -m luad.fetch` has downloaded data/raw/expression_panel.csv)
Steps:
  1. Expression signature scores: log2(expression + 1), standardized per gene, averaged within each signature.
  2. Biological consistency checks across omics layers (known biology the data should reproduce):
       - KEAP1-mutant tumors should have higher NRF2 target expression (KEAP1 normally degrades NRF2)
       - STK11-mutant tumors are reported to be less immune-infiltrated
  3. Cox models on the SAME patients: clinical -> + mutations -> + expression signatures, compared with
     cross-validated C-index.
  4. A Cox model stratified by stage, which removes the proportional-hazards assumption for stage.
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

from luad.fetch import SIGNATURES
from luad.survival import GENOMIC_COVARIATES, cox, cv_cindex, survival_cohort

RAW, PROCESSED, RESULTS = Path("data/raw"), Path("data/processed"), Path("results")
SCORES = [f"score_{name}" for name in SIGNATURES]


def signature_scores(expr):
    """expr: long table (patientId, gene, value) -> one row per patient with one score per signature."""
    wide = expr.pivot_table(index="patientId", columns="gene", values="value", aggfunc="mean")
    z = np.log2(wide.clip(lower=0) + 1)
    z = (z - z.mean()) / z.std(ddof=0)  # unsupervised per-gene standardization across the cohort
    scores = pd.DataFrame(index=z.index)
    for name, genes in SIGNATURES.items():
        present = [g for g in genes if g in z.columns]
        scores[f"score_{name}"] = z[present].mean(axis=1)
    return scores.rename_axis("patient_id").reset_index()


def consistency_checks(d):
    checks = []
    for mut, score, expected in [("mut_KEAP1", "score_nrf2_targets", "higher"),
                                 ("mut_STK11", "score_cytotoxic_immune", "lower"),
                                 ("mut_STK11", "score_immune_checkpoint", "lower")]:
        a, b = d.loc[d[mut] == 1, score].dropna(), d.loc[d[mut] == 0, score].dropna()
        diff = a.median() - b.median()
        checks.append({"comparison": f"{score.replace('score_', '')} in {mut.replace('mut_', '')}-mutant vs wild type",
                       "expected": expected, "median_difference": round(diff, 2),
                       "direction_matches": (diff > 0) == (expected == "higher"),
                       "mann_whitney_p": mannwhitneyu(a, b).pvalue})
    return pd.DataFrame(checks)


def plot_checks(d, path):
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.5))
    for ax, (mut, score, title) in zip(axes, [("mut_KEAP1", "score_nrf2_targets", "NRF2 target expression"),
                                               ("mut_STK11", "score_cytotoxic_immune", "Cytotoxic immune expression")]):
        groups = [d.loc[d[mut] == v, score].dropna() for v in (0, 1)]
        ax.boxplot(groups, tick_labels=["wild type", f"{mut.replace('mut_', '')} mutant"])
        ax.set_title(title)
        ax.set_ylabel("Signature score (z)")
    fig.suptitle("Cross-omics consistency: mutations and the expression programs they should affect")
    fig.tight_layout()
    fig.savefig(path, dpi=150)


def main():
    patients = pd.read_csv(PROCESSED / "patients.csv")
    expr = pd.read_csv(RAW / "expression_panel.csv")
    scores = signature_scores(expr)
    d = patients.merge(scores, on="patient_id", how="left")
    print(f"Expression available for {d[SCORES[0]].notna().sum()} of {len(d)} patients")

    checks = consistency_checks(d)
    checks["mann_whitney_p"] = checks.mann_whitney_p.map(lambda p: f"{p:.2e}")
    checks.to_csv(RESULTS / "multiomics_consistency.csv", index=False)
    plot_checks(d, RESULTS / "multiomics_consistency.png")
    print("\nCross-omics consistency checks:\n" + checks.to_string(index=False))

    s, _ = survival_cohort(d)
    s = s.dropna(subset=["age"] + SCORES)
    s = s[s.stage != "unknown"]  # too few for a separate category (see survival.py)
    print(f"\nSurvival models on the same {len(s)} patients ({int(s.os_event.sum())} deaths) with complete "
          "clinical, mutation, and expression data")

    clinical = ["age", "female", "stage_II", "stage_III", "stage_IV"]
    models = {"clinical": clinical, "+ mutations": clinical + GENOMIC_COVARIATES,
              "+ mutations + expression": clinical + GENOMIC_COVARIATES + SCORES}
    rows = []
    for name, cov in models.items():
        mean, sd = cv_cindex(s, cov)
        rows.append({"model": name, "covariates": len(cov), "cv_c_index": round(mean, 3), "cv_sd": round(sd, 3),
                     "apparent_c_index": round(cox(s, cov).concordance_index_, 3)})
    comp = pd.DataFrame(rows)
    comp.to_csv(RESULTS / "multiomics_c_index.csv", index=False)
    print("\nModel comparison (5 x 5-fold cross-validated C-index):\n" + comp.to_string(index=False))

    full = cox(s, models["+ mutations + expression"])
    hr = full.summary.loc[GENOMIC_COVARIATES + SCORES, ["exp(coef)", "exp(coef) lower 95%", "exp(coef) upper 95%", "p"]]
    hr.columns = ["hazard_ratio", "hr_lower_95", "hr_upper_95", "p"]
    hr.round(3).to_csv(RESULTS / "cox_multiomics.csv")
    print("\nFull model, molecular covariates (hazard ratios per alteration or per 1 SD of score):\n"
          + hr.round(3).to_string())

    # Stratified by stage: each stage gets its own baseline hazard, so stage no longer needs proportional hazards
    from lifelines import CoxPHFitter
    strat_cov = ["age", "female"] + GENOMIC_COVARIATES + SCORES
    strat = CoxPHFitter(penalizer=0.01).fit(s[["os_months", "os_event", "stage"] + strat_cov],
                                            duration_col="os_months", event_col="os_event", strata=["stage"])
    st = strat.summary.loc[["mut_STK11"] + SCORES, ["exp(coef)", "exp(coef) lower 95%", "exp(coef) upper 95%", "p"]]
    st.columns = ["hazard_ratio", "hr_lower_95", "hr_upper_95", "p"]
    st.round(3).to_csv(RESULTS / "cox_stratified_by_stage.csv")
    print("\nStratified by stage (sensitivity analysis):\n" + st.round(3).to_string())


if __name__ == "__main__":
    main()
