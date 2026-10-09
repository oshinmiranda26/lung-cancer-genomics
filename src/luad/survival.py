"""Survival analysis: do tumor alterations add prognostic information to clinical factors?

Run:  python -m luad.survival
Steps:
  1. Cohort: patients with usable overall survival (time > 0 and a known status).
  2. Kaplan-Meier curves and log-rank tests for each key alteration, with false discovery rate correction.
  3. Co-mutation question: within KRAS-mutant tumors, do STK11 or KEAP1 co-mutations change survival?
  4. Cox models: clinical only (age, sex, stage) vs clinical + genomic, compared with cross-validated C-index,
     plus a proportional-hazards check.
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import logrank_test, proportional_hazard_test
from lifelines.utils import k_fold_cross_validation

PROCESSED, RESULTS = Path("data/processed"), Path("results")
ALTERATIONS = {"mut_TP53": "TP53 mutation", "mut_KRAS": "KRAS mutation", "mut_STK11": "STK11 mutation",
               "mut_KEAP1": "KEAP1 mutation", "mut_EGFR": "EGFR mutation", "del_CDKN2A": "CDKN2A deep deletion",
               "mut_SMARCA4": "SMARCA4 mutation", "kras_g12c": "KRAS G12C"}
GENOMIC_COVARIATES = ["mut_TP53", "mut_KRAS", "mut_STK11", "mut_KEAP1", "mut_EGFR", "del_CDKN2A", "log_tmb"]


def benjamini_hochberg(p):
    """False discovery rate adjusted p-values (q-values)."""
    p = np.asarray(p, float)
    order = np.argsort(p)
    ranked = p[order] * len(p) / np.arange(1, len(p) + 1)
    q = np.minimum.accumulate(ranked[::-1])[::-1].clip(max=1)
    out = np.empty_like(q)
    out[order] = q
    return out


def survival_cohort(df):
    keep = df.os_months.notna() & df.os_event.notna() & (df.os_months > 0)
    d = df[keep].copy()
    d["os_event"] = d.os_event.astype(int)
    d["log_tmb"] = np.log1p(d.tmb)
    d["female"] = (d.sex == "Female").astype(int)
    # Stage as indicators with stage I as reference; unknown stage kept as its own explicit category
    d["stage"] = d.stage.fillna("unknown")
    for s in ["II", "III", "IV", "unknown"]:
        d[f"stage_{s}"] = (d.stage == s).astype(int)
    return d, int((~keep).sum())


def km_table(d):
    rows = []
    for col, label in ALTERATIONS.items():
        g1, g0 = d[d[col] == 1], d[d[col] == 0]
        if len(g1) < 10:
            continue
        test = logrank_test(g1.os_months, g0.os_months, g1.os_event, g0.os_event)
        med = lambda g: KaplanMeierFitter().fit(g.os_months, g.os_event).median_survival_time_
        rows.append({"alteration": label, "n_altered": len(g1), "deaths_altered": int(g1.os_event.sum()),
                     "median_os_altered": round(med(g1), 1), "median_os_wildtype": round(med(g0), 1),
                     "logrank_p": test.p_value})
    t = pd.DataFrame(rows)
    t["fdr_q"] = benjamini_hochberg(t.logrank_p)
    return t.round({"logrank_p": 4, "fdr_q": 4}).sort_values("logrank_p")


def plot_km(d, path):
    cols = list(ALTERATIONS)[:6]
    fig, axes = plt.subplots(2, 3, figsize=(13, 8), sharey=True)
    for ax, col in zip(axes.flat, cols):
        for val, color, name in [(0, "#5e6e71", "wild type"), (1, "#c0392b", "altered")]:
            g = d[d[col] == val]
            KaplanMeierFitter().fit(g.os_months, g.os_event, label=f"{name} (n={len(g)})").plot_survival_function(
                ax=ax, color=color, ci_show=False)
        ax.set_title(ALTERATIONS[col])
        ax.set_xlabel("Months")
        ax.set_xlim(0, 120)
    axes[0, 0].set_ylabel("Overall survival")
    axes[1, 0].set_ylabel("Overall survival")
    fig.suptitle("Overall survival by tumor alteration (TCGA lung adenocarcinoma)")
    fig.tight_layout()
    fig.savefig(path, dpi=150)


def kras_comutation(d, path):
    """Within KRAS-mutant tumors: STK11 or KEAP1 co-mutation vs neither."""
    k = d[d.mut_KRAS == 1].copy()
    k["co"] = ((k.mut_STK11 == 1) | (k.mut_KEAP1 == 1)).astype(int)
    a, b = k[k.co == 1], k[k.co == 0]
    test = logrank_test(a.os_months, b.os_months, a.os_event, b.os_event)
    fig, ax = plt.subplots(figsize=(6, 5))
    for g, color, name in [(b, "#5e6e71", "KRAS only"), (a, "#c0392b", "KRAS + STK11/KEAP1")]:
        KaplanMeierFitter().fit(g.os_months, g.os_event, label=f"{name} (n={len(g)})").plot_survival_function(
            ax=ax, color=color, ci_show=True)
    ax.set_xlim(0, 120)
    ax.set_xlabel("Months")
    ax.set_ylabel("Overall survival")
    ax.set_title(f"KRAS-mutant tumors: co-mutation (log-rank p = {test.p_value:.3f})")
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    return {"kras_mutant": len(k), "co_mutated": len(a), "logrank_p": round(test.p_value, 4)}


def cox(d, covariates, penalizer=0.01):
    cph = CoxPHFitter(penalizer=penalizer)
    cph.fit(d[["os_months", "os_event"] + covariates], duration_col="os_months", event_col="os_event")
    return cph


def cv_cindex(d, covariates, folds=5, repeats=5):
    scores = []
    for r in range(repeats):
        shuffled = d.sample(frac=1, random_state=r)
        scores += k_fold_cross_validation(CoxPHFitter(penalizer=0.01), shuffled[["os_months", "os_event"] + covariates],
                                          duration_col="os_months", event_col="os_event", k=folds,
                                          scoring_method="concordance_index", seed=r)
    return np.mean(scores), np.std(scores)


def main():
    df = pd.read_csv(PROCESSED / "patients.csv")
    d, excluded = survival_cohort(df)
    print(f"Survival cohort: {len(d)} patients, {int(d.os_event.sum())} deaths "
          f"({excluded} excluded for missing or non-positive survival time)")

    km = km_table(d)
    km.to_csv(RESULTS / "km_by_alteration.csv", index=False)
    plot_km(d, RESULTS / "km_by_alteration.png")
    print("\nKaplan-Meier by alteration (median overall survival in months):\n" + km.to_string(index=False))

    co = kras_comutation(d, RESULTS / "kras_comutation.png")
    print(f"\nKRAS-mutant tumors: {co['kras_mutant']}, with STK11/KEAP1 co-mutation: {co['co_mutated']}, "
          f"log-rank p = {co['logrank_p']}")

    clinical = ["age", "female", "stage_II", "stage_III", "stage_IV", "stage_unknown"]
    d = d.dropna(subset=["age"])
    models = {"clinical": clinical, "clinical + genomic": clinical + GENOMIC_COVARIATES}
    comparison = []
    for name, cov in models.items():
        cph = cox(d, cov)
        hr = cph.summary[["exp(coef)", "exp(coef) lower 95%", "exp(coef) upper 95%", "p"]].round(3)
        hr.columns = ["hazard_ratio", "hr_lower_95", "hr_upper_95", "p"]
        hr.to_csv(RESULTS / f"cox_{name.replace(' + ', '_plus_').replace(' ', '_')}.csv")
        mean, sd = cv_cindex(d, cov)
        comparison.append({"model": name, "covariates": len(cov), "cv_c_index": round(mean, 3), "cv_sd": round(sd, 3),
                           "apparent_c_index": round(cph.concordance_index_, 3)})
        print(f"\nCox model, {name} (hazard ratios):\n" + hr.to_string())
    comp = pd.DataFrame(comparison)
    comp.to_csv(RESULTS / "c_index_comparison.csv", index=False)
    print("\nModel comparison (5 x 5-fold cross-validated C-index):\n" + comp.to_string(index=False))

    ph = proportional_hazard_test(cox(d, models["clinical + genomic"]), d[["os_months", "os_event"] +
                                  models["clinical + genomic"]], time_transform="rank").summary[["p"]].round(4)
    ph.to_csv(RESULTS / "proportional_hazards_test.csv")
    violations = ph[ph.p < 0.05].index.tolist()
    print("\nProportional-hazards check (Schoenfeld residuals): " +
          (f"possible violations for {violations}" if violations else "no covariate with p < 0.05"))


if __name__ == "__main__":
    main()
