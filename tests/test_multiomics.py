import numpy as np
import pandas as pd

from luad.fetch import SIGNATURES
from luad.multiomics import consistency_checks, signature_scores


def test_signature_scores_are_standardized_averages():
    rng = np.random.default_rng(0)
    rows = [{"patientId": f"P{i}", "gene": g, "value": float(rng.exponential(100))}
            for i in range(50) for genes in SIGNATURES.values() for g in genes]
    scores = signature_scores(pd.DataFrame(rows))
    assert len(scores) == 50 and {f"score_{n}" for n in SIGNATURES} <= set(scores.columns)
    assert abs(scores.score_proliferation.mean()) < 1e-9          # average of standardized genes


def test_consistency_check_detects_a_planted_effect():
    rng = np.random.default_rng(1)
    n = 400
    d = pd.DataFrame({"mut_KEAP1": rng.integers(0, 2, n), "mut_STK11": rng.integers(0, 2, n)})
    d["score_nrf2_targets"] = rng.normal(size=n) + 1.5 * d.mut_KEAP1           # planted: higher in KEAP1-mutant
    d["score_cytotoxic_immune"] = rng.normal(size=n) - 1.0 * d.mut_STK11       # planted: lower in STK11-mutant
    d["score_immune_checkpoint"] = rng.normal(size=n)
    c = consistency_checks(d).set_index("comparison")
    assert c.iloc[0].direction_matches and c.iloc[0].mann_whitney_p < 1e-6
    assert c.iloc[1].direction_matches
