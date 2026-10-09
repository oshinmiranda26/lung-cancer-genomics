import numpy as np
import pandas as pd

from luad.survival import benjamini_hochberg, cox, survival_cohort


def test_benjamini_hochberg_matches_known_values():
    q = benjamini_hochberg([0.01, 0.04, 0.03, 0.20])
    assert np.allclose(q, [0.04, 0.0533333, 0.0533333, 0.20], atol=1e-6)


def test_survival_cohort_excludes_unusable_rows_and_codes_stage():
    df = pd.DataFrame({"os_months": [10, 0, None, 5], "os_event": [1, 0, 1, None], "tmb": [1, 2, 3, 4],
                       "sex": ["Female", "Male", "Male", "Female"], "stage": ["I", "II", None, None]})
    d, excluded = survival_cohort(df)
    assert len(d) == 1 and excluded == 3 and d.stage_unknown.sum() == 0


def test_cox_recovers_a_planted_hazard_ratio():
    rng = np.random.default_rng(0)
    n = 3000
    x = rng.integers(0, 2, n)
    t = rng.exponential(1 / (0.05 * np.exp(np.log(2.0) * x)))   # true hazard ratio 2 for x = 1
    c = rng.exponential(40, n)                                    # random censoring
    d = pd.DataFrame({"os_months": np.minimum(t, c), "os_event": (t <= c).astype(int), "x": x})
    hr = cox(d, ["x"], penalizer=0.0).hazard_ratios_["x"]
    assert 1.8 < hr < 2.2
