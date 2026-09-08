import numpy as np
import pandas as pd
import pytest
import statsmodels.api as sm

from matcheddoe.cli import example_config, example_path
from matcheddoe.common import read_table
from matcheddoe.core import analyse, design_matrix, fit_surface


def test_levo_anova_matches_independent_deposited_arithmetic():
    result = analyse(read_table(example_path()), **example_config("levo_design.csv"))
    summary = result.tables["model_summary"].iloc[0]
    assert summary.residual_ss == pytest.approx(38.8985233333333, abs=1e-9)
    assert summary.pure_error_ss == pytest.approx(21.1812666666667, abs=1e-9)
    assert summary.lack_of_fit_f == pytest.approx(0.2788195964, abs=1e-9)
    assert summary.residual_df == 8


def test_banana_two_responses_and_no_invented_replication():
    result = analyse(read_table(example_path("banana_design.csv")), **example_config("banana_design.csv"))
    summary = result.tables["model_summary"]
    np.testing.assert_allclose(summary.residual_ss, [18.6533333333, 1.1444444444], atol=1e-8)
    assert summary.pure_error_ss.isna().all() and summary.lack_of_fit_p.isna().all()
    assert len(result.tables["candidates"].filter(like="prediction::").columns) == 2


def test_coefficients_and_variance_match_independent_statsmodels_ols():
    rng = np.random.default_rng(90)
    coded = rng.uniform(-1, 1, size=(40, 3))
    x, _ = design_matrix(coded, ["A", "B", "C"])
    y = 5 + 2 * coded[:, 0] - coded[:, 1] ** 2 + rng.normal(size=40)
    ours = fit_surface(x, y)
    reference = sm.OLS(y, x).fit()
    np.testing.assert_allclose(ours[0], reference.params, atol=1e-11)
    assert ours[5] == pytest.approx(reference.mse_resid, abs=1e-11)
    np.testing.assert_allclose(ours[6], reference.normalized_cov_params, atol=1e-11)
    np.testing.assert_allclose(ours[7], reference.get_influence().hat_matrix_diag, atol=1e-11)


def test_missing_cells_are_excluded_and_rank_deficiency_gives_audit():
    data = read_table(example_path("banana_design.csv"))
    data.loc[0, "moisture_percent"] = np.nan
    result = analyse(data, **example_config("banana_design.csv"))
    assert result.tables["model_summary"].n.eq(8).all()
    sparse = pd.DataFrame({"sample_id": range(5), "A": range(5), "B": range(5), "response": range(5)})
    assert analyse(sparse, factors=["A", "B"]).status == "audit_only"


def test_box_behnken_candidates_do_not_extrapolate_to_unsupported_corners():
    result = analyse(read_table(example_path()), **example_config("levo_design.csv"))
    candidates = result.tables["candidates"]
    coded = (
        candidates[
            [result.settings["candidate_columns"]["factors"][f] for f in ["dose_g_L", "concentration_ppm", "pH"]]
        ].to_numpy()
        - [1, 7, 7]
    ) / [0.5, 3, 3]
    assert (np.abs(coded).sum(axis=1) <= 2 + 1e-8).all()
    assert not result.tables["selected_candidates"].measured_confirmation.any()


def test_unconfirmed_run_independence_preserves_fit_but_withholds_inference():
    config = example_config("levo_design.csv") | {"independent_runs": False}
    result = analyse(read_table(example_path()), **config)
    assert result.status == "partial"
    assert result.tables["coefficients"].coefficient_coded.notna().all()
    assert result.tables["coefficients"].ci95_low.isna().all()
    assert result.tables["model_summary"].pure_error_ss.isna().all()
    assert result.tables["selected_candidates"].prediction_interval95_low.isna().all()
