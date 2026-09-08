import numpy as np
import pandas as pd
import pytest

from matcheddoe import core
from matcheddoe.common import BlockPLS, InputError, metrics


@pytest.mark.parametrize("scale", [1e-12, 1e-9, 1, 1e6])
def test_metric_units_are_invariant(scale):
    a = metrics([1, 2, 3], [1.1, 1.8, 3.1], ["a", "b", "c"])
    b = metrics(np.array([1, 2, 3]) * scale, np.array([1.1, 1.8, 3.1]) * scale, ["a", "b", "c"])
    assert b["r2"] == pytest.approx(0.97, abs=1e-12)
    for key in ["rmse", "mae", "bias", "group_rmse"]:
        assert b[key] / scale == pytest.approx(a[key], abs=1e-12)


@pytest.mark.parametrize("scale", [1e-12, 1e-9, 1e6])
def test_pls_response_and_predictor_units_roundtrip(scale):
    x = np.array([[1.0, 2], [2, 1], [3, 7], [4, 3], [6, 8], [8, 10]])
    y = np.array([1.0, 2, 3, 2, 7, 8])
    a = BlockPLS(2).fit([x], y).predict([x])
    b = BlockPLS(2).fit([x * scale], y * scale).predict([x * scale]) / scale
    np.testing.assert_allclose(a, b, rtol=1e-8, atol=1e-10)


def test_metrics_translation_negative_r2_and_missing_pairs():
    y = np.array([1.0, 2, 3])
    p = np.array([5.0, 6, 7])
    assert metrics(y, p)["r2"] < 0
    assert metrics(y + 100, p + 100)["r2"] == metrics(y, p)["r2"]
    assert metrics([0, 0], [0, 0])["r2"] is None
    assert metrics([np.nan, 1], [1, np.nan])["n"] == 0


def test_api_rejects_unknown_configuration_with_field_name():
    with pytest.raises(InputError, match="unknown_option"):
        core.analyse(pd.DataFrame(), unknown_option=True)


from matcheddoe.cli import example_config, example_path
from matcheddoe.common import read_table


@pytest.mark.parametrize("value", ["false", "true", None, 0, 1])
def test_independence_requires_real_boolean(value):
    with pytest.raises(InputError, match="independent_runs"):
        core.analyse(read_table(example_path()), **(example_config("levo_design.csv") | {"independent_runs": value}))


def test_duplicate_responses_rejected_before_fit():
    with pytest.raises(InputError, match="responses"):
        core.analyse(
            read_table(example_path()), **(example_config("levo_design.csv") | {"responses": ["response", "response"]})
        )


@pytest.mark.parametrize("scale", [1e-12, 1e-9, 1e6])
def test_doe_response_units_preserve_dimensionless_inference(scale):
    data = read_table(example_path())
    config = example_config("levo_design.csv")
    a = core.analyse(data, **config)
    data["response"] = pd.to_numeric(data.response) * scale
    b = core.analyse(data, **config)
    for key in ["r2_fitted", "adjusted_r2", "lack_of_fit_f", "lack_of_fit_p"]:
        assert b.tables["model_summary"].iloc[0][key] == pytest.approx(a.tables["model_summary"].iloc[0][key], rel=1e-8)
    for key in ["coefficient_coded", "standard_error", "ci95_low", "ci95_high"]:
        np.testing.assert_allclose(
            b.tables["coefficients"][key] / scale, a.tables["coefficients"][key], rtol=1e-7, atol=1e-9
        )
