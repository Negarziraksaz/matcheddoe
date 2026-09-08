import io
import json

import numpy as np
import pandas as pd
import pytest

from matcheddoe import core
from matcheddoe.cli import example_config, example_context, example_path
from matcheddoe.common import InputError, read_table
from matcheddoe.reporting import files_for_result


def test_resolved_configuration_and_portable_public_attribution():
    path = example_path()
    data = read_table(path)
    a = core.analyse(data, **(example_config(path.name) | example_context(path.name)))
    files = files_for_result(a)
    config = json.loads(files["config.json"])
    b = core.analyse(data, **config)
    assert a.settings["input_sha256"] == b.settings["input_sha256"]
    assert json.loads(files["run.json"])["software_version"] == "0.5.2"
    assert json.loads(files["provenance.json"])["original_source"]["license"] == "CC-BY-4.0"
    assert b.settings["provenance"]["kind"] == "bundled_example"
    assert "created_utc" not in config
    for key, table in a.tables.items():
        pd.testing.assert_frame_equal(table, b.tables[key], check_exact=False, rtol=1e-7, atol=1e-9)


def test_uploaded_data_does_not_inherit_bundled_citation():
    path = example_path()
    data = read_table(path)
    result = core.analyse(data, **example_config(path.name))
    assert result.settings["provenance"]["kind"] == "user_provided"
    assert "original_source" not in result.settings["provenance"]


def test_nonobject_cli_config_returns_clear_failure(tmp_path, capsys):
    from matcheddoe.cli import main

    file = tmp_path / "bad.json"
    file.write_text("[1,2]")
    assert main(["demo", "--config", str(file), "--output", str(tmp_path / "out")]) == 2
    assert "top level must be an object" in capsys.readouterr().err


def test_duplicate_xlsx_headers_are_rejected():
    frame = pd.DataFrame([["id", "x", "x"], ["001", 1, 2]])
    source = io.BytesIO()
    frame.to_excel(source, index=False, header=False)
    source.seek(0)
    with pytest.raises(InputError, match="unique"):
        read_table(source, "duplicate.xlsx")


def test_html_escapes_untrusted_labels():
    from matcheddoe.common import Result
    from matcheddoe.reporting import html_report

    result = Result("<script>alert(1)</script>", notes=["<img src=x onerror=alert(1)>"])
    html = html_report(result, "<bad>", False)
    assert "<script>alert" not in html and "&lt;script&gt;" in html


import statsmodels.api as sm


def known_surface():
    a, b = np.meshgrid([-1.0, 0, 1], [-1.0, 0, 1])
    a = a.ravel()
    b = b.ravel()
    y = 4 + 2 * a - 3 * b + 0.5 * a * b
    return pd.DataFrame({"sample_id": [f"r{i}" for i in range(9)], "A": a, "B": b, "response": y, "reported": y})


def equation():
    return {
        "response": {
            "basis": "actual",
            "coefficients": {"Intercept": 4.0, "A": 2.0, "B": -3.0, "A:B": 0.5},
            "prediction_column": "reported",
            "atol": 1e-10,
            "rtol": 1e-10,
        }
    }


def test_stated_equation_and_deliberate_discrepancy():
    frame = known_surface()
    specification = equation()
    correct = core.analyse(frame, factors=["A", "B"], model="2fi", reported_equations=specification)
    assert correct.tables["equation_comparison"].status.eq("within_tolerance").all()
    specification["response"]["coefficients"]["A:B"] = 0.8
    wrong = core.analyse(frame, factors=["A", "B"], model="2fi", reported_equations=specification)
    assert wrong.tables["equation_comparison"].status.eq("discrepancy").sum() == 8
    assert not wrong.tables["equation_comparison"].measured_confirmation.any()


def test_coded_equation_requires_explicit_coding_and_rounding_is_bounded():
    frame = known_surface()
    specification = equation()
    specification["response"]["basis"] = "coded"
    with pytest.raises(InputError, match="coding"):
        core.analyse(frame, factors=["A", "B"], model="2fi", reported_equations=specification)
    specification["response"]["coding"] = {
        "A": {"centre": 0.0, "half_range": 1.0},
        "B": {"centre": 0.0, "half_range": 1.0},
    }
    specification["response"]["coefficients"]["Intercept"] = 4.0001
    specification["response"]["coefficient_rounding"] = 0.0001
    r = core.analyse(frame, factors=["A", "B"], model="2fi", reported_equations=specification)
    assert r.tables["equation_comparison"].status.eq("within_tolerance").all()


def test_overall_anova_matches_independent_ols():
    frame = known_surface()
    frame["response"] += np.array([0.1, -0.2, 0.3, 0.2, 0.1, -0.1, -0.3, 0.2, 0.1])
    r = core.analyse(frame, factors=["A", "B"], model="2fi", independent_runs=True)
    ref = sm.OLS(frame.response, np.column_stack([np.ones(9), frame.A, frame.B, frame.A * frame.B])).fit()
    table = r.tables["anova"].set_index("source")
    assert table.loc["model", "f"] == pytest.approx(ref.fvalue, rel=1e-10)
    assert table.loc["model", "p"] == pytest.approx(ref.f_pvalue, rel=1e-10)
    assert table.loc["model", "ss"] == pytest.approx(ref.ess, rel=1e-10)
    assert table.loc["total_corrected", "ss"] == pytest.approx(table.ss.iloc[:2].sum(), rel=1e-12)


def test_each_response_exports_a_figure():
    path = example_path("banana_design.csv")
    r = core.analyse(read_table(path), **(example_config(path.name) | example_context(path.name)))
    out = files_for_result(r)
    assert {"figure.png", "figure_response_2.png", "figure_response_2.pdf"} <= set(out)
    assert b"colour_deltaE" in out["figure_response_2_caption.txt"]


def test_perfect_fit_and_constant_response_withhold_unresolved_variance():
    for data in [known_surface(), known_surface().assign(response=4.0)]:
        r = core.analyse(data, factors=["A", "B"], model="2fi", independent_runs=True)
        assert r.status == "partial"
        assert r.tables["coefficients"].ci95_low.isna().all()
        assert r.tables["anova"].f.isna().all()
