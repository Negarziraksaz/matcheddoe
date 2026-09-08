"""Closeout regressions: report selection, cohort traces and export boundaries."""

import io
import json
from unittest.mock import patch

import numpy as np
import pandas as pd
import pytest

from matcheddoe import core, templates
from matcheddoe.presentation import select_view
from matcheddoe.reporting import files_for_result, summary_files


def test_summary_real_paragraphs_and_traceable_view():
    result = core.analyse(templates.starter_frame(), **templates.schema()["config"])
    with patch.object(core, "analyse", side_effect=AssertionError("Viewing must not refit")):
        output = summary_files(result)
        text = output["summary.txt"].decode("utf-8")
        assert "\n\n" in text and r"\n" not in text
        view = json.loads(output["summary_view.json"])
        cohort = pd.read_csv(io.BytesIO(output["summary_cohort.csv"]))
        assert cohort.included_in_primary_score_and_figure.sum() == view["n_primary_pairs"]


def test_stale_selection_resolves_and_zip_summary_matches():
    result = core.analyse(templates.starter_frame(), **templates.schema()["config"])
    default = select_view(result)
    stale = select_view(result, method="removed", domain="removed", response="removed", budget=9999)
    for key in ["method", "domain", "response", "budget", "scope"]:
        assert stale[key] == default[key]
    selection = {k: default[k] for k in ["method", "domain", "response", "budget"]}
    full = files_for_result(result, formats=("png",), png_dpi=70, **selection)
    assert json.loads(full["summary_view.json"])["scope"] == default["scope"]
    assert full["summary.txt"] == summary_files(result)["summary.txt"]


@pytest.mark.parametrize(
    "name",
    [
        "response",
        "direction",
        "predicted_response",
        "group",
        "model",
        "method",
        "residual",
        "prediction",
        "budget",
        "prediction_interval95_low",
        "inside_observed_convex_hull",
        "measured_confirmation",
        "Response__prediction",
        "Temperature (C)",
    ],
)
def test_candidate_factor_rename_preserves_numbers_and_metadata(name):
    frame = templates.starter_frame()
    before = core.analyse(frame, **templates.schema()["config"])
    frame = frame.rename(columns={"Time (min)": name})
    after = core.analyse(frame, **(templates.schema()["config"] | {"factors": [name]}))
    a = before.tables["selected_candidates"].iloc[0]
    b = after.tables["selected_candidates"].iloc[0]
    for key in ["response", "direction", "inside_observed_convex_hull", "measured_confirmation"]:
        assert a[key] == b[key]
    for key in ["predicted_response", "prediction_interval95_low", "prediction_interval95_high"]:
        assert (pd.isna(a[key]) and pd.isna(b[key])) or a[key] == pytest.approx(b[key])
    factor_key = after.settings["candidate_columns"]["factors"][name]
    assert b[factor_key] == pytest.approx(a[before.settings["candidate_columns"]["factors"]["Time (min)"]])
    np.testing.assert_allclose(before.tables["predictions"].prediction, after.tables["predictions"].prediction)
    assert len(after.tables["candidates"].columns) == len(set(after.tables["candidates"].columns))
    if name in {"response", "direction", "predicted_response"}:
        output = files_for_result(after, formats=("png",), png_dpi=70)
        csv = pd.read_csv(io.BytesIO(output["selected_candidates.csv"]))
        assert csv.iloc[0][factor_key] == pytest.approx(b[factor_key])
        book = pd.ExcelFile(io.BytesIO(output["results.xlsx"]))
        xls = pd.read_excel(book, sheet_name=next(n for n in book.sheet_names if n.endswith("selected_candidates")))
        assert xls.iloc[0][factor_key] == pytest.approx(b[factor_key])
        assert not after.settings.get("diagnostic_failures")


def test_duplicate_headers_rejected_before_fitting():
    frame = templates.starter_frame()
    frame.columns = ["duplicate"] * len(frame.columns)
    with pytest.raises(ValueError, match="distinct"):
        core.analyse(frame, **templates.schema()["config"])


def test_per_response_candidate_metadata_survives_response_factor():
    frame = templates.starter_frame().rename(columns={"Time (min)": "response"})
    frame["Other response"] = frame["Response"] * 2
    result = core.analyse(
        frame,
        **(
            templates.schema()["config"]
            | {"factors": ["response"], "responses": ["Response", "Other response"], "response_mode": "per_response"}
        ),
    )
    assert set(result.tables["candidates_by_response"].response) == {"Response", "Other response"}
    assert "factor::response" in result.tables["candidates_by_response"]
    assert not result.settings.get("diagnostic_failures")
