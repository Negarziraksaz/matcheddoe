"""Evaluate declared polynomial coefficients without parsing executable expressions."""

import numpy as np

from .common import InputError, numbers
from .contract import real_number


def compare_equation(frame, factors, model, specification, response, ids, fitted):
    """Compare a stated equation at observed rows with reported/fitted predictions.

    Coefficient rounding is an absolute bound for each stated coefficient, in
    that coefficient's declared basis. It is propagated as sum(abs(term))*bound.
    Discrepancies describe arithmetic, not a judgement about the original study.
    """
    from .core import design_matrix

    allowed = {"basis", "coding", "coefficients", "prediction_column", "atol", "rtol", "coefficient_rounding"}
    if not isinstance(specification, dict) or set(specification) - allowed:
        raise InputError(
            "reported_equations: use basis, coding, coefficients, prediction_column, atol, rtol and coefficient_rounding fields only."
        )
    basis = specification.get("basis")
    if basis not in ["actual", "coded"]:
        raise InputError("reported_equations.basis: explicitly select actual or coded factor values.")
    x = numbers(frame, factors).to_numpy()
    if basis == "coded":
        coding = specification.get("coding")
        if not isinstance(coding, dict) or set(coding) != set(factors):
            raise InputError("reported_equations.coding: provide centre and half_range for every factor.")
        for i, name in enumerate(factors):
            item = coding[name]
            if not isinstance(item, dict) or set(item) != {"centre", "half_range"}:
                raise InputError(f"reported_equations.coding.{name}: provide centre and half_range.")
            centre = item["centre"]
            width = item["half_range"]
            if isinstance(centre, bool) or not isinstance(centre, (int, float)) or not np.isfinite(centre):
                raise InputError("reported_equations.coding: centre must be finite.")
            real_number(width, "reported_equations.coding.half_range")
            if width == 0:
                raise InputError("reported_equations.coding.half_range must be positive.")
            x[:, i] = (x[:, i] - centre) / width
    elif specification.get("coding"):
        raise InputError("reported_equations.coding: omit coding when basis is actual.")
    matrix, terms = design_matrix(x, factors, model)
    coefficients = specification.get("coefficients")
    if not isinstance(coefficients, dict) or set(coefficients) != set(terms):
        raise InputError(
            f"reported_equations.coefficients: supply each term exactly once, including explicit zeroes: {terms}."
        )
    values = []
    for term in terms:
        value = coefficients[term]
        if isinstance(value, bool) or not isinstance(value, (int, float)) or not np.isfinite(value):
            raise InputError(f"reported_equations.coefficients.{term}: supply a finite number.")
        values.append(value)
    atol = real_number(specification.get("atol", 1e-8), "reported_equations.atol")
    rtol = real_number(specification.get("rtol", 1e-6), "reported_equations.rtol")
    rounding = real_number(specification.get("coefficient_rounding", 0.0), "reported_equations.coefficient_rounding")
    equation = matrix @ values
    allowance = np.sum(np.abs(matrix), axis=1) * rounding
    comparisons = {"refitted_prediction": np.asarray(fitted)}
    column = specification.get("prediction_column")
    if column is not None:
        if not isinstance(column, str) or column in factors or column == response:
            raise InputError(
                "reported_equations.prediction_column: choose a separate reported-prediction column; measured response is not a prediction."
            )
        comparisons["reported_prediction"] = numbers(frame, [column])[column].to_numpy()
    rows = []
    for origin, other in comparisons.items():
        for i in range(len(frame)):
            difference = equation[i] - other[i]
            tolerance = atol + rtol * abs(other[i]) + allowance[i]
            finite = np.isfinite(other[i]) and np.isfinite(equation[i])
            status = (
                ("within_tolerance" if abs(difference) <= tolerance else "discrepancy") if finite else "unavailable"
            )
            rows.append(
                {
                    "sample_id": ids[i],
                    "response": response,
                    "basis": basis,
                    "comparison_origin": origin,
                    "equation_prediction": equation[i],
                    "comparison_prediction": other[i],
                    "difference": difference,
                    "absolute_tolerance": atol,
                    "relative_tolerance": rtol,
                    "rounding_allowance": allowance[i],
                    "combined_tolerance": tolerance,
                    "status": status,
                    "measured_confirmation": False,
                    "reason": "Missing/nonfinite prediction"
                    if not finite
                    else "Equation/table arithmetic with declared coefficient-rounding allowance; not experimental confirmation",
                }
            )
    return rows


def coefficient_table_specification(table, basis, coding=None):
    """Convert explicit coefficient cells to the safe equation schema; no expression parser."""
    import pandas as pd

    if basis not in {"actual", "coded"}:
        raise InputError("Choose the actual or coded coefficient basis.")
    if not {"response", "term", "coefficient"} <= set(table):
        raise InputError("The coefficient table needs response, term and coefficient columns.")
    if table.empty or table[["response", "term"]].isna().any().any() or table.duplicated(["response", "term"]).any():
        raise InputError("Supply each response and polynomial term once.")
    values = pd.to_numeric(table.coefficient, errors="coerce")
    if not np.isfinite(values.to_numpy(float)).all():
        raise InputError("Fill every coefficient with a finite number, including explicit zero terms.")
    coding_map = {}
    if basis == "coded":
        if coding is None or not {"factor", "centre", "half_range"} <= set(coding):
            raise InputError("Coded coefficients require factor, centre and half_range columns.")
        if coding.factor.duplicated().any():
            raise InputError("Coding factors must be unique.")
        centres = pd.to_numeric(coding.centre, errors="coerce")
        widths = pd.to_numeric(coding.half_range, errors="coerce")
        if not np.isfinite(centres).all() or not np.isfinite(widths).all() or not (widths > 0).all():
            raise InputError("Fill each factor centre and a positive half-range.")
        coding_map = {
            str(f): {"centre": float(c), "half_range": float(w)} for f, c, w in zip(coding.factor, centres, widths)
        }
    output = {}
    for response in table.response.unique():
        selected = table.response.eq(response)
        output[str(response)] = {
            "basis": basis,
            "coefficients": dict(zip(table.loc[selected, "term"].astype(str), values[selected].astype(float))),
        }
        if basis == "coded":
            output[str(response)]["coding"] = coding_map
    return output
