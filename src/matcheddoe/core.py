"""Estimable response surfaces, replicated-design diagnostics and supported candidates."""

from itertools import combinations, product

import numpy as np
import pandas as pd
from scipy.spatial import Delaunay, QhullError
from scipy.stats import f as f_distribution
from scipy.stats import t as t_distribution

from .common import InputError, Result, audit, fingerprint, identifiers, numbers
from .contract import analysis_contract


def design_matrix(coded, factors, model="quadratic"):
    if model not in {"linear", "2fi", "quadratic", "additive_quadratic"}:
        raise InputError("Model must be linear, 2fi, quadratic or additive_quadratic.")
    if any(name == "Intercept" or ":" in name or "^" in name for name in factors):
        raise InputError(
            "factors: rename factor headers containing ':', '^' or the reserved polynomial term Intercept so structured term names are unambiguous."
        )
    columns, names = [np.ones(len(coded))], ["Intercept"]
    columns.extend(coded[:, j] for j in range(len(factors)))
    names.extend(factors)
    if model in {"2fi", "quadratic"}:
        for a, b in combinations(range(len(factors)), 2):
            columns.append(coded[:, a] * coded[:, b])
            names.append(f"{factors[a]}:{factors[b]}")
    if model in {"quadratic", "additive_quadratic"}:
        for j, name in enumerate(factors):
            columns.append(coded[:, j] ** 2)
            names.append(f"{name}^2")
    return np.column_stack(columns), names


def fit_surface(x, y):
    rank = np.linalg.matrix_rank(x)
    if rank < x.shape[1]:
        raise InputError(
            f"The selected model has {x.shape[1]} coefficients but design rank is {rank}. Add distinct factor combinations or choose a simpler model."
        )
    coefficients = np.linalg.lstsq(x, y, rcond=None)[0]
    fitted = x @ coefficients
    residual = y - fitted
    degrees = len(y) - rank
    sse = float(residual @ residual)
    mse = sse / degrees if degrees > 0 else None
    _, singular, right = np.linalg.svd(x, full_matrices=False)
    inverse = (right.T / singular**2) @ right
    leverage = np.einsum("ij,jk,ik->i", x, inverse, x)
    return coefficients, fitted, residual, degrees, sse, mse, inverse, leverage


@analysis_contract
def analyse(
    frame,
    factors,
    responses=None,
    id_col="sample_id",
    model="quadratic",
    directions=None,
    grid_points=21,
    independent_runs=False,
    reported_equations=None,
    response_mode="matched",
    max_design_cells=2_000_000,
    max_grid_candidates=100_000,
):
    responses = responses or ["response"]
    candidate_columns = {
        "factors": {name: "factor::" + name for name in factors},
        "predictions": {name: "prediction::" + name for name in responses},
    }
    if not factors or len(factors) != len(set(factors)):
        raise InputError("Choose one or more distinct numeric factors.")
    if response_mode not in {"matched", "per_response"}:
        raise InputError("response_mode: choose matched or per_response.")
    if any(type(v) is not int or v < 1 for v in (max_design_cells, max_grid_candidates)):
        raise InputError("Computation budgets must be positive integer cell/candidate counts.")
    term_count = (
        1
        + len(factors)
        + (len(factors) * (len(factors) - 1) // 2 if model in {"2fi", "quadratic"} else 0)
        + (len(factors) if model in {"quadratic", "additive_quadratic"} else 0)
    )
    if len(frame) * term_count > max_design_cells:
        raise InputError(
            f"This model needs {len(frame) * term_count:,} design cells, above max_design_cells={max_design_cells:,}. Choose fewer terms or explicitly increase the computation budget."
        )
    if response_mode == "per_response" and len(responses) > 1:
        result = Result(
            "MatchedDoE — independent response analyses",
            {"input_audit": audit(frame)},
            settings={
                "input_sha256": fingerprint(frame),
                "factors": factors,
                "candidate_columns": candidate_columns,
                "responses": responses,
                "model": model,
                "response_mode": response_mode,
                "independent_runs_confirmed": independent_runs,
                "response_models": {},
            },
        )
        combined = {}
        for response in responses:
            child = analyse.__wrapped__(
                frame,
                factors=factors,
                responses=[response],
                id_col=id_col,
                model=model,
                directions={response: (directions or {}).get(response, "maximise")},
                grid_points=grid_points,
                independent_runs=independent_runs,
                reported_equations={response: reported_equations[response]}
                if reported_equations and response in reported_equations
                else None,
                response_mode="matched",
                max_design_cells=max_design_cells,
                max_grid_candidates=max_grid_candidates,
            )
            result.settings["response_models"][response] = child.settings
            if child.status != "complete":
                result.status = "partial"
            for name, table in child.tables.items():
                if name == "input_audit":
                    continue
                table = table.copy()
                if "response" not in table:
                    table["response"] = response
                combined.setdefault("candidates_by_response" if name == "candidates" else name, []).append(table)
            result.notes.extend(f"{response}: {note}" for note in child.notes)
        result.tables.update(
            {name: pd.concat(parts, ignore_index=True, sort=False) for name, parts in combined.items()}
        )
        result.notes.insert(
            0,
            "Responses use their own complete rows and coding ranges. These fits are not matched-cohort comparisons; no joint Pareto frontier is calculated.",
        )
        return result
    if set(factors) & set(responses) or id_col in set(factors + responses):
        raise InputError("Factors, responses and observation ID must use different columns.")
    if not 5 <= grid_points <= 41:
        raise InputError("Grid resolution must be between 5 and 41 levels per factor.")
    ids, _ = identifiers(frame, id_col)
    result = Result(
        "MatchedDoE — response surfaces at matched conditions",
        {"input_audit": audit(frame)},
        settings={
            "input_sha256": fingerprint(frame),
            "factors": factors,
            "candidate_columns": candidate_columns,
            "responses": responses,
            "model": model,
            "response_mode": response_mode,
            "grid_points_per_factor": grid_points,
            "independent_runs_confirmed": independent_runs,
        },
    )
    result.notes.extend(
        [
            "All responses use the same complete observed runs and the same candidate conditions. Missing responses and missing factor settings are never imputed.",
            "Coefficient inference assumes independent errors with constant variance and a suitable response-surface form. Repeated factor settings are treated as independent experimental replicates only if that is how they were collected.",
            "Candidate maxima/minima are predictions on a finite grid inside the observed convex hull. They are not measured confirmations or guaranteed continuous optima.",
        ]
    )
    xf = numbers(frame, factors)
    if any(r not in frame for r in responses):
        result.status = "audit_only"
        result.notes.append(
            "One or more response columns are absent. The design audit is available; supply measured responses to fit surfaces."
        )
        return result
    yf = numbers(frame, responses)
    complete = xf.notna().all(axis=1) & yf.notna().all(axis=1)
    result.tables["excluded_runs"] = pd.DataFrame(
        {"sample_id": ids[~complete], "reason": "incomplete factor or response for matched comparison"}
    )
    result.tables["observed_factor_values"] = pd.DataFrame(
        [
            {"sample_id": ids[i], "factor": factor, "value": xf.iloc[i][factor]}
            for i in np.flatnonzero(complete)
            for factor in factors
        ]
    )
    x_raw, observed = xf[complete].to_numpy(), yf[complete]
    if len(x_raw) == 0:
        result.status = "audit_only"
        result.notes.append(
            "No complete run remains for this response selection. Choose per-response mode for incomplete secondary responses."
        )
        return result
    low, high = x_raw.min(axis=0), x_raw.max(axis=0)
    centre, half_range = (low + high) / 2, (high - low) / 2
    if (half_range == 0).any():
        result.status = "audit_only"
        result.notes.append(
            "At least one selected factor has only one observed level. Remove that factor or add varying settings."
        )
        return result
    coded = (x_raw - centre) / half_range
    x, terms = design_matrix(coded, factors, model)
    rank = int(np.linalg.matrix_rank(x))
    result.tables["factor_coding"] = pd.DataFrame(
        {"factor": factors, "low": low, "centre": centre, "high": high, "half_range": half_range}
    )
    result.tables["design_check"] = pd.DataFrame(
        [
            {
                "complete_runs": len(x),
                "unique_settings": len(np.unique(x_raw, axis=0)),
                "coefficients": x.shape[1],
                "rank": rank,
                "residual_df": len(x) - rank,
                "condition_number": float(np.linalg.cond(x)),
            }
        ]
    )
    support = []
    for family in ["linear", "2fi", "additive_quadratic", "quadratic"]:
        candidate_design, _ = design_matrix(coded, factors, family)
        family_rank = int(np.linalg.matrix_rank(candidate_design))
        support.append(
            {
                "model": family,
                "terms": candidate_design.shape[1],
                "rank": family_rank,
                "estimable": family_rank == candidate_design.shape[1],
            }
        )
    result.tables["model_support"] = pd.DataFrame(support)
    if rank < x.shape[1]:
        result.status = "audit_only"
        result.notes.append(
            f"The {model} model needs {x.shape[1]} independent coefficients but rank is {rank}. Choose a simpler model or add distinct design settings."
        )
        return result
    if not independent_runs:
        result.status = "partial"
        result.notes.append(
            "Independence of experimental runs was not confirmed. Fitted surfaces remain descriptive; coefficient uncertainty, prediction intervals, pure error and lack-of-fit inference are withheld. Repeated instrument readings do not establish independent runs."
        )
    if len(factors) > 3:
        grid = np.empty((0, len(factors)))
        result.notes.append(
            "More than three factors: coefficients, fitted values and residuals are available. Automatic multidimensional candidate optimization is disabled to avoid an exponential grid; no global optimum is claimed."
        )
        result.settings["candidate_status"] = "unavailable: more than three factors; no Cartesian grid generated"
    elif grid_points ** len(factors) > max_grid_candidates:
        grid = np.empty((0, len(factors)))
        result.notes.append(
            "Candidate grid exceeds the selected computation budget. Fitted-model calculations remain available; reduce grid resolution or raise max_grid_candidates explicitly."
        )
        result.settings["candidate_status"] = "unavailable: grid budget exceeded"
    else:
        grid = np.array(list(product(np.linspace(-1, 1, grid_points), repeat=len(factors))))
        if len(factors) == 1:
            supported = (grid[:, 0] >= coded[:, 0].min()) & (grid[:, 0] <= coded[:, 0].max())
        else:
            try:
                supported = Delaunay(np.unique(coded, axis=0)).find_simplex(grid, tol=1e-9) >= 0
            except QhullError:
                supported = np.zeros(len(grid), dtype=bool)
                result.notes.append(
                    "Observed settings do not span a full-dimensional convex hull. Candidate optimization is unavailable."
                )
        grid = grid[supported]
        result.settings["candidate_status"] = (
            "supported finite grid" if len(grid) else "unavailable: no supported candidate"
        )
    grid_x, _ = design_matrix(grid, factors, model)
    grid_real = grid * half_range + centre
    candidate_table = pd.DataFrame(grid_real, columns=[candidate_columns["factors"][f] for f in factors])
    result.tables["candidate_column_labels"] = pd.DataFrame(
        [
            {"role": role, "original_label": name, "export_column": column}
            for role, mapping in candidate_columns.items()
            for name, column in mapping.items()
        ]
    )
    coefficients, summaries, predictions, chosen, anova_rows, equation_rows = [], [], [], [], [], []
    _, setting_code = np.unique(x_raw, axis=0, return_inverse=True)
    if reported_equations is not None and (
        not isinstance(reported_equations, dict) or set(reported_equations) - set(responses)
    ):
        raise InputError("reported_equations: map selected response names to structured equations.")
    for response in responses:
        y = observed[response].to_numpy()
        coef, fitted, residual, df, sse, mse, inverse, leverage = fit_surface(x, y)
        sst = float(np.sum((y - y.mean()) ** 2))
        pe = sum(
            float(np.sum((y[setting_code == k] - y[setting_code == k].mean()) ** 2)) for k in np.unique(setting_code)
        )
        pe_df = len(y) - len(np.unique(setting_code))
        lof_df = len(np.unique(setting_code)) - rank
        lof_ss = max(0.0, sse - pe)
        lof_f = (lof_ss / lof_df) / (pe / pe_df) if pe_df > 0 and lof_df > 0 and pe > 0 else None
        # A residual at arithmetic roundoff does not provide a resolved variance estimate.
        arithmetic_scale = np.linalg.norm(x) * np.linalg.norm(coef) + np.linalg.norm(y)
        resolved_variance = df > 0 and mse is not None and sse > (8 * np.finfo(float).eps * arithmetic_scale) ** 2
        critical = float(t_distribution.ppf(0.975, df)) if resolved_variance and independent_runs else None
        if independent_runs and not resolved_variance:
            result.status = "partial"
            result.notes.append(
                f"{response}: residual variance is unavailable or indistinguishable from arithmetic roundoff. F tests and uncertainty intervals are withheld; coefficients and descriptive residuals remain available."
            )
        standard_errors = (
            np.sqrt(np.maximum(0, np.diag(inverse) * mse))
            if resolved_variance and independent_runs
            else np.full(len(coef), np.nan)
        )
        for term, value, se in zip(terms, coef, standard_errors):
            coefficients.append(
                {
                    "response": response,
                    "term": term,
                    "coefficient_coded": value,
                    "standard_error": se,
                    "ci95_low": value - critical * se if critical is not None else np.nan,
                    "ci95_high": value + critical * se if critical is not None else np.nan,
                }
            )
        press = float(np.sum((residual / (1 - leverage)) ** 2)) if np.all(leverage < 1 - 1e-10) else None
        summaries.append(
            {
                "response": response,
                "n": len(y),
                "r2_fitted": 1 - sse / sst if sst > 0 else None,
                "adjusted_r2": 1 - (sse / df) / (sst / (len(y) - 1)) if df > 0 and sst > 0 else None,
                "residual_ss": sse,
                "residual_df": df,
                "pure_error_ss": pe if pe_df > 0 else None,
                "pure_error_df": pe_df if independent_runs else None,
                "lack_of_fit_ss": lof_ss if pe_df > 0 and independent_runs else None,
                "lack_of_fit_df": lof_df if pe_df > 0 and independent_runs else None,
                "lack_of_fit_f": lof_f if independent_runs else None,
                "lack_of_fit_p": float(f_distribution.sf(lof_f, lof_df, pe_df))
                if lof_f is not None and independent_runs
                else None,
                "leave_one_run_out_press": press,
            }
        )
        model_df = rank - 1
        model_ss = max(0.0, sst - sse)
        model_f = (
            (model_ss / model_df) / mse if independent_runs and model_df > 0 and resolved_variance and sst > 0 else None
        )
        model_p = float(f_distribution.sf(model_f, model_df, df)) if model_f is not None else None
        for component_name, ss_value, degrees_value, ms_value, f_value, p_value in [
            ("model", model_ss, model_df, model_ss / model_df if model_df else None, model_f, model_p),
            ("residual", sse, df, mse, None, None),
            ("total_corrected", sst, len(y) - 1, sst / (len(y) - 1) if len(y) > 1 else None, None, None),
        ]:
            anova_rows.append(
                {
                    "response": response,
                    "source": component_name,
                    "ss": ss_value,
                    "df": degrees_value,
                    "ms": ms_value,
                    "f": f_value,
                    "p": p_value,
                    "ss_convention": "overall regression with intercept; no term-level SS",
                    "inference": "available under independent-run assumptions"
                    if model_f is not None
                    else "unavailable: independence or positive residual variance/df required",
                }
            )
        if reported_equations and response in reported_equations:
            from .equations import compare_equation

            equation_rows.extend(
                compare_equation(
                    frame.loc[complete].reset_index(drop=True),
                    factors,
                    model,
                    reported_equations[response],
                    response,
                    ids[complete],
                    fitted,
                )
            )
        if not independent_runs:
            summaries[-1]["pure_error_ss"] = None
        predictions.extend(
            {
                "sample_id": ids[np.flatnonzero(complete)[i]],
                "response": response,
                "reference": y[i],
                "prediction": fitted[i],
                "residual": fitted[i] - y[i],
                "leverage": leverage[i],
                "prediction_kind": "fitted_on_same_runs",
            }
            for i in range(len(y))
        )
        candidate_prediction = grid_x @ coef
        candidate_table[candidate_columns["predictions"][response]] = candidate_prediction
        if len(grid):
            direction = (directions or {}).get(response, "maximise")
            if direction not in {"maximise", "minimise"}:
                raise InputError("Response direction must be maximise or minimise.")
            position = int(
                np.argmax(candidate_prediction) if direction == "maximise" else np.argmin(candidate_prediction)
            )
            point = grid_x[position]
            se_pred = np.sqrt(mse * (1 + point @ inverse @ point)) if mse is not None else np.nan
            chosen.append(
                {
                    "response": response,
                    "direction": direction,
                    **dict(zip(candidate_columns["factors"].values(), grid_real[position])),
                    "predicted_response": candidate_prediction[position],
                    "prediction_interval95_low": candidate_prediction[position] - critical * se_pred
                    if critical is not None
                    else np.nan,
                    "prediction_interval95_high": candidate_prediction[position] + critical * se_pred
                    if critical is not None
                    else np.nan,
                    "inside_observed_convex_hull": True,
                    "measured_confirmation": False,
                }
            )
    result.tables.update(
        coefficients=pd.DataFrame(coefficients),
        model_summary=pd.DataFrame(summaries),
        anova=pd.DataFrame(anova_rows),
        equation_comparison=pd.DataFrame(equation_rows),
        predictions=pd.DataFrame(predictions),
        candidates=candidate_table,
        selected_candidates=pd.DataFrame(chosen),
    )
    if len(x) - rank == 0:
        result.notes.append(
            "The model is saturated: coefficients are estimable but no residual variance or prediction intervals can be estimated."
        )
    if len(np.unique(setting_code)) == len(x):
        result.notes.append("No repeated setting is present, so pure error and lack-of-fit tests are unavailable.")
    result.notes.append(
        "Pointwise prediction intervals are conditional on the fitted polynomial and normal-error assumptions; intervals at selected optima do not account for the search or model selection."
    )
    return result
