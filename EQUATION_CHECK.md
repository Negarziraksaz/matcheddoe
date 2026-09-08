# Checking a stated equation

Open **Optional reported-equation check** or set `reported_equations` in the JSON configuration. This example explains the format with illustrative coefficients; it does not attribute an equation to a paper:

```json
{
  "response": {
    "basis": "coded",
    "coding": {
      "dose_g_L": {"centre": 1.0, "half_range": 0.5},
      "concentration_ppm": {"centre": 20.0, "half_range": 10.0},
      "pH": {"centre": 7.0, "half_range": 2.0}
    },
    "coefficients": {
      "Intercept": 50.0, "dose_g_L": 2.0, "concentration_ppm": -1.0, "pH": 3.0,
      "dose_g_L:concentration_ppm": 0.0, "dose_g_L:pH": 0.0, "concentration_ppm:pH": 0.0
    },
    "atol": 1e-8, "rtol": 1e-6, "coefficient_rounding": 0.0005
  }
}
```

Enter the actual factor names, stated coefficients and coding from the equation you intend to check. Use `basis: "actual"` and omit `coding` for an equation expressed directly in physical factor values. The selected model determines the exact required terms, including explicit zero coefficients. An optional `prediction_column` names a separate published/reported prediction column, never the measured response.

At each observed run, the allowed difference is `atol + rtol × abs(comparison prediction) + coefficient_rounding × sum(abs(polynomial terms))`. The rounding bound is absolute per coefficient in the stated basis. Missing predictions are unavailable. `equation_comparison.csv` retains equation values, comparison origin, differences, allowances and discrepancy flags. A numerical discrepancy asks for investigation of basis, rounding, model specification or transcription; agreement does not verify experimental confirmation.
