# Recorded GitHub Actions result

Observed in Chrome on 2026-09-09: [scientific-checks run 34273166625](https://github.com/Negarziraksaz/matcheddoe/actions/runs/34273166625) shows **Success**, with four completed jobs and four evidence artifacts.

- Repository: `Negarziraksaz/matcheddoe`
- Exact commit: `28797d7f79023ecfd330d60e57111909f7aca18f`
- Matrix: `ubuntu-latest` / `windows-latest`, each with Python `3.12` / `3.14`.
- Workflow: dependency installation, `pip check`, Ruff lint/format, complete pytest suite, package build and installed-package smoke script.

This is evidence for that committed snapshot. Exact resolved runner/interpreter/dependency versions and test counts belong to its logs/artifacts; they are not inferred from the matrix labels. CI is not an interactive browser-hosting test, native macOS validation or independent scientific validation.

The earlier local closeout records correctly marked remote CI and local Linux as NOT_RUN at that time. This later cloud-run observation adds evidence without changing those historical records. The README badge links to current workflow status rather than embedding a permanent PASS label.
