# Methods and interpretation

## Question and experimental unit

Which response surfaces are estimable from the observed runs, and how do different responses compare at the same supported factor settings? The program fits one or more factors with one or several continuous responses. Matched mode uses the same complete observed runs for every selected response [P26–P30, P33].

An observation ID identifies one row. A repeated factor setting is a genuine experimental replicate only if a separate experimental run was performed. Repeated readings of one preparation do not establish pure error between experimental runs. In the guided interface, independence is initially unconfirmed for uploaded files. Descriptive fits remain available, but coefficient uncertainty, prediction intervals and pure-error/lack-of-fit inference are withheld until independent runs are explicitly identified. The bundled design examples use the deposited run structure.

## Capability rules

One or more factors are accepted. The selected model must have full column rank; model_support.csv lists existing model forms without silently choosing another. Saturated designs retain coefficients and fitted values without residual-variance inference. More than three factors retain fits and residuals, with automatic Cartesian candidate search disabled. Grid and design-matrix budgets are explicit. Matched mode uses the same complete runs for every response; per-response mode fits each available response on its own cohort and coding and disables joint Pareto claims. Reported equations can be entered as finite coefficient tables in actual or explicitly coded units. Expressions are never executed.

## Model construction

For each factor, coding is (observed value − midpoint)/half-range, with midpoint and half-range calculated from the complete matched runs. Physical units and coding parameters are retained. A factor with no variation yields an audit. Available forms are linear, 2FI (main effects and all pairwise interactions), additive quadratic (main effects and squares), and full quadratic (main effects, pairwise interactions and squares). Main effects are retained whenever interactions or squares are present. The model form is an explicit choice; no automatic term-by-term p-value selection is performed.

Coefficients are calculated by least squares. Rank is checked against the number of coefficients before inference. Covariance geometry and leverage are calculated from a singular-value decomposition, avoiding inversion of the squared-condition normal matrix. Residual sum of squares (SSE) is the sum of squared measured-minus-fitted differences. Residual degrees of freedom are n − rank. Coefficient standard errors use residual mean square and assume independent errors of constant variance. Student-t intervals additionally assume an appropriate normal-error model. A saturated model has coefficients but no residual-variance estimate.

When independent replication is confirmed, pure-error SSE is the sum of squared deviations from each identical nominal factor setting's response mean. Pure-error degrees of freedom equal n − number of distinct settings. Lack-of-fit SSE is residual SSE minus pure-error SSE; lack-of-fit degrees of freedom equal number of distinct settings minus model rank. An F statistic is reported only when both relevant degrees of freedom and a positive pure-error variance exist. Exact numeric nominal settings define replicates; nearly equal measured values are not silently rounded into replication [P26–P27].

The analytic PRESS statistic uses residual/(1 − leverage) for each omitted run and is unavailable at unit leverage. It is a leave-one-run-out diagnostic for the fixed model, not a new-setting validation when other replicates of that setting remain in training. Fitted and adjusted R² are descriptive fit statistics, not independently validated accuracy.

## Matched candidate conditions

A fixed grid, by default 21 levels per factor, is generated in coded coordinates and restricted to the convex hull of observed complete settings. Every response is predicted at the same retained conditions. The table gives each response's highest or lowest grid prediction according to its declared objective. Separate response optima are not automatically a joint optimum; the common candidate table makes tradeoffs visible.

For Box–Behnken designs, cube corners outside the observed hull are excluded. A finite grid supplies candidate conditions, not a guaranteed continuous optimum. A pointwise 95% prediction interval, when estimable, combines new-observation variance and coefficient uncertainty. Intervals at selected extrema do not account for the search or for model-form selection. A candidate is never labelled a measured confirmation [P26–P29].

## Examples and limits

The 15-run LEVO example uses a full 2FI surface. Our residual SSE, pure-error SSE and lack-of-fit arithmetic agree with the deposited design table. The paper's separate confirmation claim conflicts with the deposited confirmation sheet and is excluded. Its reported optimum is at an unsupported cube corner; it is not imported as our confirmed optimum [P30].

The banana example contains nine distinct drying settings and two responses. It uses additive quadratic surfaces and supplies no independent replication at identical settings, so pure error is unavailable. Detailed deposited fit calculations are used rather than the inconsistent summary R² labels [P33]. Missing cells exclude a run from the affected complete cohort: all selected responses in matched mode, or each response separately in per-response mode. They are never imputed. A polynomial may be estimable yet scientifically unsuitable. Inspect residuals, condition number and the recorded design before acting on candidate settings.

## Full model ANOVA, stated equations and influence

The ANOVA export partitions corrected total SS into model and residual SS, with degrees of freedom, mean squares and the overall model F test where its assumptions and variance support it. It does not silently select type-I/II/III term-wise sums of squares. A saturated or roundoff-perfect fit withholds residual-variance inference.

A stated equation is a structured coefficient object in explicitly actual or coded coordinates. Complete term names, coding constants, numeric tolerances and a coefficient-rounding bound are required as appropriate. The equation is compared at matched observed rows with the refitted predictions and, optionally, a separate reported-prediction column. No executable expression is evaluated. See [Equation check](EQUATION_CHECK.md).

Cook's distance uses e²/(p MSE) × h/(1−h)²; internal studentization uses e/sqrt(MSE(1−h)). These are descriptive review diagnostics and are withheld when residual variance or leverage cannot support them. For exactly two responses, the nondominated candidate set uses the declared objective directions on the existing supported grid. It does not select a single joint optimum or create confirmation data.


## Candidate-table schema

User factor labels occupy `factor::<original label>` columns and predicted responses occupy `prediction::<original label>` columns. These disjoint namespaces prevent metadata collisions, including a factor literally called `response`, `direction` or `predicted_response`. The exact original label and its role are retained in `candidate_column_labels.csv` and `settings.candidate_columns`. The interface labels these as Factor and Prediction without changing the original spelling.

`selected_candidates` retains response identity, direction, physical factor settings, fitted prediction, interval bounds, hull support and confirmation status. All candidate tables, per-response tables and Pareto diagnostics use the same mapping. No input header is forbidden because it matches output metadata. Duplicate input headers are rejected before analysis. The structured polynomial equation syntax still requires unambiguous term names as documented in EQUATION_CHECK.md.
