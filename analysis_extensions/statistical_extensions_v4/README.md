# Statistical extensions v4

This directory was generated from the schema-locked, project-coded AADR-
derived analytical input included with the release. The input SHA-256 is
`692a69cf38cc736f96ea5aa6b3b15024a9a49d50c03ae3e1470a8dc475504cc3`.

This run used seed 20260821, 9,999 cluster-wild resamples,
10,000 parametric draws per transition, and
10,000 paired-cluster bootstrap replicates. The frozen defaults
are seed 20260821, 9,999, 10,000, and
10,000, respectively. Runs with fewer replicates are for
execution checks and should not replace the frozen analysis for inference.

The input is project-coded rather than anonymous because combinations of
public AADR-derived attributes may remain linkable to the source resource. No
source person identifier, locality string, skeletal field, exact coordinate or
terminal haplogroup call is required by this extension.

## Decision rules

* **Revised primary profile estimand:** the country-adjusted difference among
  equally weighted site-period lineage profiles containing at least two calls.
  This excludes indicator-vector singleton profiles while retaining a usable
  number of localities. It is an archive/site-profile estimand, not a Central
  Asian population-frequency estimand.
* **Stringent sensitivity:** profiles containing at least three calls.
* The all-profile analysis remains an archival-scope context analysis so the
  revision does not erase the originally analysed record.
* Cluster-wild tests use a null-imposed HC2 Rademacher sign shared across all
  time profiles of a country-locality. Holm adjustment is applied across the
  two markers within each threshold. Because the thresholds are nested, the
  >=3 result is interpreted as sensitivity rather than a second independent
  confirmatory test.

## Stable-profile results

| marker | minimum_calls | analysis_role | n_analyzed_profiles | n_calls_in_analyzed_profiles | n_sites | partial_r2 | raw_cluster_wild_p | holm_within_threshold_p |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| mtDNA | 1 | archival-scope context | 170 | 438 | 132 | 0.082533 | 0.0001 | 0.0002 |
| mtDNA | 2 | revised primary | 71 | 339 | 60 | 0.120017 | 0.2224 | 0.2224 |
| mtDNA | 3 | stringent sensitivity | 40 | 277 | 38 | 0.232995 | 0.4203 | 0.4203 |
| Y | 1 | archival-scope context | 118 | 229 | 101 | 0.130205 | 0.0009 | 0.0009 |
| Y | 2 | revised primary | 48 | 156 | 41 | 0.287234 | 0.007 | 0.014 |
| Y | 3 | stringent sensitivity | 22 | 104 | 21 | 0.570548 | 0.036 | 0.072 |

The all-profile Y row is an extension rerun with a different fixed seed. The
frozen 25 July 2026 analysis reported raw/Holm *P*=0.0012; the extension gave
0.0009 under the frozen defaults. The result of the current run is shown above.
The discussion below describes the frozen-default analysis; changing the seed
or resampling counts requires reviewing the current tables before using it.

The mtDNA period association in the all-profile archive is not reproduced
after singleton exclusion. The Y result persists at >=2 calls, while the >=3
analysis has few residual degrees of freedom and should not be overinterpreted.

## Why no claimed hierarchical Dirichlet-multinomial result

The requested site/publication multilevel specification is structurally weak:
most sites occur in only one period, so a full site intercept is almost nested
with period. The identifiability audit is:

| marker | minimum_calls | n_profiles | n_sites | n_multiperiod_sites | proportion_single_period_sites | n_profile_publication_combinations | n_publications_after_splitting_combinations | rank_nuisance_country_publication_site | rank_increment_period_after_site_publication | nominal_period_df | interpretation |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| mtDNA | 1 | 170 | 132 | 31 | 0.7652 | 19 | 13 | 138 | 7 | 7 | A full site-intercept model identifies period chiefly from the small set of multi-period sites; site and period are otherwise nested. |
| mtDNA | 2 | 71 | 60 | 8 | 0.8667 | 13 | 10 | 63 | 3 | 7 | A full site-intercept model identifies period chiefly from the small set of multi-period sites; site and period are otherwise nested. |
| mtDNA | 3 | 40 | 38 | 1 | 0.9737 | 12 | 10 | 39 | 1 | 7 | A full site-intercept model identifies period chiefly from the small set of multi-period sites; site and period are otherwise nested. |
| Y | 1 | 118 | 101 | 14 | 0.8614 | 18 | 12 | 105 | 6 | 7 | A full site-intercept model identifies period chiefly from the small set of multi-period sites; site and period are otherwise nested. |
| Y | 2 | 49 | 42 | 5 | 0.881 | 13 | 10 | 44 | 2 | 7 | A full site-intercept model identifies period chiefly from the small set of multi-period sites; site and period are otherwise nested. |
| Y | 3 | 23 | 22 | 1 | 0.9545 | 10 | 9 | 22 | 1 | 7 | A full site-intercept model identifies period chiefly from the small set of multi-period sites; site and period are otherwise nested. |

Fitting a nominal mixed model would not manufacture independent longitudinal
information. Instead, `multinomial_predictive_sensitivity.csv` reports a
nested site-group cross-validation analysis using L2-shrunk country,
publication and site intercept terms, with and without period. A ridge penalty
is a Gaussian-prior/MAP analogue, but this remains a predictive sensitivity,
not a Dirichlet-multinomial hypothesis test.

| marker | minimum_calls | valid | n_individual_calls | n_site_period_profiles | n_sites | n_lineage_categories | reduced_heldout_site_log_loss | full_period_heldout_site_log_loss | delta_log_loss_full_minus_reduced | site_bootstrap_delta_ci_low | site_bootstrap_delta_ci_high | outer_group_folds | model_interpretation |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| mtDNA | 2 | True | 339 | 71 | 60 | 17 | 2.6698 | 2.6931 | 0.023296 | 0.00077581 | 0.046271 | 4 | Negative delta favors period for held-out-site prediction. Ridge-penalized site/publication terms are nuisance intercept proxies, not a fitted random-effects sampling model; interval is descriptive site bootstrap, not a P value. |
| mtDNA | 3 | True | 277 | 40 | 38 | 17 | 2.6437 | 2.6473 | 0.003678 | -0.016797 | 0.026931 | 4 | Negative delta favors period for held-out-site prediction. Ridge-penalized site/publication terms are nuisance intercept proxies, not a fitted random-effects sampling model; interval is descriptive site bootstrap, not a P value. |
| Y | 2 | True | 160 | 49 | 42 | 9 | 1.66 | 1.6079 | -0.052084 | -0.10693 | 0.0044344 | 4 | Negative delta favors period for held-out-site prediction. Ridge-penalized site/publication terms are nuisance intercept proxies, not a fitted random-effects sampling model; interval is descriptive site bootstrap, not a P value. |
| Y | 3 | True | 108 | 23 | 22 | 9 | 1.7542 | 1.6534 | -0.1008 | -0.16252 | -0.038601 | 4 | Negative delta favors period for held-out-site prediction. Ridge-penalized site/publication terms are nuisance intercept proxies, not a fitted random-effects sampling model; interval is descriptive site bootstrap, not a P value. |

Negative full-minus-reduced log loss means period improved prediction for
held-out sites. The site-bootstrap interval is descriptive and is not a P
value.

## Finite-sample TV sensitivity

For each adjacent transition and marker, two corrections are supplied:

1. `null_expected_tv_noise_floor`: expected positive TV when the two periods
   share one Jeffreys-smoothed composition but retain their observed per-site
   sample sizes.
2. `working_model_estimated_bias`: parametric plug-in bias when the periods
   have their own Jeffreys-smoothed compositions.

Neither correction is assumption-free. In particular, multinomial simulation
does not include extra site heterogeneity, and subtracting a null noise floor
does not generally produce an unbiased distance. Therefore corrected values
are sensitivities, not replacements for the raw estimator.

| y_encoding | adjustment | mt_point_mean_tv | y_point_mean_tv | delta_point_y_minus_mt | bootstrap_delta_median_diagnostic | bootstrap_delta_mean_minus_point | percentile_delta_ci_low | percentile_delta_ci_high | basic_delta_ci_low | basic_delta_ci_high | median_centered_delta_ci_low | median_centered_delta_ci_high | cluster_bootstrap_replicates | rejected_empty_period_draws | total_attempted_cluster_draws | rejected_empty_period_fraction | interval_caveat |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| broad_L1 | raw | 0.5174 | 0.41102 | -0.10638 | -0.1314 | -0.024889 | -0.26366 | 0.0016838 | -0.21445 | 0.050899 | -0.23864 | 0.0267 | 10000 | 45 | 10045 | 0.0044798 | Observed-data values are the point estimates. Percentile, basic and median-centered intervals are all reported because sparse-bin ratio statistics make the bootstrap distribution off-center. Empty-period draws are rejected; the reported fraction quantifies this conditioning. Correction terms are fixed, so their model uncertainty is not propagated. |
| broad_L1 | noise_floor_subtracted | 0.096028 | 0.11975 | 0.023726 | 0.008279 | -0.014912 | -0.105 | 0.12652 | -0.079066 | 0.15245 | -0.089549 | 0.14197 | 10000 | 45 | 10045 | 0.0044798 | Observed-data values are the point estimates. Percentile, basic and median-centered intervals are all reported because sparse-bin ratio statistics make the bootstrap distribution off-center. Empty-period draws are rejected; the reported fraction quantifies this conditioning. Correction terms are fixed, so their model uncertainty is not propagated. |
| broad_L1 | parametric_bias_corrected | 0.2762 | 0.29799 | 0.021784 | -0.0033134 | -0.025164 | -0.13456 | 0.12869 | -0.085119 | 0.17813 | -0.10946 | 0.15379 | 10000 | 45 | 10045 | 0.0044798 | Observed-data values are the point estimates. Percentile, basic and median-centered intervals are all reported because sparse-bin ratio statistics make the bootstrap distribution off-center. Empty-period draws are rejected; the reported fraction quantifies this conditioning. Correction terms are fixed, so their model uncertainty is not propagated. |
| AADR_ISOGG_prefix_family | raw | 0.5174 | 0.50611 | -0.011292 | -0.037328 | -0.02621 | -0.15349 | 0.076679 | -0.099264 | 0.1309 | -0.12745 | 0.10271 | 10000 | 46 | 10046 | 0.0045789 | Observed-data values are the point estimates. Percentile, basic and median-centered intervals are all reported because sparse-bin ratio statistics make the bootstrap distribution off-center. Empty-period draws are rejected; the reported fraction quantifies this conditioning. Correction terms are fixed, so their model uncertainty is not propagated. |
| AADR_ISOGG_prefix_family | noise_floor_subtracted | 0.096028 | 0.1459 | 0.04987 | 0.042956 | -0.0066929 | -0.062461 | 0.14795 | -0.048207 | 0.1622 | -0.055547 | 0.15486 | 10000 | 46 | 10046 | 0.0045789 | Observed-data values are the point estimates. Percentile, basic and median-centered intervals are all reported because sparse-bin ratio statistics make the bootstrap distribution off-center. Empty-period draws are rejected; the reported fraction quantifies this conditioning. Correction terms are fixed, so their model uncertainty is not propagated. |
| AADR_ISOGG_prefix_family | parametric_bias_corrected | 0.2762 | 0.36992 | 0.093719 | 0.067358 | -0.026565 | -0.048531 | 0.18089 | 0.0065428 | 0.23597 | -0.02217 | 0.20726 | 10000 | 46 | 10046 | 0.0045789 | Observed-data values are the point estimates. Percentile, basic and median-centered intervals are all reported because sparse-bin ratio statistics make the bootstrap distribution off-center. Empty-period draws are rejected; the reported fraction quantifies this conditioning. Correction terms are fixed, so their model uncertainty is not propagated. |

The previously displayed surrogate ordinary-bootstrap probabilities are
intentionally not reproduced because they are not null-hypothesis P values.
Interpretation should be based on effect estimates, intervals and sensitivity
across Y encodings.

## Time-grid sensitivity

The manuscript bins span unequal durations. `equal_500y_*` recomputes
site-balanced composition and adjacent TV in ten equal 500-y bins. Summary:

| marker | minimum_calls | mean_tv | max_tv |
| --- | --- | --- | --- |
| Y | 1 | 0.41031 | 0.7 |
| Y | 2 | 0.54642 | 0.95 |
| mtDNA | 1 | 0.43374 | 0.63492 |
| mtDNA | 2 | 0.43582 | 0.64339 |

`original_bin_linear_scaled_tv.csv` also exposes TV divided by midpoint gap,
but labels it explicitly as a linear-scaling diagnostic rather than a
biological rate. `continuous_time_cluster_tests.csv` tests a country-adjusted
linear and flexible spline association between profile year and composition:

| marker | minimum_calls | effect | n_analyzed_profiles | partial_r2 | raw_cluster_wild_p | holm_within_minimum_and_time_basis_p | valid |
| --- | --- | --- | --- | --- | --- | --- | --- |
| mtDNA | 1 | linear year | 170 | 0.029013 | 0.0001 | 0.0002 | True |
| mtDNA | 1 | cubic B-spline year | 170 | 0.074262 | 0.0001 | 0.0002 | True |
| mtDNA | 2 | linear year | 71 | 0.034875 | 0.0062 | 0.0124 | True |
| mtDNA | 2 | cubic B-spline year | 71 | 0.10359 | 0.0211 | 0.0422 | True |
| mtDNA | 3 | linear year | 40 | 0.04091 | 0.325 | 0.65 | True |
| mtDNA | 3 | cubic B-spline year | 40 | 0.18727 | 0.2553 | 0.3968 | True |
| Y | 1 | linear year | 118 | 0.017103 | 0.1574 | 0.1574 | True |
| Y | 1 | cubic B-spline year | 118 | 0.080026 | 0.0388 | 0.0388 | True |
| Y | 2 | linear year | 48 | 0.026027 | 0.5534 | 0.5534 | True |
| Y | 2 | cubic B-spline year | 48 | 0.20434 | 0.0257 | 0.0422 | True |
| Y | 3 | linear year | 22 | 0.095974 | 0.3266 | 0.65 | True |
| Y | 3 | cubic B-spline year | 22 | 0.41454 | 0.1984 | 0.3968 | True |

These checks support only statements about archive-indexed compositional
dissimilarity/association. They do not identify a constant lineage turnover
velocity, migration rate, or sex-biased demographic mechanism.

## Files

* `stable_profile_cluster_tests.csv` and `stable_profile_vectors.csv`
* `hierarchical_model_identifiability_audit.csv`
* `multinomial_predictive_sensitivity.csv`
* `paired_male_tv_finite_sample_adjustments.csv`
* `paired_male_tv_adjusted_cluster_bootstrap_summary.csv`
* `paired_male_tv_adjusted_cluster_bootstrap_draws_*.csv`
* `equal_500y_composition.csv`, `equal_500y_turnover.csv`
* `original_bin_linear_scaled_tv.csv`
* `continuous_time_cluster_tests.csv`
* `run_manifest.json`
