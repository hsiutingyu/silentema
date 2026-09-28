# silentema

<!-- badges: start -->
[![R-CMD-check](https://github.com/hsiutingyu/silentema/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/hsiutingyu/silentema/actions/workflows/R-CMD-check.yaml)
[![CRAN status](https://www.r-pkg.org/badges/version/silentema)](https://CRAN.R-project.org/package=silentema)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
<!-- badges: end -->

Dynamic missingness graphs, recoverability checks, and sensitivity analysis for informative nonresponse in ecological momentary assessment (EMA) and other experience-sampling data.

Participants skip prompts, and the reasons are rarely unrelated to the states being measured. `silentema` lets an analyst

1. **declare** the assumed nonresponse mechanism as a dynamic missingness graph (`dm_graph()`), built from a taxonomy of seven motifs (completely random, lagged-state dependence, self-censoring, burden, person propensity, context confounding, and reactivity);
2. **check** by d-separation which estimands of a two-level VAR(1) (temporal network, contemporaneous network, person means, between-person law) remain structurally recoverable and by which estimator (`recoverability()`);
3. **test** whether skipped prompts were informative, using the state after a skipped prompt (`silence_test()`) or an always-observed sensor (`sensor_gap_test()`), and describe response persistence (`fatigue_check()`);
4. **estimate** from answered adjacent prompts with person-specific intercepts and a half-panel jackknife (`fit_pairs()`), by inverse-probability weighting on an observed context (`fit_ipw()`), or by full-information maximum likelihood under missing at random (MAR) (`fit_fiml()`);
5. **profile** the estimates over a self-censoring sensitivity value by inverse-probability weighting with a fixed probit selection model whose intercept is calibrated to the response rate ("tilting") (`fit_tilt()`, `tilt_profile()`, `break_even()`);
6. **calibrate** the sensitivity value from a passive sensor, randomized probes, or the post-skip contrast (`calibrate_delta()`), and compare with worst-case bounds (`bounds_support()`);
7. **report** with a preregistration-ready missingness declaration (`missingness_declaration()`).

A simulator for the whole taxonomy (`simulate_ema()`, `simulate_from_fit()`) supports design planning and replication.

## Installation

From CRAN (once the package is accepted):

```r
install.packages("silentema")
```

The development version from GitHub:

```r
# install.packages("remotes")
remotes::install_github("hsiutingyu/silentema")
```

The package contains C++ code (the state-space EM algorithm), so a compiler toolchain is needed to install from source: Rtools on Windows, Xcode command-line tools on macOS, `r-base-dev` or equivalent on Linux.

## Minimal example

```r
library(silentema)
sim <- simulate_ema(N = 100, n_prompts = 56, motifs = "M2", compliance = 0.7, delta = -1,
                    sensor_cor = 0.6, seed = 1)
g <- dm_graph("M2", sensor = TRUE)
recoverability(g)                       # what can be recovered under the declared mechanism?
silence_test(sim$data, sim$vars)        # was silence informative?
prof <- tilt_profile(sim$data, sim$vars, delta_grid = seq(-2, 0.5, by = 0.5))
cal <- calibrate_delta(prof, sim$data, method = "sensor")
break_even(prof, delta_max = 1.5)
cat(missingness_declaration(g, profile = prof, calibration = cal, plausible = c(-1.5, 0)), sep = "\n")
```

## Documentation

* `vignette("silentema-workflow")`: a complete analysis from declaration to report.
* `vignette("dm-graphs")`: the motif taxonomy, d-separation and the recoverability report.
* `vignette("testing-informativeness")`: the silence test, the sensor-gap test and the fatigue check.
* `vignette("sensitivity-analysis")`: tilting, break-even values, calibration and reporting.
* `vignette("simulation-and-design")`: the simulator and design planning.
* Documentation site: https://hsiutingyu.github.io/silentema/

## Citation

Yu, H.-T. (2026). What skipped prompts hide: Detecting, diagnosing, and correcting informative nonresponse in ecological momentary assessment. Manuscript under review. Materials: https://osf.io/x6d2t/

Yu, H.-T. (2026). silentema: Dynamic Missingness Graphs and Sensitivity Analysis for EMA Data. R package version 1.0.0. https://CRAN.R-project.org/package=silentema

Run `citation("silentema")` in R for BibTeX entries.

## License

GPL (>= 3)
