# silentema 1.0.0

First CRAN release. The statistical code is unchanged from version 0.2.1, so
every estimate produced with 0.2.0 or 0.2.1 is reproduced exactly; the release
prepares the package for CRAN, completes its documentation and fixes one
display issue.

* New vignettes: `dm-graphs` (the motif taxonomy, d-separation and
  recoverability), `testing-informativeness` (the silence test, the
  sensor-gap test and the fatigue check), `sensitivity-analysis` (tilting,
  break-even values, the three calibration designs and the missingness
  declaration) and `simulation-and-design` (the simulator, every motif's
  parameters, design planning). The workflow vignette is revised.
* Every exported function documents its return value in full and has
  runnable examples; the documentation cites the methods (Mohan & Pearl,
  2021; Dhaene & Jochmans, 2015; Shumway & Stoffer, 1982; Cameron & Miller,
  2015; Manski, 2003) instead of the numbered propositions of the
  accompanying manuscript.
* Package overview (`?silentema`) describing the workflow, the data format
  and the simulator.
* `DESCRIPTION`: method references with DOIs, `URL` and `BugReports` fields
  pointing to the GitHub repository and the documentation site; the title is
  shortened.
* `citation("silentema")` now points to the CRAN page of the package.
* American spelling throughout the documentation.
* Subsetting a `recoverability`, `silence_test` or `sensor_gap_test` object
  with `[` now returns a plain data frame; previously the class was kept and
  the print method showed an incomplete report. No estimate is affected.
* `fit_fiml()` stops with a message when persons have fewer than two prompts;
  `fatigue_check()` no longer fails when the prompt index has fewer than four
  distinct values.

# silentema 0.2.1

Metadata release; no change to any R or C++ code, so every result produced
with 0.2.0 is reproduced exactly (checked bit-for-bit on six simulation cells).

* `URL` points to the OSF project; the placeholder `BugReports` field is
  removed (contact the maintainer by email).
* `citation("silentema")` uses the revised title of the accompanying
  manuscript, "What skipped prompts hide: Detecting, diagnosing, and
  correcting informative nonresponse in ecological momentary assessment".

# silentema 0.2.0

Revision after peer review of the accompanying manuscript.

* `calibrate_delta(method = "postskip")` gains `burden = "fit"`: the simulated
  model includes a burden term (motif M3) calibrated to the observed response
  persistence, so that the post-skip calibration is valid under a declared
  M2 + M3 graph. The number of simulated data sets is now `n_sim` (default 20;
  formerly `B = 10`). The result reports the number of crossings, and no
  crossing is reported as such.
* `simulate_from_fit()` gains `days` and `kappa_R`, and stops with a clear
  message when the fitted dynamics are not stable.
* `recoverability()` distinguishes, under reactivity (M6), the recoverable
  assessment-conditioned kernel (Phi and Psi) from the dynamics-recovered
  person mean, which is biased; the observed person mean is judged by
  d-separation as for the other motifs.
* `fit_tilt()` returns every quantity from one undamped weighted fit at the
  converged iterate, reports `max_weight`, passes `min_pairs` through, and
  returns class `tilt_fit` also at `delta = 0`.
* `silence_test()` gains `poly` (polynomial degree in the lagged states) and
  `seed`; `fatigue_check()` gains `day` and a print method; `summary.pairs_fit()`
  uses a t(G - 1) reference; `coef()`, `vcov()`, `confint()`, `nobs()` methods
  for fits; `plot()` method for calibrations; `plot.dm_graph()` redrawn.
* `missingness_declaration()` fills its sections from the test, profile and
  calibration objects when they are supplied.
* Input validation with informative messages (motif codes, column names,
  integer prompt index, `delta` length, weights); functions that simulate or
  bootstrap restore the caller's random-number state.
* Examples for every exported function; American spelling throughout.

# silentema 0.1.0

* First release accompanying the manuscript "When silence is informative:
  Dynamic missingness graphs, recoverability, and sensitivity analysis for
  ecological momentary assessment" (Yu, 2026).
