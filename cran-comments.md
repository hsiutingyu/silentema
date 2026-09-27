## silentema 1.0.0: first submission to CRAN

This is a new package. It implements the methods of a manuscript under review (Yu, 2026, "What skipped prompts hide: Detecting, diagnosing, and correcting informative nonresponse in ecological momentary assessment"); the DESCRIPTION cites the published methods the package builds on (Mohan and Pearl, 2021; Cameron and Miller, 2015; Dhaene and Jochmans, 2015; Shumway and Stoffer, 1982; Manski, 2003) with DOIs.

## Test environments

* local: Ubuntu 24.04, R 4.3.3 (x86_64-pc-linux-gnu), gcc 13.3; `R CMD check --as-cran`
* win-builder: R-devel, R-release (https://win-builder.r-project.org) [UPDATE: date and result]
* macOS builder: R-release, arm64 (https://mac.r-project.org/macbuilder/) [UPDATE: date and result]
* R-hub v2 (GitHub Actions): ubuntu-release, windows, macos, ubuntu-gcc12 [UPDATE: date and result]
* GitHub Actions (r-lib/actions): macOS release, Windows release, Ubuntu devel/release/oldrel-1 [UPDATE: date and result]

## R CMD check results

0 errors | 0 warnings | 2 notes (local)

* `checking CRAN incoming feasibility ... NOTE: New submission`. This is the first release.

* `checking for future file timestamps ... NOTE: unable to verify current time`. The local check machine has no access to the time server; the note does not occur on win-builder or R-hub.

* `checking compilation flags used ... NOTE: Compilation used the following non-portable flag(s): '-mno-omit-leaf-frame-pointer'`. The flag comes from the local R installation's site `Makeconf`, not from the package: `src/Makevars` sets only `PKG_LIBS = $(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS)`. The note does not occur on win-builder or R-hub.

## Comments on the package

* The package contains C++ code (an EM algorithm for a state-space model) through Rcpp and RcppArmadillo; no C++ standard is specified in `Makevars` (the default of the R version is used).
* Every exported function has a documented `\value` section and runnable examples; the one `\donttest{}` block (`calibrate_delta`, post-skip calibration) is only slower, not unrunnable, and passes `--run-donttest`.
* Functions that simulate or bootstrap accept an optional `seed` argument. When it is supplied, the caller's random-number state (`.Random.seed`) is saved and restored on exit, so the functions never change the user's random-number stream; when it is not supplied, the current stream is used and advanced as usual.
* Nothing is written to the user's file system except by `missingness_declaration(file = )` when the user gives a path; the examples and tests write only to `tempfile()`.
* Examples run in about 10 seconds in total, the tests in about 20 seconds and the five vignettes in about 60 seconds on the local machine.
* Spelling: `inst/WORDLIST` lists the domain terms (EMA, nonresponse, recoverability, estimand, motif codes M0-M6, author names).

## Downstream dependencies

There are none; this is a new package.
