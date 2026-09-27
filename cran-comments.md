## silentema 1.0.0: first submission to CRAN

This is a new package. It implements the methods of a manuscript under review (Yu, 2026, "What skipped prompts hide: Detecting, diagnosing, and correcting informative nonresponse in ecological momentary assessment"); the DESCRIPTION cites the published methods the package builds on (Mohan and Pearl, 2021; Cameron and Miller, 2015; Dhaene and Jochmans, 2015; Shumway and Stoffer, 1982; Manski, 2003) with DOIs.

## Test environments

* local: Ubuntu 24.04, R 4.3.3 (x86_64-pc-linux-gnu), gcc 13.3; `R CMD check --as-cran`
* win-builder: R-devel, R-release (https://win-builder.r-project.org) [UPDATE: date and result]
* macOS builder: R-release, arm64 (https://mac.r-project.org/macbuilder/) [UPDATE: date and result]
* R-hub v2 (GitHub Actions): linux, windows, macos, ubuntu-gcc12, clang-asan [UPDATE: date and result]
* GitHub Actions (r-lib/actions): macOS release, Windows release, Ubuntu devel/release/oldrel-1 [UPDATE: date and result]

## R CMD check results

0 errors | 0 warnings | 1 note (New submission) is the expected result on the external platforms. The local check reports 0 errors, 0 warnings and the following notes, all of which are environmental:

* `checking CRAN incoming feasibility ... NOTE: New submission`. This is the first release.

* `checking for future file timestamps ... NOTE: unable to verify current time`. The local check machine has no access to the time server.

* `checking compilation flags used ... NOTE: Compilation used the following non-portable flag(s): '-mno-omit-leaf-frame-pointer'`. The flag comes from the local R installation's site `Makeconf`, not from the package: `src/Makevars` sets only `PKG_LIBS = $(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS)`.

Should a platform report `checking installed package size ... NOTE` (the compiled Armadillo code is about 3.6 MB in `libs`), this is the usual size of a small RcppArmadillo package.

## Comments on the package

* The package contains C++ code (an EM algorithm for a state-space model) through Rcpp and RcppArmadillo; no C++ standard is specified in `Makevars` (the default of the R version is used), and the code compiles without warnings under gcc 13 and clang with `-Wall -Wextra -pedantic`.
* Every exported function has a documented `\value` section and runnable examples; there is no `\dontrun{}` or `\donttest{}`. Examples run in about 7 seconds in total (none above 1.2 s), the tests in about 25 seconds and the five vignettes in about 60 seconds on the local machine.
* Functions that simulate or bootstrap accept an optional `seed` argument. When it is supplied, the caller's random-number state (`.Random.seed`) is saved and restored on exit, so the functions never change the user's random-number stream; when it is not supplied (the default), the current stream is used and advanced as usual. The package never sets a seed on its own.
* Nothing is written to the user's file system except by `missingness_declaration(file = )` when the user gives a path; the examples and tests write only to `tempfile()`.
* Spelling: the words that the DESCRIPTION spell check flags (Dhaene, Jochmans, Manski, Mohan, Shumway, Stoffer, EMA, estimands, nonresponse) are author names of the cited references and standard terms of the missing-data and experience-sampling literatures; `inst/WORDLIST` lists them.

## Downstream dependencies

There are none; this is a new package.
