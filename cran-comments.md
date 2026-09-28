## Submission

This is the first submission of silentema to CRAN.

## Test environments

* Local: Ubuntu 24.04, R 4.3.3, gcc 13.3 (R CMD check --as-cran)
* win-builder: R-devel and R-release
* macOS builder: R-release (arm64)
* GitHub Actions: macOS, Windows, Ubuntu (R-devel, R-release, R-oldrel-1)

## R CMD check results

0 errors | 0 warnings | 1 note

* checking CRAN incoming feasibility ... NOTE
  New submission
  Possibly misspelled words in DESCRIPTION: Dhaene, EMA, Jochmans, Manski,
  Mohan, Shumway, Stoffer, Yu, estimands, nonresponse

This is a first submission. None of the flagged words is a misspelling:
Dhaene, Jochmans, Manski, Mohan, Shumway, Stoffer and Yu are the surnames of
the authors of the cited references; EMA is defined in the Description
(ecological momentary assessment); "estimands" and "nonresponse" are standard
terms of the missing-data literature.

## Notes on the package

* The Description cites the published methods the package implements, in the
  form authors (year) <doi:...>.
* All exported functions have a \value section and runnable examples; there
  is no \dontrun or \donttest. Examples run in about 7 s in total.
* C++ code (a Kalman-filter EM algorithm) uses Rcpp and RcppArmadillo;
  Makevars sets only the LAPACK/BLAS libraries.
* Functions with a seed argument restore the user's .Random.seed on exit;
  with the default seed = NULL the current stream is used. The package never
  writes to the user's file space except when a file path is supplied to
  missingness_declaration().

## Downstream dependencies

There are currently no downstream dependencies for this package.
