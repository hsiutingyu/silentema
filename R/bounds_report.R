#' Worst-case (support) bounds for person means under arbitrary nonresponse
#'
#' Worst-case bounds in the spirit of Manski (2003): every skipped prompt
#' could have taken any value on the response scale \code{[L, U]}. The width,
#' \eqn{(1 - \pi_i)(U - L)}, depends only on the response rate and does not
#' shrink with autocorrelation, in contrast to the tilt-bounded identified
#' sets of \code{\link{break_even}}.
#'
#' @inheritParams make_pairs
#' @param id Name of the person column.
#' @param L,U Lower and upper limits of the response scale (scalars or one
#'   value per variable).
#' @return A data frame with one row per person and variable and the columns
#'   \code{id}, \code{variable}, \code{response_rate}, \code{mean_observed}
#'   (\code{NA} for a person without any answered prompt, whose bounds are then
#'   the whole scale), \code{lower} and \code{upper}.
#' @references Manski, C. F. (2003). \emph{Partial identification of
#'   probability distributions}. Springer. \doi{10.1007/b97478}
#' @examples
#' sim <- simulate_ema(N = 5, n_prompts = 20, motifs = "M2", seed = 1)
#' b <- bounds_support(sim$data, sim$vars, L = -3, U = 8)
#' b[b$variable == "NegA", ]
#' # the width is (1 - response rate) x (U - L)
#' with(b, all.equal(upper - lower, (1 - response_rate) * 11))
#' @export
bounds_support <- function(data, vars, L, U, id = "id", R = "R") {
  p <- length(vars)
  check_columns(data, c(vars, id))
  if (length(L) == 1L) L <- rep(L, p); if (length(U) == 1L) U <- rep(U, p)
  if (length(L) != p || length(U) != p) stop("'L' and 'U' must be scalars or have one value per variable")
  Rv <- response_indicator(data, vars, R)
  ids <- unique(data[[id]])
  out <- list()
  for (i in ids) {
    sel <- data[[id]] == i
    pi_i <- mean(Rv[sel])
    for (j in seq_len(p)) {
      xo <- data[[vars[j]]][sel & Rv == 1L]
      m <- if (length(xo)) mean(xo) else NA_real_       # no answered prompt: the bounds are the whole scale
      out[[length(out) + 1L]] <- data.frame(id = i, variable = vars[j], response_rate = pi_i, mean_observed = m,
                                            lower = if (pi_i > 0) pi_i * m + (1 - pi_i) * L[j] else L[j],
                                            upper = if (pi_i > 0) pi_i * m + (1 - pi_i) * U[j] else U[j])
    }
  }
  do.call(rbind, out)
}

#' Missingness declaration for preregistrations and reports
#'
#' Writes a Markdown declaration that records the declared dm-graph, the
#' recoverability verdicts, the tests that were run, the sensitivity
#' interval and its source, and the reported band, in the order that Yu
#' (2026) recommends for preregistrations and reports. Sections for which the
#' corresponding object is supplied are filled from it; the others are left as
#' prompts to complete.
#'
#' @param g A \code{dm_graph}.
#' @param rec Optional \code{\link{recoverability}} report for \code{g}.
#' @param silence Optional \code{\link{silence_test}} result.
#' @param sensor_gap Optional \code{\link{sensor_gap_test}} result.
#' @param fatigue Optional \code{\link{fatigue_check}} result.
#' @param profile Optional \code{\link{tilt_profile}}.
#' @param calibration Optional \code{\link{calibrate_delta}} result (or a list
#'   of them).
#' @param plausible Optional numeric vector of length 2: the plausible interval
#'   used for the band.
#' @param file Optional path; if \code{NULL} the declaration is returned as a
#'   character vector.
#' @return Invisibly, a character vector with one element per line of the
#'   Markdown declaration; when \code{file} is given the lines are also
#'   written to that file.
#' @examples
#' g <- dm_graph(c("M2", "M4"), sensor = TRUE)
#' sim <- simulate_ema(N = 40, n_prompts = 30, motifs = c("M2", "M4"), sensor_cor = 0.6, seed = 1)
#' st <- silence_test(sim$data, sim$vars)
#' prof <- tilt_profile(sim$data, sim$vars, delta_grid = c(-1, -0.5, 0), propensity = "person")
#' cal <- calibrate_delta(prof, sim$data, method = "sensor")
#' cat(missingness_declaration(g, silence = st, profile = prof, calibration = cal,
#'                             plausible = c(-1, 0)), sep = "\n")
#' # a bare template for a preregistration, written to a file
#' tf <- tempfile(fileext = ".md")
#' missingness_declaration(dm_graph(c("M1", "M3")), file = tf)
#' readLines(tf)[1:12]
#' @export
missingness_declaration <- function(g, rec = NULL, silence = NULL, sensor_gap = NULL, fatigue = NULL,
                                    profile = NULL, calibration = NULL, plausible = NULL, file = NULL) {
  if (is.null(rec)) rec <- recoverability(g)
  f3 <- function(v) formatC(v, digits = 3, format = "f")
  sec_tests <- c("## 4. Tests of informativeness")
  sec_tests <- c(sec_tests, if (is.null(silence)) "- Silence test (coefficient of R_t, SE, p, per variable): [to complete]" else
    c(sprintf("- Silence test (%d triples, %d with a skipped middle prompt; SE: %s):", attr(silence, "n_triples"), attr(silence, "n_skipped"), attr(silence, "se")),
      sprintf("    - %s: coefficient of R_t = %s (SE %s), p = %s", silence$variable, f3(silence$coef_R), f3(silence$se), formatC(silence$p, digits = 2, format = "g")),
      sprintf("    - joint test: statistic %s on %d df, p = %s", f3(attr(silence, "joint")["chisq"]), as.integer(attr(silence, "joint")["df"]), formatC(attr(silence, "joint")["p"], digits = 2, format = "g"))))
  sec_tests <- c(sec_tests, if (is.null(sensor_gap)) "- Sensor-gap test (if a sensor is available): [to complete or state that no sensor exists]" else
    sprintf("- Sensor-gap test: coefficient of R_t = %s (SE %s), p = %s", f3(sensor_gap$coef_R), f3(sensor_gap$se), formatC(sensor_gap$p, digits = 2, format = "g")))
  sec_tests <- c(sec_tests, if (is.null(fatigue)) "- Fatigue check (response rate after answered vs skipped prompt; by study quarter): [to complete]" else
    sprintf("- Fatigue check: response rate %s after an answered and %s after a skipped prompt; by quarter %s",
            f3(fatigue$after_answered), f3(fatigue$after_skipped), paste(f3(fatigue$by_quarter), collapse = ", ")))
  sec_sens <- c("## 5. Sensitivity analysis")
  if (is.null(profile)) sec_sens <- c(sec_sens, "- Sensitivity parameter, its scale, and the grid: [to complete]") else {
    sec_sens <- c(sec_sens, sprintf("- Sensitivity parameter: delta, probit units per unit of %s; grid {%s}; response-model intercept: %s",
                                    profile$vars[profile$which], paste(profile$delta_grid, collapse = ", "), profile$propensity),
                  sprintf("- Effective number of complete pairs along the grid: %s (of %d)", paste(round(profile$ess), collapse = ", "), profile$n_pairs))
  }
  cals <- if (is.null(calibration)) NULL else if (inherits(calibration, "delta_calibration")) list(calibration) else calibration
  if (is.null(cals)) sec_sens <- c(sec_sens, "- Plausible interval and its source (sensor calibration / probe calibration / post-skip calibration / elicitation): [to complete]") else
    for (cl in cals) sec_sens <- c(sec_sens, sprintf("- Calibration by the %s method%s: delta = %s, %d%% interval [%s, %s]%s",
                                                     cl$method, if (identical(cl$burden, "fit")) " (burden term fitted)" else "",
                                                     if (is.na(cl$delta_hat)) "not identified on the grid" else f3(cl$delta_hat), round(100 * cl$level), f3(cl$interval[1]), f3(cl$interval[2]),
                                                     if (!is.null(cl$boot_interval)) sprintf("; bootstrap percentile interval [%s, %s]", f3(cl$boot_interval[1]), f3(cl$boot_interval[2])) else ""))
  if (!is.null(profile) && !is.null(plausible)) {
    be <- break_even(profile, delta_max = max(abs(plausible)))
    foc <- be[grepl(paste0("<-", profile$vars[profile$which], "$"), be$coef), ]
    sec_sens <- c(sec_sens, sprintf("- Plausible interval used for the band: [%s, %s]", f3(plausible[1]), f3(plausible[2])),
                  "- Identified set and band over the plausible interval (effects of the self-censoring variable):",
                  sprintf("    - %s: estimate at delta = 0: %s; set [%s, %s]; band [%s, %s]; sign change at %s; significance change at %s",
                          foc$coef, f3(foc$estimate_0), f3(foc$set_lower), f3(foc$set_upper), f3(foc$band_lower), f3(foc$band_upper),
                          ifelse(is.na(foc$sign_flip_delta), "none on the grid", f3(foc$sign_flip_delta)),
                          ifelse(is.na(foc$significance_flip_delta), "none on the grid", f3(foc$significance_flip_delta))))
  } else sec_sens <- c(sec_sens, "- Break-even values and the band over the plausible interval: [to complete]")
  lines <- c(
    "# Missingness declaration (silentema)",
    "",
    "## 1. Design facts",
    "- Prompts per day, days, scheduling, and the definition of a scheduled prompt: [to complete]",
    "- Response rate overall and by person (median, range): [to complete]",
    "- Whether prompt-level (whole prompt) or item-level missingness occurs: [to complete]",
    "",
    "## 2. Declared dynamic missingness graph",
    paste0("- Motifs: ", paste(g$motifs, collapse = " + ")),
    paste0("- Context: ", g$context, if (g$context_persistent) " (persistent)" else ""),
    paste0("- Passive sensor available: ", g$sensor, "; randomized probes: ", g$probe),
    "- Substantive justification for each declared edge (cite compliance evidence or pilot data): [to complete]",
    "",
    "## 3. Recoverability verdicts",
    paste0("- ", rec$estimand, ": ", ifelse(is.na(rec$recoverable), "", ifelse(rec$recoverable, "recoverable; ", "NOT recoverable; ")), rec$estimator),
    "",
    sec_tests,
    "",
    sec_sens,
    "",
    "## 6. What is reported in the paper",
    "- Estimates under the declared graph, the sensitivity band, and this declaration."
  )
  if (!is.null(file)) writeLines(lines, file)
  invisible(lines)
}
