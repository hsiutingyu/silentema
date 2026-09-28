# Internal: weighted within (person fixed-effects) regression of a Y matrix on an X matrix
# with cluster-robust covariance of the stacked coefficients (outcome by outcome).
within_reg <- function(Y, X, g, w = NULL) {
  Y <- as.matrix(Y); X <- as.matrix(X)
  n <- nrow(Y); q <- ncol(Y); k <- ncol(X)
  if (is.null(w)) w <- rep(1, n)
  gf <- factor(g); G <- nlevels(gf); gi <- as.integer(gf)
  sw <- as.numeric(tapply(w, gf, sum))
  Ybar <- matrix(apply(Y, 2, function(v) tapply(w * v, gf, sum)) / sw, ncol = q)
  Xbar <- matrix(apply(X, 2, function(v) tapply(w * v, gf, sum)) / sw, ncol = k)
  Yd <- Y - Ybar[gi, , drop = FALSE]; Xd <- X - Xbar[gi, , drop = FALSE]
  Sxx <- crossprod(Xd * w, Xd)
  B <- solve(Sxx, crossprod(Xd * w, Yd))          # k x q
  res <- Yd - Xd %*% B
  Sxx_inv <- solve(Sxx)
  scores <- matrix(0, G, k * q)
  for (j in seq_len(q)) {
    sj <- Xd * (w * res[, j])
    scores[, (j - 1) * k + seq_len(k)] <- apply(sj, 2, function(v) tapply(v, gf, sum))
  }
  meat <- crossprod(scores)
  bread <- kronecker(diag(q), Sxx_inv)
  adj <- if (G > 1) G / (G - 1) else 1
  V <- adj * bread %*% meat %*% bread
  df_res <- max(n - G - k, 1)
  V_iid <- kronecker(crossprod(res * w, res) / df_res, Sxx_inv)
  list(coef = B, vcov = V, vcov_iid = V_iid, residuals = res, G = G, n = n, Xd = Xd, Yd = Yd, gf = gf)
}

# Internal: block bootstrap over person-day blocks (resampled with replacement) for a
# single-person (or few-person) regression; returns the covariance of the coefficients.
block_boot_coef <- function(Y, X, g, block, B = 500, seed = NULL) {
  with_seed(seed, {
    bl <- factor(block); ub <- levels(bl); nb <- length(ub)
    idx_by_block <- split(seq_len(nrow(as.matrix(Y))), bl)
    out <- matrix(NA_real_, B, ncol(as.matrix(X)) * ncol(as.matrix(Y)))
    for (b in seq_len(B)) {
      take <- sample(ub, nb, replace = TRUE)
      idx <- unlist(idx_by_block[take], use.names = FALSE)
      fr <- tryCatch(within_reg(as.matrix(Y)[idx, , drop = FALSE], as.matrix(X)[idx, , drop = FALSE], g[idx]), error = function(e) NULL)
      if (!is.null(fr)) out[b, ] <- as.vector(fr$coef)
    }
    stats::cov(out, use = "complete.obs")
  })
}

#' The silence test: does the state after a skipped prompt differ?
#'
#' Regresses the state at prompt \eqn{t+1} on the state at \eqn{t-1} and the
#' response indicator at \eqn{t}, within person, using all triples in which
#' prompts \eqn{t-1} and \eqn{t+1} were answered (prompt \eqn{t} may or may
#' not have been). Under every combination of the recoverable motifs
#' (M0, M1, M3, M4) that does not contain both M1 and M3, \eqn{R_t} is
#' conditionally independent of \eqn{X_{t+1}} given \eqn{X_{t-1}} and the
#' person (Yu, 2026), so the population coefficient of \eqn{R_t} is zero
#' when the conditional mean of \eqn{X_{t+1}} given \eqn{X_{t-1}} is linear
#' (exactly under M0, M3 and M4; under M1 the selection on \eqn{X_t} through
#' \eqn{R_{t+1}} makes it slightly nonlinear, and \code{poly = 2} adds squared
#' lagged states to absorb the curvature). Under self-censoring (M2), latent
#' context (M5) or reactivity (M6) the coefficient is not zero. When burden
#' (M3) is combined with a state-dependent motif (M1 + M3, M2 + M3), the
#' conditioning on \eqn{R_{t+1} = 1} opens a collider at \eqn{R_{t+1}}, whose
#' parents are \eqn{R_t} (burden) and either \eqn{X_t} (under M1, and
#' \eqn{X_t} drives \eqn{X_{t+1}}) or \eqn{X_{t+1}} itself (under M2): under
#' M1 + M3 the test rejects in large samples although the kernel is
#' recoverable (with a coefficient of the opposite sign to the self-censoring
#' signature), and under M2 + M3 the collider works against the self-censoring
#' signal and the test loses power. The sign of the
#' coefficient on the self-censoring variable is the sign of the tilt: a
#' negative value means that the states hidden by skips were higher than the
#' states that were reported.
#'
#' @inheritParams make_pairs
#' @param se \code{"cluster"} (cluster-robust by person with \eqn{t(G-1)} and
#'   \eqn{F(q, G-1)} reference distributions; the default when there are at
#'   least 10 persons) or \code{"dayblock"} (block bootstrap that resamples
#'   person-day blocks with replacement; the default with fewer than 10
#'   persons; requires \code{day}).
#' @param B Number of bootstrap resamples for \code{se = "dayblock"}.
#' @param poly Degree of the polynomial in the lagged states (1 = linear).
#' @param seed Optional seed for the bootstrap (the caller's random-number
#'   state is restored).
#' @return An object of class \code{c("silence_test", "data.frame")}: a data
#'   frame with one row per state variable and the columns \code{variable},
#'   \code{coef_R} (coefficient of \eqn{R_t}), \code{se}, \code{z} (the ratio
#'   of the coefficient to its standard error, referred to \eqn{t(G - 1)} under
#'   \code{se = "cluster"}) and \code{p}, and the attributes \code{joint} (a named vector with the Wald
#'   statistic over all variables, its degrees of freedom and its p value),
#'   \code{n_triples}, \code{n_skipped} (triples with a skipped middle
#'   prompt) and \code{se} (the method used). It has a \code{print} method.
#' @references Cameron, A. C., & Miller, D. L. (2015). A practitioner's guide
#'   to cluster-robust inference. \emph{Journal of Human Resources, 50},
#'   317-372. \doi{10.3368/jhr.50.2.317}
#'
#'   Yu, H.-T. (2026). What skipped prompts hide: Detecting, diagnosing, and
#'   correcting informative nonresponse in ecological momentary assessment.
#'   Manuscript under review.
#' @seealso \code{\link{sensor_gap_test}}, \code{\link{fatigue_check}},
#'   \code{\link{recoverability}} (which reports whether the null is expected
#'   under a declared graph)
#' @examples
#' sim <- simulate_ema(N = 40, n_prompts = 30, motifs = "M2", delta = -1, seed = 1)
#' st <- silence_test(sim$data, sim$vars)
#' st
#' attr(st, "joint")
#' sim0 <- simulate_ema(N = 40, n_prompts = 30, motifs = "M1", seed = 1)
#' silence_test(sim0$data, sim0$vars)
#' # quadratic terms in the lagged states, and a day structure
#' simd <- simulate_ema(N = 40, n_prompts = 30, motifs = "M2", delta = -1, days = 5, seed = 1)
#' silence_test(simd$data, simd$vars, day = "day", poly = 2)
#' @export
silence_test <- function(data, vars, id = "id", time = "time", day = NULL, R = "R",
                         se = NULL, B = 500, poly = 1L, seed = NULL) {
  check_columns(data, c(vars, id, time, day))
  check_time(data, time)
  Rv <- response_indicator(data, vars, R)
  if (!is.null(se)) se <- match.arg(se, c("cluster", "dayblock"))
  ord <- order(data[[id]], data[[time]]); data <- data[ord, , drop = FALSE]; Rv <- Rv[ord]
  rid <- data[[id]]; rt <- data[[time]]
  n <- nrow(data)
  i <- seq_len(n - 2L)
  ok <- rid[i] == rid[i + 2L] & rt[i + 2L] == rt[i] + 2L & Rv[i] == 1L & Rv[i + 2L] == 1L
  if (!is.null(day)) ok <- ok & data[[day]][i] == data[[day]][i + 2L]
  i <- i[ok]
  if (length(i) < 10)
    stop(sprintf("Too few triples for the silence test: %d of %d prompts start a triple with both outer prompts answered (at least 10 are needed); check that '%s' is a consecutive integer prompt index within person%s",
                 length(i), n, time, if (!is.null(day)) " and day" else ""))
  Xm <- as.matrix(data[, vars, drop = FALSE])
  Y <- Xm[i + 2L, , drop = FALSE]; Xlag <- Xm[i, , drop = FALSE]; Rt <- Rv[i + 1L]
  if (poly >= 2) for (dg in 2:poly) Xlag <- cbind(Xlag, Xlag[, seq_along(vars), drop = FALSE]^dg)
  X <- cbind(Xlag, R_t = Rt)
  g <- rid[i + 2L]
  G <- length(unique(g))
  if (is.null(se)) se <- if (G >= 10) "cluster" else "dayblock"
  fr <- within_reg(Y, X, g)
  k <- ncol(X); q <- ncol(Y)
  pos <- (seq_len(q) - 1) * k + k                # positions of the R_t coefficient
  est <- fr$coef[k, ]
  if (se == "cluster") {
    V <- fr$vcov[pos, pos, drop = FALSE]
  } else {
    if (is.null(day)) stop("the day-block bootstrap (the default with fewer than 10 persons) needs a 'day' column; give 'day', or set se = \"cluster\"")
    Vall <- block_boot_coef(Y, X, g, block = paste(g, data[[day]][i + 2L]), B = B, seed = seed)
    V <- Vall[pos, pos, drop = FALSE]
  }
  s <- sqrt(diag(V))
  z <- est / s
  # small-sample reference distributions for cluster-robust inference: t(G - 1) and F(q, G - 1)
  # (Cameron & Miller, 2015); the normal and chi-square are used for the day-block bootstrap
  pval <- if (se == "cluster") 2 * stats::pt(-abs(z), df = max(G - 1, 1)) else 2 * stats::pnorm(-abs(z))
  out <- data.frame(variable = vars, coef_R = est, se = s, z = z, p = pval, row.names = NULL)
  joint <- tryCatch(as.numeric(t(est) %*% solve(V, est)), error = function(e) NA_real_)
  joint_p <- if (is.na(joint)) NA_real_ else if (se == "cluster") stats::pf(joint / q, q, max(G - 1, 1), lower.tail = FALSE) else stats::pchisq(joint, q, lower.tail = FALSE)
  attr(out, "joint") <- c(chisq = joint, df = q, p = joint_p)
  attr(out, "n_triples") <- length(i); attr(out, "n_skipped") <- sum(Rt == 0L)
  attr(out, "se") <- se
  class(out) <- c("silence_test", "data.frame")
  out
}

#' @export
`[.silence_test` <- function(x, ...) {
  out <- NextMethod()
  if (is.data.frame(out)) class(out) <- "data.frame"
  out
}

#' @export
`[.sensor_gap_test` <- function(x, ...) {
  out <- NextMethod()
  if (is.data.frame(out)) class(out) <- "data.frame"
  out
}

#' @export
print.silence_test <- function(x, digits = 3, ...) {
  cat("Silence test: X_{t+1} ~ X_{t-1} + R_t (within person),", attr(x, "n_triples"), "triples,",
      attr(x, "n_skipped"), "with a skipped middle prompt\n")
  print(data.frame(variable = x$variable, coef_R = round(x$coef_R, digits), se = round(x$se, digits),
                   z = round(x$z, 2), p = signif(x$p, 2)), row.names = FALSE)
  j <- attr(x, "joint")
  cat(sprintf("Joint Wald test: statistic = %.2f on %d df, p = %.3g   [SE: %s; %s reference]\n", j["chisq"], j["df"], j["p"], attr(x, "se"),
              if (attr(x, "se") == "cluster") "t / F with G - 1 df" else "normal / chi-square"))
  cat("Interpretation: coefficients near zero are expected under the recoverable motifs;",
      "a negative coefficient on a variable means the states hidden by skips were higher.\n")
  invisible(x)
}

#' The sensor-gap test: does an always-observed sensor differ at skipped prompts?
#'
#' Regresses a passive-sensor channel at prompt \eqn{t} on the state at
#' \eqn{t-1} and the response indicator at \eqn{t}, within person, over all
#' prompts whose predecessor was answered. Under the recoverable motifs the
#' coefficient of \eqn{R_t} is zero, including under M1 + M3 (no collider is
#' conditioned on, because the sensor is always observed). Under
#' self-censoring the gap identifies the direction of the tilt and, through
#' \code{\link{calibrate_delta}}, its magnitude.
#'
#' @inheritParams silence_test
#' @param sensor Name of the sensor column.
#' @return A one-row data frame of class \code{c("sensor_gap_test",
#'   "data.frame")} with the columns \code{coef_R} (coefficient of
#'   \eqn{R_t}), \code{se}, \code{z} and \code{p}, and the attributes
#'   \code{loading} (the regression coefficients of the sensor on the state
#'   variables, estimated within person from answered prompts), \code{n} (the
#'   number of prompts used) and \code{se} (the method used). It has a
#'   \code{print} method.
#' @examples
#' sim <- simulate_ema(N = 40, n_prompts = 30, motifs = "M2", delta = -1, sensor_cor = 0.6, seed = 1)
#' sg <- sensor_gap_test(sim$data, sim$vars, sensor = "S")
#' sg
#' attr(sg, "loading")
#' # under a recoverable mechanism the gap is null
#' sim1 <- simulate_ema(N = 40, n_prompts = 30, motifs = "M1", sensor_cor = 0.6, seed = 1)
#' sensor_gap_test(sim1$data, sim1$vars, sensor = "S")
#' @export
sensor_gap_test <- function(data, vars, sensor = "S", id = "id", time = "time", day = NULL, R = "R",
                            se = NULL, B = 500, seed = NULL) {
  check_columns(data, sensor)
  if (!is.null(se)) se <- match.arg(se, c("cluster", "dayblock"))
  pr <- make_pairs(data, vars, id, time, day, R)
  lo <- pr$lag_obs
  S <- data[[sensor]][lo$row_t]
  X <- cbind(lo$x, R_t = lo$R)
  g <- lo$id
  G <- length(unique(g))
  if (is.null(se)) se <- if (G >= 10) "cluster" else "dayblock"
  fr <- within_reg(matrix(S, ncol = 1), X, g)
  k <- ncol(X)
  est <- fr$coef[k, 1]
  V <- if (se == "cluster") fr$vcov[k, k] else {
    if (is.null(day)) stop("the day-block bootstrap (the default with fewer than 10 persons) needs a 'day' column; give 'day', or set se = \"cluster\"")
    block_boot_coef(matrix(S, ncol = 1), X, g, block = paste(g, data[[day]][lo$row_t]), B = B, seed = seed)[k, k]
  }
  s <- sqrt(V); z <- est / s
  pval <- if (se == "cluster") 2 * stats::pt(-abs(z), df = max(G - 1, 1)) else 2 * stats::pnorm(-abs(z))
  # sensor loading from answered prompts (S ~ X_t within person)
  ok <- which(pr$Rv == 1L)
  fl <- within_reg(matrix(data[[sensor]][ok], ncol = 1), as.matrix(data[ok, vars, drop = FALSE]), data[[id]][ok])
  out <- data.frame(coef_R = est, se = s, z = z, p = pval)
  attr(out, "loading") <- stats::setNames(as.vector(fl$coef), vars)
  attr(out, "n") <- length(S); attr(out, "se") <- se
  class(out) <- c("sensor_gap_test", "data.frame")
  out
}

#' @export
print.sensor_gap_test <- function(x, digits = 3, ...) {
  cat("Sensor-gap test: S_t ~ X_{t-1} + R_t (within person),", attr(x, "n"), "prompts with answered predecessor\n")
  cat(sprintf("  coefficient of R_t: %.3f (SE %.3f), z = %.2f, p = %.3g   [SE: %s]\n", x$coef_R, x$se, x$z, x$p, attr(x, "se")))
  cat("  sensor loading on states (answered prompts):", round(attr(x, "loading"), digits), "\n")
  invisible(x)
}

#' Descriptive check for burden or fatigue in the response sequence
#'
#' Reports the response rate after an answered and after a skipped prompt
#' (pooled over consecutive prompts, within day when \code{day} is given, and
#' as the mean over persons of the person-specific rates) and the response
#' rate by study quarter. Serial dependence in \eqn{R} is expected under M3,
#' but some serial dependence also arises under M1, M2 and M4 through the
#' autocorrelated states and the person propensities, so this is a
#' descriptive check, not a test of M3 against the other motifs.
#' @inheritParams make_pairs
#' @param R Name of the response-indicator column (0/1, no \code{NA});
#'   required, because the check has no state columns from which to infer it.
#' @return An object of class \code{"fatigue_check"}: a list with
#'   \code{after_answered} and \code{after_skipped} (pooled response rates
#'   after an answered and after a skipped prompt), \code{difference}
#'   (their difference, pooled), \code{difference_person} (mean of the
#'   person-specific differences), \code{by_quarter} (response rate by study
#'   quarter) and \code{n_transitions}. It has a \code{print} method.
#' @examples
#' sim <- simulate_ema(N = 40, n_prompts = 30, motifs = "M3", seed = 1)
#' fatigue_check(sim$data)
#' # compare with a mechanism without burden
#' fatigue_check(simulate_ema(N = 40, n_prompts = 30, motifs = "M0", seed = 1)$data)
#' @export
fatigue_check <- function(data, id = "id", time = "time", day = NULL, R = "R") {
  check_columns(data, c(id, time, day, R))
  check_time(data, time)
  Rv <- response_indicator(data, character(0), R)
  ord <- order(data[[id]], data[[time]]); data <- data[ord, , drop = FALSE]; Rv <- Rv[ord]
  rid <- data[[id]]; rt <- data[[time]]
  n <- nrow(data); i <- seq_len(n - 1L)
  ok <- rid[i] == rid[i + 1L] & rt[i + 1L] == rt[i] + 1L
  if (!is.null(day)) ok <- ok & data[[day]][i] == data[[day]][i + 1L]
  i <- i[ok]
  d <- data.frame(id = rid[i + 1L], prev = Rv[i], cur = Rv[i + 1L])
  by_person <- lapply(split(d, d$id), function(s) c(a = mean(s$cur[s$prev == 1L]), s = mean(s$cur[s$prev == 0L])))
  m <- do.call(rbind, by_person)
  br <- unique(stats::quantile(rt, c(0, .25, .5, .75, 1)))
  q <- if (length(br) >= 2L) cut(rt, breaks = br, include.lowest = TRUE, labels = if (length(br) == 5L) paste0("Q", 1:4) else NULL) else factor(rep("Q1", length(rt)))
  structure(list(after_answered = mean(d$cur[d$prev == 1L]), after_skipped = mean(d$cur[d$prev == 0L]),
                 difference = mean(d$cur[d$prev == 1L]) - mean(d$cur[d$prev == 0L]),
                 difference_person = mean(m[, "a"] - m[, "s"], na.rm = TRUE),
                 by_quarter = tapply(Rv, q, mean), n_transitions = nrow(d)), class = "fatigue_check")
}

#' @export
print.fatigue_check <- function(x, digits = 3, ...) {
  cat("Fatigue check over", x$n_transitions, "consecutive prompts\n")
  cat(sprintf("  response rate after an answered prompt: %.3f; after a skipped prompt: %.3f; difference %.3f (person-mean difference %.3f)\n",
              x$after_answered, x$after_skipped, x$difference, x$difference_person))
  cat("  response rate by study quarter:", round(x$by_quarter, digits), "\n")
  invisible(x)
}
