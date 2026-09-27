# Internal: argument checks shared by the estimators
check_columns <- function(data, cols, what = "data") {
  if (!is.data.frame(data)) stop("'", what, "' must be a data frame")
  miss <- setdiff(cols, names(data))
  if (length(miss)) stop("column(s) not found in '", what, "': ", paste(miss, collapse = ", "))
  invisible(TRUE)
}

# Internal: the response indicator as a validated integer 0/1 vector. When the named column exists it
# must be 0/1 without NA and the state columns must be complete at answered prompts; when the default
# name "R" is absent, a prompt counts as answered when all state columns are non-missing.
response_indicator <- function(data, vars, R = "R", what = "data") {
  if (!is.null(R) && R %in% names(data)) {
    Rv <- data[[R]]
    if (anyNA(Rv) || !all(Rv %in% c(0, 1))) stop("'", R, "' must contain only 0 and 1, without NA")
    Rv <- as.integer(Rv)
    if (length(vars) && !all(stats::complete.cases(data[Rv == 1L, vars, drop = FALSE])))
      stop("state columns contain NA at prompts with ", R, " = 1 (item-level missingness): set ", R,
           " to 0 at those prompts, or drop the column so that a prompt counts as answered when all state columns are non-missing")
    return(Rv)
  }
  if (!is.null(R) && !identical(R, "R")) stop("column '", R, "' not found in '", what, "'")
  if (!length(vars)) stop("a response-indicator column '", R, "' is needed")
  as.integer(stats::complete.cases(data[, vars, drop = FALSE]))
}

# Internal: the prompt index must be an integer index
check_time <- function(data, time) {
  if (!is.numeric(data[[time]]) || anyNA(data[[time]]) || any(abs(data[[time]] - round(data[[time]])) > 1e-8))
    stop("'", time, "' must be an integer prompt index (1, 2, 3, ... within person, without NA); see ?make_pairs")
  invisible(TRUE)
}

#' Build adjacent prompt pairs from long experience-sampling data
#'
#' A complete pair is two adjacent scheduled prompts of one person (within a
#' day, when a day column is given) that were both answered. Complete pairs are
#' the building block of every estimator in the package: \code{\link{fit_pairs}},
#' \code{\link{fit_tilt}} and \code{\link{fit_ipw}} regress the states at the
#' second prompt on the states at the first. This function forms the pairs once
#' and also collects the prompts whose predecessor was answered, which the
#' response models of the weighted estimators use.
#'
#' @param data Long data frame with one row per scheduled prompt.
#' @param vars Names of the state columns.
#' @param id,time Names of the person and prompt-index columns.
#' @param day Optional name of a day column; pairs are formed only within a
#'   day (the overnight gap is not a lag-1 transition).
#' @param R Name of the response-indicator column (0/1, no \code{NA}). If the
#'   default name is not among the columns of \code{data}, a prompt counts as
#'   answered when all \code{vars} are non-missing; any other name must exist.
#' @details The \code{time} column must be an integer prompt index that
#'   increases by one from one scheduled prompt to the next within a person
#'   (and within a day when \code{day} is given); two rows whose indices differ
#'   by more than one are not treated as adjacent. Every scheduled prompt,
#'   answered or not, needs a row.
#' @return A list with \code{y} (matrix of states at \eqn{t}), \code{x}
#'   (states at \eqn{t - 1}), \code{id}, \code{time}, \code{row_t},
#'   \code{row_tm1} (row indices in \code{data}) for the complete pairs;
#'   \code{lag_obs}, a list describing every prompt whose predecessor is
#'   answered (its row \code{row_t}, its predecessor's row \code{row_tm1}, the
#'   predecessor's states \code{x}, its own response indicator \code{R},
#'   \code{id} and \code{time}), which response models on the observed history
#'   use; and \code{vars}, \code{n_prompts}, \code{response_rate} and
#'   \code{Rv} (the response indicator in the row order of \code{data}).
#' @examples
#' sim <- simulate_ema(N = 10, T = 12, motifs = "M1", seed = 1)
#' pr <- make_pairs(sim$data, sim$vars)
#' nrow(pr$y); pr$response_rate
#' # pairs are formed within days only when a day column is given
#' simd <- simulate_ema(N = 10, T = 12, days = 4, seed = 1)
#' nrow(make_pairs(simd$data, simd$vars)$y)
#' nrow(make_pairs(simd$data, simd$vars, day = "day")$y)
#' @export
make_pairs <- function(data, vars, id = "id", time = "time", day = NULL, R = "R") {
  check_columns(data, c(vars, id, time, day))
  check_time(data, time)
  Rv <- response_indicator(data, vars, R)
  ord <- order(data[[id]], data[[time]])
  data <- data[ord, , drop = FALSE]; Rv <- Rv[ord]
  rid <- data[[id]]; rt <- data[[time]]
  n <- nrow(data)
  same_person <- rid[-1] == rid[-n]
  consecutive <- rt[-1] == rt[-n] + 1
  same_day <- if (!is.null(day)) data[[day]][-1] == data[[day]][-n] else TRUE
  adj <- which(same_person & consecutive & same_day)     # pair (adj, adj + 1)
  lag_obs_rows <- adj[Rv[adj] == 1L]                       # predecessor answered
  both <- lag_obs_rows[Rv[lag_obs_rows + 1L] == 1L]
  Xm <- as.matrix(data[, vars, drop = FALSE])
  list(y = Xm[both + 1L, , drop = FALSE], x = Xm[both, , drop = FALSE],
       id = rid[both + 1L], time = rt[both + 1L], row_t = ord[both + 1L], row_tm1 = ord[both],
       lag_obs = list(row_t = ord[lag_obs_rows + 1L], row_tm1 = ord[lag_obs_rows],
                      x = Xm[lag_obs_rows, , drop = FALSE], R = Rv[lag_obs_rows + 1L],
                      id = rid[lag_obs_rows + 1L], time = rt[lag_obs_rows + 1L]),
       vars = vars, n_prompts = n, response_rate = mean(Rv), Rv = Rv[order(ord)])
}

# Internal: weighted within (person fixed effects) least squares of y on x.
within_fit <- function(y, x, g, w, levels_g) {
  p <- ncol(y); k <- ncol(x)
  gf <- factor(g, levels = levels_g); gi <- as.integer(gf)
  sw <- as.numeric(tapply(w, gf, sum)); sw[is.na(sw)] <- 0
  ybar <- matrix(apply(y, 2, function(v) tapply(w * v, gf, sum)), ncol = p) / sw
  xbar <- matrix(apply(x, 2, function(v) tapply(w * v, gf, sum)), ncol = k) / sw
  yd <- y - ybar[gi, , drop = FALSE]; xd <- x - xbar[gi, , drop = FALSE]
  Sxx <- crossprod(xd * w, xd); Sxy <- crossprod(xd * w, yd)
  Phi <- t(solve(Sxx, Sxy))
  res <- yd - xd %*% t(Phi)
  Sxx_inv <- solve(Sxx)
  G <- length(levels_g)
  scores <- matrix(0, G, p * k)
  for (j in seq_len(p)) {
    sj <- xd * (w * res[, j])
    scores[, (j - 1) * k + seq_len(k)] <- apply(sj, 2, function(v) tapply(v, gf, sum))
  }
  scores[is.na(scores)] <- 0
  infl <- scores %*% kronecker(diag(p), Sxx_inv)   # G x (p k) influence contributions to the stacked rows of the coefficient matrix
  list(Phi = Phi, ybar = ybar, xbar = xbar, res = res, infl = infl, sw = sw, n_i = as.numeric(table(gf)))
}

#' Within-person VAR(1) estimated from complete adjacent pairs
#'
#' Estimates the lagged-coefficient matrix \eqn{\Phi}, the innovation
#' covariance \eqn{\Psi}, person intercepts and person means from complete
#' adjacent prompt pairs, with person-specific intercepts (the within
#' estimator). This is the estimator that recovers the transition kernel under
#' the recoverable motifs (M0, M1, M3, M4 and their combinations, and M5 with
#' an observed context; see \code{\link{recoverability}}). Optional weights
#' turn it into the inverse-probability-weighted or tilted estimator. By
#' default the half-panel (split-panel) jackknife of Dhaene and Jochmans (2015)
#' removes the finite-\eqn{T} (Nickell) bias of the within estimator of
#' \eqn{\Phi}.
#'
#' @inheritParams make_pairs
#' @param weights Optional non-negative weights, one per complete pair (in the
#'   order returned by \code{\link{make_pairs}}).
#' @param pairs Optional precomputed result of \code{\link{make_pairs}}.
#' @param min_pairs Persons with fewer complete pairs are dropped from the
#'   between-person summaries (they still contribute to \eqn{\Phi}). A person
#'   can have at most \eqn{T - 1} pairs, so \code{min_pairs} must be smaller
#'   than the number of prompts per person.
#' @param bias_correct Logical; apply the half-panel jackknife to \eqn{\Phi},
#'   \eqn{\Psi}, the person intercepts and means, and the between-person
#'   covariance (all of which carry \eqn{O(1/T)} incidental-parameter bias).
#' @param covariates Optional names of observed context columns measured at
#'   prompt \eqn{t} that enter the within regression as covariates; the
#'   returned \eqn{\Phi} is then the context-conditional transition kernel
#'   (the recoverable estimand under M5 with an observed context), and
#'   \code{B_cov} holds the covariate coefficients. Person means are evaluated
#'   at each person's mean covariate level.
#' @return An object of class \code{"pairs_fit"}: a list with \code{Phi}
#'   (lagged coefficients, rows = outcome at \eqn{t}, columns = predictor at
#'   \eqn{t-1}), \code{Psi} (innovation covariance), \code{pcor}
#'   (contemporaneous partial correlations), \code{B_cov} (covariate
#'   coefficients, or \code{NULL}), \code{intercepts} (person intercepts
#'   \eqn{c_i}), \code{mu_dyn} (dynamics-recovered person means
#'   \eqn{(I - \Phi)^{-1} c_i}), \code{mu_obs} (observed person means),
#'   \code{mu} and \code{Sigma_mu} (person-weighted between-person law based
#'   on \code{mu_dyn}), \code{mu_obs_person} (person-weighted mean of the
#'   observed means), \code{mu_obs_pooled} (prompt-weighted mean of answered
#'   prompts), \code{Sigma_mu_obs}, the uncorrected (least-squares) versions
#'   \code{Phi_ls}, \code{Psi_ls}, \code{mu_ls}, \code{Sigma_mu_ls},
#'   \code{vcov_Phi} (cluster-robust covariance of the stacked rows of
#'   \code{Phi}), \code{n_pairs}, \code{n_persons} (persons with at least one
#'   complete pair), \code{n_pairs_person},
#'   \code{n_dropped}, \code{weights}, \code{pairs} (the
#'   \code{\link{make_pairs}} object), \code{vars}, \code{residuals} and
#'   \code{bias_correct}. Methods: \code{print}, \code{\link{summary.pairs_fit}},
#'   \code{coef}, \code{vcov}, \code{confint}, \code{nobs}.
#' @details Person intercepts are handled by weighted within-person demeaning
#'   on the sample of complete pairs (equivalent to person dummies), not by
#'   centering on the observed person means as in the two-step
#'   \code{mlVAR} approach; the two coincide for \eqn{\Phi} only in large
#'   samples. Persons contribute to \eqn{\Phi} and \eqn{\Psi} whatever their
#'   number of pairs; the between-person summaries use persons with at least
#'   \code{min_pairs} pairs, and the number dropped is reported in
#'   \code{n_dropped}.
#' @references Dhaene, G., & Jochmans, K. (2015). Split-panel jackknife
#'   estimation of fixed-effect models. \emph{The Review of Economic Studies,
#'   82}, 991-1030. \doi{10.1093/restud/rdv007}
#'
#'   Nickell, S. (1981). Biases in dynamic models with fixed effects.
#'   \emph{Econometrica, 49}, 1417-1426. \doi{10.2307/1911408}
#' @examples
#' sim <- simulate_ema(N = 40, T = 30, motifs = "M1", seed = 1)
#' f <- fit_pairs(sim$data, sim$vars)
#' f
#' summary(f)
#' coef(f)[1:4]; confint(f)[1:2, ]
#' # compare with the population values used by the simulator
#' round(f$Phi - default_params()$Phi, 2)
#' # covariate adjustment for an observed context (motif M5)
#' sim5 <- simulate_ema(N = 40, T = 30, motifs = "M5", seed = 2)
#' f5 <- fit_pairs(sim5$data, sim5$vars, covariates = "C")
#' f5$B_cov
#' @export
fit_pairs <- function(data, vars, id = "id", time = "time", day = NULL, R = "R",
                      weights = NULL, pairs = NULL, min_pairs = 5L, bias_correct = TRUE, covariates = NULL) {
  pr <- if (is.null(pairs)) make_pairs(data, vars, id, time, day, R) else pairs
  p <- length(vars)
  y <- pr$y; x <- pr$x; g <- pr$id
  q <- 0L
  if (!is.null(covariates)) {   # context-conditional kernel: covariates measured at prompt t enter the within regression
    Cm <- as.matrix(data[pr$row_t, covariates, drop = FALSE]); q <- ncol(Cm)
    x <- cbind(x, Cm)
  }
  n <- nrow(y)
  if (n < 2 * p + 2)
    stop(sprintf("Too few complete pairs: %d of %d prompts form a complete adjacent pair (at least %d are needed); check that '%s' is a consecutive integer prompt index within person%s",
                 n, pr$n_prompts, 2L * p + 2L, time, if (!is.null(day)) " and day" else ""))
  w <- if (is.null(weights)) rep(1, n) else as.numeric(weights)
  if (length(w) != n) stop("weights must have one entry per complete pair")
  if (any(!is.finite(w)) || any(w < 0)) stop("weights must be finite and non-negative")
  levels_g <- levels(factor(g)); G <- length(levels_g)
  # one fixed-effects fit and its derived quantities
  derive <- function(f, wsub, keep) {
    Psi <- crossprod(f$res * wsub, f$res) / max(sum(wsub) - G - p - q, 1)
    Psi <- (Psi + t(Psi)) / 2
    Phi_full <- f$Phi                      # p x (p + q): lagged states first, then covariates
    Phi <- Phi_full[, seq_len(p), drop = FALSE]
    B_cov <- if (q > 0) Phi_full[, p + seq_len(q), drop = FALSE] else NULL
    cI <- f$ybar - f$xbar %*% t(Phi_full)   # intercept net of lagged states and covariates
    if (q > 0) cI <- cI + f$xbar[, p + seq_len(q), drop = FALSE] %*% t(B_cov)   # evaluate at the person's mean covariate level
    mu_i <- cI %*% t(solve(diag(p) - Phi))
    mu <- colMeans(mu_i[keep, , drop = FALSE], na.rm = TRUE)
    Smu <- if (sum(keep) > 2) stats::cov(mu_i[keep, , drop = FALSE], use = "complete.obs") else matrix(NA_real_, p, p)
    infl <- f$infl
    if (q > 0) {   # keep the influence columns of the lagged-state coefficients only
      keep_cols <- as.vector(sapply(seq_len(p), function(k) (k - 1) * (p + q) + seq_len(p)))
      infl <- infl[, keep_cols, drop = FALSE]
    }
    list(Phi = Phi, Psi = Psi, cI = cI, mu_i = mu_i, mu = mu, Smu = Smu, infl = infl, B_cov = B_cov)
  }
  full <- within_fit(y, x, g, w, levels_g)
  n_i <- full$n_i
  keep <- n_i >= min_pairs
  if (!any(keep))
    warning("no person has at least min_pairs = ", min_pairs, " complete pairs (the maximum is ", max(n_i),
            "): the between-person summaries are NA; lower 'min_pairs'")
  est <- derive(full, w, keep)
  est_ls <- est
  if (bias_correct) {
    med <- stats::ave(pr$time, g, FUN = stats::median)
    A <- pr$time <= med; B <- !A
    okA <- sum(A) > 2 * p + 2; okB <- sum(B) > 2 * p + 2
    if (okA && okB) {
      fa <- within_fit(y[A, , drop = FALSE], x[A, , drop = FALSE], g[A], w[A], levels_g)
      fb <- within_fit(y[B, , drop = FALSE], x[B, , drop = FALSE], g[B], w[B], levels_g)
      ea <- derive(fa, w[A], keep); eb <- derive(fb, w[B], keep)
      jk <- function(nm) 2 * est[[nm]] - (ea[[nm]] + eb[[nm]]) / 2
      est <- list(Phi = jk("Phi"), Psi = jk("Psi"), cI = jk("cI"), mu_i = jk("mu_i"), mu = jk("mu"), Smu = jk("Smu"), infl = jk("infl"),
                  B_cov = if (q > 0) jk("B_cov") else NULL)
      # persons absent from one half keep their full-sample intercepts (no correction possible)
      bad <- !stats::complete.cases(est$cI)
      if (any(bad)) { est$cI[bad, ] <- est_ls$cI[bad, ]; est$mu_i[bad, ] <- est_ls$mu_i[bad, ] }
      est$infl[is.na(est$infl)] <- 0
      # keep Psi positive definite after the jackknife
      ev <- eigen(est$Psi, symmetric = TRUE)
      if (min(ev$values) < 1e-6) est$Psi <- ev$vectors %*% diag(pmax(ev$values, 1e-6), p) %*% t(ev$vectors)
    }
  }
  adj <- if (G > 1) G / (G - 1) else 1
  vcov_Phi <- adj * crossprod(est$infl)
  Phi <- est$Phi; Psi <- est$Psi; intercepts <- est$cI; mu_dyn <- est$mu_i
  Xall <- as.matrix(data[, vars, drop = FALSE])
  ok <- stats::complete.cases(Xall)
  idall <- factor(data[[id]][ok], levels = levels_g)
  mu_obs <- matrix(apply(Xall[ok, , drop = FALSE], 2, function(v) tapply(v, idall, mean)), ncol = p)
  mu_obs_pooled <- colMeans(Xall[ok, , drop = FALSE])
  pn <- as.vector(t(outer(vars, vars, function(a, b) paste0(a, "<-", b))))
  dimnames(vcov_Phi) <- list(pn, pn)
  dimnames(Phi) <- dimnames(Psi) <- list(vars, vars)
  rownames(intercepts) <- rownames(mu_dyn) <- rownames(mu_obs) <- levels_g
  colnames(intercepts) <- colnames(mu_dyn) <- colnames(mu_obs) <- vars
  structure(list(Phi = Phi, Phi_ls = est_ls$Phi, Psi = Psi, Psi_ls = est_ls$Psi, pcor = partial_cors(Psi), B_cov = est$B_cov,
                 intercepts = intercepts, mu_dyn = mu_dyn, mu_obs = mu_obs, mu = est$mu, Sigma_mu = est$Smu,
                 mu_ls = est_ls$mu, Sigma_mu_ls = est_ls$Smu,
                 mu_obs_pooled = mu_obs_pooled, mu_obs_person = colMeans(mu_obs[keep, , drop = FALSE], na.rm = TRUE),
                 Sigma_mu_obs = if (sum(keep) > 2) stats::cov(mu_obs[keep, , drop = FALSE], use = "complete.obs") else NULL,
                 vcov_Phi = vcov_Phi, n_pairs = n, n_persons = G, n_pairs_person = n_i, n_dropped = sum(!keep),
                 weights = w, pairs = pr, vars = vars, residuals = full$res, bias_correct = bias_correct),
            class = "pairs_fit")
}

#' @export
print.pairs_fit <- function(x, digits = 3, ...) {
  cat("Within-person VAR(1) from", x$n_pairs, "complete adjacent pairs,", x$n_persons, "persons",
      if (isTRUE(x$bias_correct)) "(half-panel jackknife)" else "", "\n")
  cat("Lagged coefficients Phi (rows = outcome at t, columns = predictor at t-1):\n")
  print(round(x$Phi, digits))
  cat("Contemporaneous partial correlations (innovations):\n")
  print(round(x$pcor, digits))
  cat("Between-person means (dynamics-recovered):", round(x$mu, digits), "\n")
  invisible(x)
}

#' Standard errors and confidence intervals for the lagged coefficients
#'
#' Cluster-robust (by person) standard errors with a \eqn{t(G - 1)} reference
#' distribution, \eqn{G} being the number of persons, as in
#' \code{\link{silence_test}}. With few persons (single-case or small
#' samples) the cluster-robust intervals are unreliable; a bootstrap over days
#' of the estimates (refitting \code{\link{fit_pairs}} on resampled day blocks)
#' is the alternative, and \code{\link{silence_test}} implements it for the
#' test.
#' @param object A \code{pairs_fit}.
#' @param level Confidence level.
#' @param ... Ignored.
#' @return \code{summary} returns a data frame with one row per lagged
#'   coefficient (\code{"Y<-X"} names the effect of \code{X} at \eqn{t-1} on
#'   \code{Y} at \eqn{t}) and the columns \code{coef}, \code{estimate},
#'   \code{se}, \code{lower}, \code{upper}, \code{z} (the ratio of the
#'   estimate to its standard error, referred to \eqn{t(G - 1)}) and \code{p}.
#' @references Cameron, A. C., & Miller, D. L. (2015). A practitioner's guide
#'   to cluster-robust inference. \emph{Journal of Human Resources, 50},
#'   317-372. \doi{10.3368/jhr.50.2.317}
#' @examples
#' sim <- simulate_ema(N = 30, T = 20, seed = 2)
#' f <- fit_pairs(sim$data, sim$vars)
#' summary(f)
#' coef(f)
#' confint(f, parm = "NegA<-NegA", level = 0.9)
#' nobs(f)
#' @export
summary.pairs_fit <- function(object, level = 0.95, ...) {
  est <- as.vector(t(object$Phi))
  se <- sqrt(diag(object$vcov_Phi))
  df <- max(object$n_persons - 1, 1)
  q <- stats::qt(1 - (1 - level) / 2, df)
  data.frame(coef = names(se), estimate = est, se = se, lower = est - q * se, upper = est + q * se,
             z = est / se, p = 2 * stats::pt(-abs(est / se), df), row.names = NULL)
}

#' @rdname summary.pairs_fit
#' @return \code{coef} returns the stacked rows of \eqn{\Phi} as a named
#'   vector, \code{vcov} their cluster-robust covariance matrix,
#'   \code{confint} a two-column matrix of limits and \code{nobs} the number of
#'   complete pairs.
#' @export
coef.pairs_fit <- function(object, ...) stats::setNames(as.vector(t(object$Phi)), rownames(object$vcov_Phi))

#' @rdname summary.pairs_fit
#' @export
vcov.pairs_fit <- function(object, ...) object$vcov_Phi

#' @rdname summary.pairs_fit
#' @param parm Coefficient names or indices (default: all).
#' @export
confint.pairs_fit <- function(object, parm, level = 0.95, ...) {
  s <- summary.pairs_fit(object, level = level)
  m <- as.matrix(s[, c("lower", "upper")]); rownames(m) <- s$coef
  colnames(m) <- paste0(format(100 * c((1 - level) / 2, 1 - (1 - level) / 2), trim = TRUE), " %")
  if (!missing(parm)) m <- m[parm, , drop = FALSE]
  m
}

#' @rdname summary.pairs_fit
#' @export
nobs.pairs_fit <- function(object, ...) object$n_pairs
