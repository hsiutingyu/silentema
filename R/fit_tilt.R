#' Tilted (self-censoring-adjusted) within-person VAR(1) at a fixed sensitivity value
#'
#' Sensitivity analysis for motif M2 (self-censoring). The selection model
#' \eqn{P(R_t = 1 \mid X_t = x) = \Phi_N(\alpha + \delta' x)} is held at a fixed
#' sensitivity vector \eqn{\delta} that is never estimated from the fit. The
#' response intercept is calibrated so that the model-implied response rate
#' given the observed history matches the observed one, and complete pairs
#' are reweighted by \eqn{w_t = P(R_t = 1 \mid x_{t-1}) / P(R_t = 1 \mid x_t)}
#' (inverse-probability weighting with a stabilized numerator). Weighted
#' within-person regression and intercept calibration are iterated to a
#' fixed point; the population parameters are a fixed point when
#' \eqn{\delta} is the true selection vector.
#'
#' @inheritParams fit_pairs
#' @param delta Numeric vector (length = number of variables) of
#'   sensitivity coefficients in probit units per unit of the state; a scalar
#'   is applied to the first variable. Negative values encode "high states
#'   are skipped".
#' @param propensity \code{"common"} (one intercept) or \code{"person"}
#'   (one intercept per person; use when a person-propensity motif M4 is
#'   declared together with M2).
#' @param min_pairs Passed to \code{\link{fit_pairs}}.
#' @param probe Optional name of a randomized-probe indicator column; probe
#'   prompts are answered by design, receive weight one and are excluded from
#'   the intercept calibration.
#' @param max_iter,tol Fixed-point iteration control.
#' @param start Optional \code{pairs_fit} (typically the fit at a neighboring
#'   sensitivity value) used as the starting point of the iteration.
#' @param damping Step size of the damped fixed-point update (1 = undamped).
#' @return An object of class \code{c("tilt_fit", "pairs_fit")}: a
#'   \code{\link{fit_pairs}} object with additional elements \code{delta},
#'   \code{alpha} (calibrated intercept(s)), \code{iterations},
#'   \code{converged}, \code{weights}, \code{max_weight}, and \code{ess}, the
#'   effective number of complete pairs \eqn{(\sum w)^2 / \sum w^2} (Kish).
#'   The weights are not truncated: profiles in which \code{ess} falls far
#'   below the number of pairs, or in which a few pairs carry very large
#'   weights, rest on few observations and should be read with caution. The
#'   cluster-robust covariance treats the weights as fixed (the uncertainty of
#'   the calibrated intercept and of the parameters entering the weights is
#'   not propagated).
#' @seealso \code{\link{tilt_profile}} fits a grid of sensitivity values;
#'   \code{\link{calibrate_delta}} chooses one from a design feature.
#' @examples
#' sim <- simulate_ema(N = 40, T = 30, motifs = "M2", delta = -1, seed = 1)
#' f0 <- fit_pairs(sim$data, sim$vars)      # treats skips as ignorable
#' f <- fit_tilt(sim$data, sim$vars, delta = -1)
#' f
#' f$ess / f$n_pairs                          # effective share of pairs
#' c(untilted = f0$Phi[1, 1], tilted = f$Phi[1, 1], truth = default_params()$Phi[1, 1])
#' # person-specific response intercepts (motif M4 declared with M2)
#' fp <- fit_tilt(sim$data, sim$vars, delta = -1, propensity = "person")
#' summary(fp$alpha)
#' @export
fit_tilt <- function(data, vars, delta, id = "id", time = "time", day = NULL, R = "R",
                     propensity = c("common", "person"), pairs = NULL, probe = NULL, max_iter = 200L, tol = 1e-5,
                     start = NULL, damping = 0.5, min_pairs = 5L) {
  propensity <- match.arg(propensity)
  p <- length(vars)
  if (!is.numeric(delta) || !(length(delta) %in% c(1L, p)))
    stop("'delta' must be a single number (applied to the first variable) or a vector with one entry per variable (", p, ")")
  if (length(delta) == 1L) { d <- rep(0, p); d[1] <- delta; delta <- d }
  if (!is.null(probe) && !(probe %in% names(data))) stop("probe column '", probe, "' not found in 'data'")
  pr <- if (is.null(pairs)) make_pairs(data, vars, id, time, day, R) else pairs
  if (all(delta == 0) || is.null(start)) fit0 <- fit_pairs(data, vars, id, time, day, R, pairs = pr, min_pairs = min_pairs)
  fit <- if (is.null(start)) fit0 else start
  if (all(delta == 0)) { fit0$delta <- delta; fit0$alpha <- NA; fit0$iterations <- 0L; fit0$converged <- TRUE
    fit0$ess <- fit0$n_pairs; fit0$max_weight <- 1; fit0$propensity <- propensity
    class(fit0) <- c("tilt_fit", "pairs_fit"); return(fit0) }
  lo <- pr$lag_obs                       # prompts with answered predecessor (response known)
  Zlo <- if (!is.null(probe) && probe %in% names(data)) data[[probe]][lo$row_t] == 1L else rep(FALSE, length(lo$R))
  Zpair <- if (!is.null(probe) && probe %in% names(data)) data[[probe]][pr$row_t] == 1L else rep(FALSE, nrow(pr$y))
  if (any(Zlo)) {   # drop probe prompts from the calibration set
    lo <- list(row_t = lo$row_t[!Zlo], row_tm1 = lo$row_tm1[!Zlo], x = lo$x[!Zlo, , drop = FALSE],
               R = lo$R[!Zlo], id = lo$id[!Zlo], time = lo$time[!Zlo])
  }
  glev <- rownames(fit$intercepts)
  gi_lo <- match(as.character(lo$id), glev)
  if (anyNA(gi_lo)) {   # persons without any complete pair have no intercept: drop from calibration
    keep <- !is.na(gi_lo)
    lo <- list(row_t = lo$row_t[keep], row_tm1 = lo$row_tm1[keep], x = lo$x[keep, , drop = FALSE],
               R = lo$R[keep], id = lo$id[keep], time = lo$time[keep])
    gi_lo <- gi_lo[keep]
  }
  gi_pair <- match(as.character(pr$id), glev)
  conv <- FALSE
  for (it in seq_len(max_iter)) {
    Phi <- fit$Phi; Psi <- fit$Psi; cI <- fit$intercepts
    kappa <- sqrt(1 + as.numeric(t(delta) %*% Psi %*% delta))
    # model-implied conditional means for prompts with answered predecessor
    m_lo <- cI[gi_lo, , drop = FALSE] + lo$x %*% t(Phi)
    dm_lo <- as.vector(m_lo %*% delta)
    # calibrate alpha: mean(R) = mean(pnorm((alpha + delta'm)/kappa))
    if (propensity == "common") {
      f <- function(a) mean(stats::pnorm((a + dm_lo) / kappa)) - mean(lo$R)
      alpha <- stats::uniroot(f, c(-30, 30), tol = 1e-8)$root
      alpha_vec_lo <- rep(alpha, length(dm_lo))
    } else {
      # person-specific intercepts by vectorized Newton steps on the calibration equations
      G <- length(glev)
      n_g <- tabulate(gi_lo, G); r_g <- as.numeric(rowsum(lo$R, gi_lo, reorder = TRUE)); r_g <- r_g / pmax(n_g, 1)
      r_g <- pmin(pmax(r_g, 0.5 / pmax(n_g, 1)), 1 - 0.5 / pmax(n_g, 1))      # keep the equation solvable
      alpha <- if (exists("alpha_prev", inherits = FALSE) && length(alpha_prev) == G) alpha_prev else
        rep(kappa * stats::qnorm(mean(lo$R)) - mean(dm_lo), G)
      for (nw in 1:25) {
        zz <- (alpha[gi_lo] + dm_lo) / kappa
        f <- as.numeric(rowsum(stats::pnorm(zz), gi_lo, reorder = TRUE)) / pmax(n_g, 1) - r_g
        fp <- as.numeric(rowsum(stats::dnorm(zz), gi_lo, reorder = TRUE)) / pmax(n_g, 1) / kappa
        step <- ifelse(fp > 1e-10, f / fp, 0)
        step <- pmax(pmin(step, 1), -1)
        alpha <- alpha - step
        if (max(abs(step[n_g > 0])) < 1e-8) break
      }
      alpha[n_g == 0] <- NA_real_
      alpha_prev <- alpha
      alpha_vec_lo <- alpha[gi_lo]
    }
    # weights on complete pairs
    m_pair <- cI[gi_pair, , drop = FALSE] + pr$x %*% t(Phi)
    a_pair <- if (propensity == "common") rep(alpha, nrow(m_pair)) else alpha[gi_pair]
    num <- stats::pnorm((a_pair + as.vector(m_pair %*% delta)) / kappa)
    den <- stats::pnorm(a_pair + as.vector(pr$y %*% delta))
    w <- num / pmax(den, 1e-12)
    w[Zpair] <- 1
    new <- fit_pairs(data, vars, id, time, day, R, weights = w, pairs = pr, min_pairs = min_pairs)
    change <- max(abs(new$Phi - fit$Phi), abs(new$Psi - fit$Psi), abs(new$intercepts - fit$intercepts), na.rm = TRUE)
    last_undamped <- new
    if (change < tol) { conv <- TRUE; fit <- new; break }
    if (damping < 1 && it > 1) {   # damped update of the quantities that feed the weights
      new$Phi <- damping * new$Phi + (1 - damping) * fit$Phi
      new$Psi <- damping * new$Psi + (1 - damping) * fit$Psi
      new$intercepts <- damping * new$intercepts + (1 - damping) * fit$intercepts
      new$mu_dyn <- new$intercepts %*% t(solve(diag(p) - new$Phi))
      keep <- new$n_pairs_person >= min_pairs
      new$mu <- colMeans(new$mu_dyn[keep, , drop = FALSE], na.rm = TRUE)
      new$pcor <- partial_cors(new$Psi)
    }
    fit <- new
  }
  if (!conv) fit <- last_undamped
  # every returned quantity (Phi, Psi, intercepts, means, vcov) comes from one undamped weighted fit,
  # whose weights were computed from the converged (or, without convergence, the last) iterate
  fit$delta <- delta; fit$alpha <- alpha; fit$iterations <- it; fit$converged <- conv
  fit$ess <- sum(fit$weights)^2 / sum(fit$weights^2)
  fit$max_weight <- max(fit$weights)
  fit$propensity <- propensity
  class(fit) <- c("tilt_fit", "pairs_fit")
  fit
}

#' Sensitivity profile over a grid of self-censoring values
#'
#' Fits \code{\link{fit_tilt}} over a grid of sensitivity values and collects
#' the lagged coefficients, the contemporaneous partial correlations, the
#' between-person means and (optionally) cluster-robust confidence limits.
#'
#' @inheritParams fit_tilt
#' @param delta_grid Numeric vector of sensitivity values for the
#'   self-censoring variable (the first of \code{vars} unless \code{which}
#'   says otherwise).
#' @param which Index (in \code{vars}) of the variable that drives
#'   self-censoring.
#' @param level Confidence level for the limits.
#' @param probe Optional probe column name passed to \code{\link{fit_tilt}}.
#' @return An object of class \code{"tilt_profile"}: a list with
#'   \code{delta_grid}, \code{which}, \code{vars}, a data frame \code{Phi}
#'   (one row per grid value and lagged coefficient, with \code{coef},
#'   \code{estimate}, \code{se}, \code{lower}, \code{upper}, \code{z},
#'   \code{p} and \code{delta}), a data frame \code{pcor} (contemporaneous
#'   partial correlations by grid value and edge), a matrix \code{mu}
#'   (between-person means by grid value), \code{fits} (the list of
#'   \code{\link{fit_tilt}} objects), \code{level}, \code{propensity},
#'   \code{converged} and \code{ess} (one entry per grid value) and
#'   \code{n_pairs}. Methods: \code{print} and
#'   \code{\link[=plot.tilt_profile]{plot}}.
#' @seealso \code{\link{break_even}}, \code{\link{calibrate_delta}},
#'   \code{\link{plot.tilt_profile}}
#' @examples
#' sim <- simulate_ema(N = 40, T = 30, motifs = "M2", delta = -1, seed = 1)
#' prof <- tilt_profile(sim$data, sim$vars, delta_grid = c(-1.5, -1, -0.5, 0, 0.5))
#' prof
#' head(prof$Phi)
#' prof$mu
#' break_even(prof, delta_max = 1)
#' plot(prof, plausible = c(-1.5, 0))
#' @export
tilt_profile <- function(data, vars, delta_grid = seq(-2, 2, by = 0.5), which = 1L,
                         id = "id", time = "time", day = NULL, R = "R",
                         propensity = c("common", "person"), level = 0.95, pairs = NULL, probe = NULL, min_pairs = 5L) {
  propensity <- match.arg(propensity)
  p <- length(vars)
  if (!is.numeric(which) || length(which) != 1L || which < 1 || which > p) stop("'which' must index one of 'vars'")
  pr <- if (is.null(pairs)) make_pairs(data, vars, id, time, day, R) else pairs
  fits <- vector("list", length(delta_grid))
  rows <- list(); prow <- list(); mus <- matrix(NA_real_, length(delta_grid), p, dimnames = list(NULL, vars))
  # warm starts: fit from the grid value closest to zero outward on each side
  ord <- order(abs(delta_grid), delta_grid)
  last_neg <- NULL; last_pos <- NULL
  for (k in ord) {
    d <- rep(0, p); d[which] <- delta_grid[k]
    st <- if (delta_grid[k] < 0) last_neg else if (delta_grid[k] > 0) last_pos else NULL
    f <- fit_tilt(data, vars, d, id, time, day, R, propensity = propensity, pairs = pr, probe = probe, start = st, min_pairs = min_pairs)
    if (delta_grid[k] < 0) last_neg <- f else if (delta_grid[k] > 0) last_pos <- f else { last_neg <- f; last_pos <- f }
    fits[[k]] <- f
    s <- summary.pairs_fit(f, level = level)
    s$delta <- delta_grid[k]
    rows[[k]] <- s
    pc <- f$pcor; ut <- upper.tri(pc)
    prow[[k]] <- if (p > 1) data.frame(delta = delta_grid[k],
                            edge = paste(vars[row(pc)[ut]], vars[col(pc)[ut]], sep = "--"),
                            estimate = pc[ut]) else data.frame(delta = numeric(0), edge = character(0), estimate = numeric(0))
    mus[k, ] <- f$mu
  }
  Phi_df <- do.call(rbind, rows); Phi_df <- Phi_df[order(Phi_df$delta, Phi_df$coef), ]
  pc_df <- do.call(rbind, prow); if (nrow(pc_df)) pc_df <- pc_df[order(pc_df$delta, pc_df$edge), ]
  structure(list(delta_grid = delta_grid, which = which, vars = vars, Phi = Phi_df,
                 pcor = pc_df, mu = mus, fits = fits, level = level,
                 propensity = propensity, converged = vapply(fits, function(f) isTRUE(f$converged), TRUE),
                 ess = vapply(fits, function(f) if (is.null(f$ess)) NA_real_ else f$ess, 0), n_pairs = nrow(pr$y)),
            class = "tilt_profile")
}

#' @export
print.tilt_profile <- function(x, digits = 3, ...) {
  cat("Tilt profile over delta in {", paste(x$delta_grid, collapse = ", "), "} for", x$vars[x$which], "\n")
  ar <- x$Phi[x$Phi$coef == paste0(x$vars[x$which], "<-", x$vars[x$which]), c("delta", "estimate", "lower", "upper")]
  ar[] <- lapply(ar, round, digits)
  ar$ess <- round(x$ess[match(ar$delta, x$delta_grid)]); ar$converged <- x$converged[match(ar$delta, x$delta_grid)]
  cat("Autoregressive coefficient of the self-censoring variable along the grid (", x$n_pairs, "complete pairs):\n")
  print(ar, row.names = FALSE)
  if (!all(x$converged)) cat("Warning: fixed-point iteration did not converge for some grid values\n")
  invisible(x)
}

#' Break-even sensitivity values and identified sets from a tilt profile
#'
#' For each lagged coefficient reports the smallest \eqn{|\delta|} on the grid
#' at which the sign of the estimate changes or the confidence interval starts
#' (or stops) covering zero, and the identified set (range of the estimate)
#' over \eqn{|\delta| \le \delta_{max}}.
#'
#' @param profile A \code{tilt_profile}.
#' @param delta_max Half-width of the plausible sensitivity interval.
#' @return A data frame with one row per lagged coefficient and the columns
#'   \code{coef}, \code{estimate_0} (the estimate at \eqn{\delta = 0}),
#'   \code{significant_0} (whether its confidence interval excludes zero
#'   there), \code{sign_flip_delta} and \code{significance_flip_delta} (the
#'   grid values nearest zero at which the sign and the significance status
#'   change; \code{NA} = never on the grid), \code{set_lower} and
#'   \code{set_upper} (the identified set: the range of the estimate over
#'   \eqn{|\delta| \le \delta_{max}}) and \code{band_lower} and
#'   \code{band_upper} (the band: the range of the confidence limits over the
#'   same interval). A significance change is a dichotomous criterion; the sign
#'   change and the band are the primary quantities.
#' @examples
#' sim <- simulate_ema(N = 40, T = 30, motifs = "M2", delta = -1, seed = 1)
#' prof <- tilt_profile(sim$data, sim$vars, delta_grid = c(-1, -0.5, 0, 0.5))
#' be <- break_even(prof, delta_max = 1)
#' be[be$coef == "NegA<-NegA", ]
#' # a narrower plausible interval gives a narrower identified set
#' break_even(prof, delta_max = 0.5)[1, c("set_lower", "set_upper")]
#' @export
break_even <- function(profile, delta_max = max(abs(profile$delta_grid))) {
  d <- profile$Phi
  out <- lapply(split(d, d$coef), function(s) {
    s <- s[order(s$delta), ]
    at0 <- s[which.min(abs(s$delta)), ]
    sign0 <- sign(at0$estimate); sig0 <- !(at0$lower <= 0 & at0$upper >= 0)
    sign_flip <- s$delta[sign(s$estimate) != sign0 & s$estimate != 0]
    sig_flip <- s$delta[(!(s$lower <= 0 & s$upper >= 0)) != sig0]
    inset <- abs(s$delta) <= delta_max
    data.frame(coef = s$coef[1], estimate_0 = at0$estimate, significant_0 = sig0,
               sign_flip_delta = if (length(sign_flip)) sign_flip[which.min(abs(sign_flip))] else NA_real_,
               significance_flip_delta = if (length(sig_flip)) sig_flip[which.min(abs(sig_flip))] else NA_real_,
               set_lower = min(s$estimate[inset]), set_upper = max(s$estimate[inset]),
               band_lower = min(s$lower[inset]), band_upper = max(s$upper[inset]))
  })
  out <- do.call(rbind, out); rownames(out) <- NULL
  out
}

#' Plot a tilt profile
#'
#' Draws the estimates of selected lagged coefficients (with confidence
#' bands) against the sensitivity value, marks the plausible interval when
#' supplied, and marks a calibrated value when supplied.
#'
#' @param x A \code{tilt_profile}.
#' @param coefs Character vector of coefficient names (\code{"Y<-X"} form);
#'   default: the row of the self-censoring variable.
#' @param plausible Optional numeric vector of length 2 (plausible interval).
#' @param calibrated Optional numeric value(s) to mark.
#' @param ... Passed to \code{plot}.
#' @return Invisibly, the data frame of plotted estimates and limits (the
#'   rows of \code{x$Phi} for the selected coefficients). Called for its side
#'   effect, the plot.
#' @examples
#' sim <- simulate_ema(N = 40, T = 30, motifs = "M2", delta = -1, seed = 1)
#' prof <- tilt_profile(sim$data, sim$vars, delta_grid = c(-1, -0.5, 0, 0.5))
#' plot(prof, plausible = c(-1, 0), calibrated = -0.6)
#' # cross-lagged effects on stress instead of the row of the self-censoring variable
#' plot(prof, coefs = c("Stress<-NegA", "Stress<-Fatigue"), main = "Effects on Stress")
#' @export
plot.tilt_profile <- function(x, coefs = NULL, plausible = NULL, calibrated = NULL, ...) {
  d <- x$Phi
  if (is.null(coefs)) coefs <- paste0(x$vars, "<-", x$vars[x$which])
  d <- d[d$coef %in% coefs, ]
  cols <- grDevices::hcl.colors(length(coefs), "Dark 3")
  graphics::plot(range(d$delta), range(c(d$lower, d$upper)), type = "n",
                 xlab = expression(paste("sensitivity value ", delta, " (probit units per unit state)")),
                 ylab = "lagged coefficient", ...)
  if (!is.null(plausible)) graphics::rect(plausible[1], graphics::par("usr")[3], plausible[2], graphics::par("usr")[4], col = "grey92", border = NA)
  graphics::abline(h = 0, lty = 3); graphics::abline(v = 0, lty = 3)
  for (k in seq_along(coefs)) {
    s <- d[d$coef == coefs[k], ]; s <- s[order(s$delta), ]
    graphics::polygon(c(s$delta, rev(s$delta)), c(s$lower, rev(s$upper)), col = grDevices::adjustcolor(cols[k], 0.15), border = NA)
    graphics::lines(s$delta, s$estimate, col = cols[k], lwd = 2)
  }
  if (!is.null(calibrated)) graphics::abline(v = calibrated, lty = 2, lwd = 1.5)
  graphics::legend("topleft", legend = coefs, col = cols, lwd = 2, bty = "n", cex = 0.8, bg = NA)
  invisible(d)
}
