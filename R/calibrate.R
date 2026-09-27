# Internal: run an expression with a temporary seed, restoring the caller's RNG state afterwards
with_seed <- function(seed, expr) {
  if (is.null(seed)) return(expr)
  had <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  old <- if (had) get(".Random.seed", envir = globalenv(), inherits = FALSE) else NULL
  on.exit({
    if (had) assign(".Random.seed", old, envir = globalenv()) else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) rm(".Random.seed", envir = globalenv())
  }, add = TRUE)
  set.seed(seed)
  expr
}

# Internal: model-implied conditional means of X_t given x_{t-1} and R_t under a tilt fit
implied_means <- function(fit, x_lag, gi) {
  Phi <- fit$Phi; Psi <- fit$Psi; delta <- fit$delta; cI <- fit$intercepts
  kappa <- sqrt(1 + as.numeric(t(delta) %*% Psi %*% delta))
  m <- cI[gi, , drop = FALSE] + x_lag %*% t(Phi)
  a <- if (length(fit$alpha) == 1L) rep(fit$alpha, nrow(m)) else fit$alpha[gi]
  z <- (a + as.vector(m %*% delta)) / kappa
  lam <- function(v) exp(stats::dnorm(v, log = TRUE) - stats::pnorm(v, log.p = TRUE))
  shift <- as.vector(Psi %*% delta) / kappa
  list(m = m, m_R1 = m + outer(lam(z), shift), m_R0 = m - outer(lam(-z), shift), z = z, pi = stats::pnorm(z))
}

# Internal: pooled response persistence, P(R_t = 1 | R_{t-1} = 1) - P(R_t = 1 | R_{t-1} = 0),
# over consecutive prompts (within day when a day column is given)
persistence_contrast <- function(data, id = "id", time = "time", day = NULL, R = "R") {
  ord <- order(data[[id]], data[[time]]); data <- data[ord, , drop = FALSE]
  rid <- data[[id]]; rt <- data[[time]]; Rv <- data[[R]]
  n <- nrow(data); i <- seq_len(n - 1L)
  ok <- rid[i] == rid[i + 1L] & rt[i + 1L] == rt[i] + 1L
  if (!is.null(day)) ok <- ok & data[[day]][i] == data[[day]][i + 1L]
  i <- i[ok]
  prev <- Rv[i]; cur <- Rv[i + 1L]
  if (!any(prev == 0L) || !any(prev == 1L)) return(NA_real_)
  mean(cur[prev == 1L]) - mean(cur[prev == 0L])
}

#' Calibrate the self-censoring sensitivity parameter from a design feature
#'
#' Maps an observed statistic that is sensitive to self-censoring onto the
#' tilt profile and returns the value of the sensitivity parameter at which
#' the fitted model reproduces it. Three statistics are supported: the sensor gap (the
#' coefficient of \eqn{R_t} in \code{\link{sensor_gap_test}}; needs an
#' always-observed sensor), the probe contrast (the within-person difference
#' between the states reported at randomized forced-response probes and at
#' ordinary answered prompts; needs probes), and the post-skip contrast of the
#' \code{\link{silence_test}} (needs no extra channel but relies entirely on
#' the functional form). All three are parametric calibrations, not
#' identification results (Yu, 2026).
#'
#' The post-skip calibration compares the observed post-skip contrast with the
#' contrast implied by data simulated from the fitted tilt model at every grid
#' value. When burden (motif M3) is declared together with self-censoring, the
#' post-skip contrast also carries the collider contribution of
#' \eqn{R_t \to R_{t+1} \leftarrow X_{t+1}}{R_t -> R_{t+1} <- X_{t+1}}, which works against the
#' self-censoring signal; with \code{burden = "fit"} the simulated model
#' includes a burden term whose size is calibrated, at every grid value, so
#' that the simulated response persistence (the response rate after an
#' answered minus after a skipped prompt) matches the observed persistence.
#' Without a burden term the model-implied curve is the M2-only curve, which
#' understates \eqn{|\delta|} when burden is present.
#'
#' @param profile A \code{\link{tilt_profile}}.
#' @param data The data used for the profile.
#' @param method One of \code{"sensor"}, \code{"probe"}, \code{"postskip"}.
#' @param sensor,probe Names of the sensor and probe columns.
#' @param id,time,day,R Column names as in \code{\link{make_pairs}}.
#' @param n_sim Number of simulated data sets per grid value (and per burden
#'   value) for \code{method = "postskip"}; the Monte Carlo error of the
#'   model-implied curve is roughly the standard error of the post-skip
#'   coefficient divided by \code{sqrt(n_sim)}.
#' @param burden For \code{method = "postskip"}: \code{"none"} simulates from
#'   the fitted self-censoring model alone; \code{"fit"} adds a burden term
#'   calibrated to the observed response persistence (use when M3 is declared).
#' @param kappa_grid Grid of burden values (probit drop in the response
#'   propensity after a skipped prompt) searched by \code{burden = "fit"}.
#' @param seed Optional seed for the post-skip simulation and the bootstrap;
#'   the caller's random-number state is restored afterwards. \code{NULL}
#'   uses (and advances) the current stream.
#' @param level Confidence level for the interval obtained by inverting the
#'   observed statistic's confidence limits through the predicted curve (the
#'   uncertainty of the curve itself is ignored; Yu, 2026, reports the
#'   empirical coverage of this interval).
#' @param boot Number of cluster-bootstrap resamples for a percentile
#'   interval that refits the profile in every resample (0 = none).
#' @param boot_unit Resampling unit: persons, or days (for single-person data).
#' @param boot_grid Optional grid for the bootstrap profiles; by default seven
#'   points around the point estimate.
#' @return An object of class \code{"delta_calibration"}: a list with
#'   \code{method}, \code{delta_hat} (the calibrated value; \code{NA} when
#'   the model-implied curve does not cross the observed value on the grid),
#'   \code{interval} (the \code{level} interval obtained by inverting the
#'   confidence limits of the observed statistic; an endpoint outside the range
#'   of the curve leaves that side open, \code{-Inf} or \code{Inf}),
#'   \code{n_crossings} (when the curve crosses more than once the crossing
#'   nearest zero is used), \code{observed} and \code{se} (the observed
#'   statistic and its standard error), \code{curve} (a data frame of the
#'   grid values and the model-implied statistics, with the calibrated burden
#'   value \code{kappa_hat} per grid point when \code{burden = "fit"} and a
#'   \code{note} column when a grid point could not be evaluated),
#'   \code{variable}, \code{level}, \code{burden},
#'   \code{persistence_observed}, and, when \code{boot > 0}, \code{boot}
#'   (the bootstrap replicates of \code{delta_hat}), \code{boot_se} and
#'   \code{boot_interval} (percentile interval). Methods: \code{print} and
#'   \code{\link[=plot.delta_calibration]{plot}}.
#' @seealso \code{\link{tilt_profile}}, \code{\link{sensor_gap_test}},
#'   \code{\link{silence_test}}, \code{\link{simulate_from_fit}}
#' @examples
#' sim <- simulate_ema(N = 40, T = 30, motifs = "M2", delta = -1, sensor_cor = 0.6,
#'                     p_probe = 0.05, seed = 3)
#' vars <- sim$vars
#' prof <- tilt_profile(sim$data, vars, delta_grid = c(-1.5, -1, -0.5, 0), probe = "Z")
#' cs <- calibrate_delta(prof, sim$data, method = "sensor")
#' cs
#' cs$curve
#' calibrate_delta(prof, sim$data, method = "probe")
#' # the post-skip calibration simulates from every fitted model
#' calibrate_delta(prof, sim$data, method = "postskip", n_sim = 5, seed = 1)
#' @export
calibrate_delta <- function(profile, data, method = c("sensor", "probe", "postskip"),
                            sensor = "S", probe = "Z", id = "id", time = "time", day = NULL, R = "R",
                            n_sim = 20L, burden = c("none", "fit"), kappa_grid = c(0, 0.5, 1, 1.5, 2.5),
                            seed = NULL, level = 0.95, boot = 0L, boot_unit = c("person", "day"),
                            boot_grid = NULL) {
  method <- match.arg(method); boot_unit <- match.arg(boot_unit); burden <- match.arg(burden)
  if (!inherits(profile, "tilt_profile")) stop("'profile' must be a tilt_profile object")
  vars <- profile$vars; wh <- profile$which
  grid <- profile$delta_grid
  pr <- make_pairs(data, vars, id, time, day, R)
  glev <- rownames(profile$fits[[1]]$intercepts)
  pred <- rep(NA_real_, length(grid)); kappa_hat <- rep(NA_real_, length(grid)); note <- rep("", length(grid))
  if (method == "sensor") {
    if (!(sensor %in% names(data))) stop("sensor column '", sensor, "' not found in 'data'")
    obs <- sensor_gap_test(data, vars, sensor, id, time, day, R)
    lo <- pr$lag_obs; gi <- match(as.character(lo$id), glev)
    keep <- !is.na(gi)
    lo <- lapply(lo, function(v) if (is.matrix(v)) v[keep, , drop = FALSE] else v[keep]); gi <- gi[keep]
    X <- cbind(lo$x, R_t = lo$R)
    for (k in seq_along(grid)) {
      f <- profile$fits[[k]]
      if (all(f$delta == 0)) { pred[k] <- 0; next }      # no tilt: the model-implied gap is exactly zero
      b <- attr(obs, "loading")
      im <- implied_means(f, lo$x, gi)
      Stil <- ifelse(lo$R == 1L, as.vector(im$m_R1 %*% b), as.vector(im$m_R0 %*% b))
      fr <- within_reg(matrix(Stil, ncol = 1), X, lo$id)
      pred[k] <- fr$coef[ncol(X), 1]
    }
    observed <- obs$coef_R; se <- obs$se
  } else if (method == "probe") {
    if (!(probe %in% names(data))) stop("probe column '", probe, "' not found in 'data'")
    Z <- data[[probe]][pr$row_t]
    gi <- match(as.character(pr$id), glev)
    X <- cbind(pr$x, Z_t = Z)
    fr_obs <- within_reg(pr$y[, wh, drop = FALSE], X, pr$id)
    observed <- fr_obs$coef[ncol(X), 1]; se <- sqrt(fr_obs$vcov[ncol(X), ncol(X)])
    for (k in seq_along(grid)) {
      f <- profile$fits[[k]]
      if (all(f$delta == 0)) { pred[k] <- 0; next }      # no tilt: the model-implied contrast is exactly zero
      im <- implied_means(f, pr$x, gi)
      Ytil <- ifelse(Z == 1L, im$m[, wh], im$m_R1[, wh])
      fr <- within_reg(matrix(Ytil, ncol = 1), X, pr$id)
      pred[k] <- fr$coef[ncol(X), 1]
    }
  } else {
    st <- silence_test(data, vars, id, time, day, R)
    observed <- st$coef_R[wh]; se <- st$se[wh]
    N <- length(unique(data[[id]])); Tn <- max(table(data[[id]]))
    days <- if (!is.null(day)) as.integer(round(stats::median(table(data[[id]], data[[day]])[table(data[[id]], data[[day]]) > 0]))) else NULL
    persist_obs <- if (burden == "fit") persistence_contrast(data, id, time, day, R) else NA_real_
    if (burden == "fit" && is.na(persist_obs)) stop("burden = 'fit' needs both answered and skipped predecessors in the data")
    kg <- if (burden == "fit") sort(unique(c(0, kappa_grid))) else 0
    sim_stats <- function(f, kappa, k) {   # mean post-skip coefficient and persistence over n_sim simulated data sets
      vals <- rep(NA_real_, n_sim); pers <- rep(NA_real_, n_sim)
      for (b in seq_len(n_sim)) {
        sim <- tryCatch(simulate_from_fit(f, N = N, T = Tn, days = days, kappa_R = kappa,
                                          seed = if (is.null(seed)) NULL else seed + 100000L * k + 1000L * match(kappa, kg) + b),
                        error = function(e) { note[k] <<- conditionMessage(e); NULL })
        if (is.null(sim)) next
        stb <- tryCatch(silence_test(sim, vars, "id", "time", if (is.null(days)) NULL else "day", "R", se = "cluster"),
                        error = function(e) NULL)
        if (!is.null(stb)) vals[b] <- stb$coef_R[wh]
        if (burden == "fit") pers[b] <- persistence_contrast(sim, "id", "time", if (is.null(days)) NULL else "day", "R")
      }
      c(post = mean(vals, na.rm = TRUE), persist = mean(pers, na.rm = TRUE))
    }
    for (k in seq_along(grid)) {
      f <- profile$fits[[k]]
      if (burden == "none") {
        pred[k] <- sim_stats(f, 0, k)["post"]
      } else {
        S <- t(sapply(kg, function(kap) sim_stats(f, kap, k)))     # rows: kappa grid; columns: post, persist
        ok <- is.finite(S[, "persist"]) & is.finite(S[, "post"])
        if (sum(ok) < 2) next
        pk <- S[ok, "persist"]; po <- S[ok, "post"]; kk <- kg[ok]
        # the burden value at which the simulated persistence matches the observed one (linear interpolation;
        # clamped to the grid when the observed persistence lies outside the simulated range)
        if (persist_obs <= min(pk)) { kappa_hat[k] <- kk[which.min(pk)]; pred[k] <- po[which.min(pk)]; note[k] <- "persistence below simulated range" }
        else if (persist_obs >= max(pk)) { kappa_hat[k] <- kk[which.max(pk)]; pred[k] <- po[which.max(pk)]; note[k] <- "persistence above simulated range" }
        else {
          kappa_hat[k] <- stats::approx(pk, kk, xout = persist_obs, ties = "ordered")$y
          pred[k] <- stats::approx(kk, po, xout = kappa_hat[k], ties = "ordered")$y
        }
      }
    }
  }
  invert <- function(target) {
    ok <- is.finite(pred); o <- order(grid[ok]); g <- grid[ok][o]; pv <- pred[ok][o]
    if (length(g) < 2) return(c(NA_real_, 0))
    s <- sign(pv - target)
    cross <- which(s[-1] != s[-length(s)])
    if (length(cross) == 0) return(c(NA_real_, 0))
    j <- cross[which.min(abs(g[cross]))]     # the crossing nearest zero is used when there are several
    c(g[j] + (target - pv[j]) * (g[j + 1] - g[j]) / (pv[j + 1] - pv[j]), length(cross))
  }
  zq <- stats::qnorm(1 - (1 - level) / 2)
  iv <- invert(observed); dh <- iv[1]; ncross <- iv[2]
  ci <- c(invert(observed - zq * se)[1], invert(observed + zq * se)[1])
  # an endpoint that falls outside the range of the model-implied curve leaves that side open
  if (all(is.na(ci))) ci <- c(-Inf, Inf) else {
    ci <- sort(ci, na.last = TRUE)
    if (is.na(ci[2])) { ci <- if (!is.na(dh) && ci[1] > dh) c(-Inf, ci[1]) else c(ci[1], Inf) }
  }
  curve <- data.frame(delta = grid, predicted = pred)
  if (method == "postskip" && burden == "fit") curve$kappa_hat <- kappa_hat
  if (any(nzchar(note))) curve$note <- note
  out <- structure(list(method = method, delta_hat = dh, interval = ci, n_crossings = ncross, observed = observed, se = se,
                        curve = curve, variable = vars[wh], level = level, burden = if (method == "postskip") burden else NA,
                        persistence_observed = if (method == "postskip" && burden == "fit") persist_obs else NA_real_),
                   class = "delta_calibration")
  if (boot > 0 && !is.na(dh)) {
    # cluster bootstrap: resample persons (or days) with replacement, refit the profile on a
    # coarse grid around the point estimate and recalibrate; percentile interval
    if (boot_unit == "day" && is.null(day)) stop("boot_unit = 'day' needs a day column")
    if (is.null(boot_grid)) { step <- if (length(grid) > 1) min(diff(sort(grid))) else 0.25; boot_grid <- dh + step * (-3:3) }
    units <- if (boot_unit == "person") data[[id]] else data[[day]]
    u <- unique(units); dh_b <- rep(NA_real_, boot)
    with_seed(if (is.null(seed)) NULL else seed + 99L, {
      for (b in seq_len(boot)) {
        take <- sample(u, length(u), replace = TRUE)
        parts <- lapply(seq_along(take), function(k) { s <- data[units == take[k], , drop = FALSE]
          if (boot_unit == "person") s[[id]] <- paste0(s[[id]], "_", k) else s[[day]] <- k
          s })
        db <- do.call(rbind, parts)
        if (boot_unit == "day") { db <- db[order(db[[id]], db[[day]], db[[time]]), ]; db[[time]] <- seq_len(nrow(db)) }
        pb <- tryCatch(tilt_profile(db, vars, delta_grid = boot_grid, which = wh, id = id, time = time, day = day, R = R,
                                    propensity = profile$propensity, probe = if (method == "probe") probe else NULL), error = function(e) NULL)
        if (is.null(pb)) next
        cb <- tryCatch(calibrate_delta(pb, db, method = method, sensor = sensor, probe = probe, id = id, time = time, day = day, R = R,
                                       n_sim = n_sim, burden = burden, kappa_grid = kappa_grid,
                                       seed = if (is.null(seed)) NULL else seed + b, level = level, boot = 0L), error = function(e) NULL)
        if (!is.null(cb)) dh_b[b] <- cb$delta_hat
      }
    })
    out$boot <- dh_b
    out$boot_interval <- unname(stats::quantile(dh_b, c((1 - level) / 2, 1 - (1 - level) / 2), na.rm = TRUE))
    out$boot_se <- stats::sd(dh_b, na.rm = TRUE)
  }
  out
}

#' @export
print.delta_calibration <- function(x, digits = 3, ...) {
  cat(sprintf("Calibration of delta (%s) by the %s method%s\n", x$variable, x$method,
              if (identical(x$burden, "fit")) " (burden term calibrated to the observed response persistence)" else ""))
  cat(sprintf("  observed statistic: %.3f (SE %.3f)\n", x$observed, x$se))
  cat(sprintf("  calibrated delta: %.3f   %d%% interval: [%.3f, %.3f]\n", x$delta_hat, round(100 * x$level),
              x$interval[1], x$interval[2]))
  if (is.na(x$delta_hat)) cat("  (no crossing on the grid: the model-implied curve does not reach the observed value; widen delta_grid or report the profile band)\n")
  if (!is.na(x$n_crossings) && x$n_crossings > 1) cat(sprintf("  (%d crossings on the grid; the one nearest zero is reported)\n", x$n_crossings))
  if (!is.null(x$boot_interval)) cat(sprintf("  bootstrap (%d resamples%s): SE %.3f, percentile interval [%.3f, %.3f]\n",
                                             sum(!is.na(x$boot)), if (anyNA(x$boot)) sprintf(", %d failed or not identified", sum(is.na(x$boot))) else "",
                                             x$boot_se, x$boot_interval[1], x$boot_interval[2]))
  invisible(x)
}

#' Plot the model-implied curve of a calibration
#'
#' Draws the model-implied statistic against the sensitivity value (points
#' and line), the observed statistic (solid horizontal line) with its
#' confidence limits (dashed), and the calibrated value (dotted vertical
#' line).
#' @param x A \code{delta_calibration}.
#' @param ... Passed to \code{plot}.
#' @return Invisibly, the curve data frame (\code{x$curve} sorted by the
#'   sensitivity value). Called for its side effect, the plot.
#' @examples
#' sim <- simulate_ema(N = 40, T = 30, motifs = "M2", delta = -1, sensor_cor = 0.6, seed = 3)
#' prof <- tilt_profile(sim$data, sim$vars, delta_grid = c(-1.5, -1, -0.5, 0))
#' cs <- calibrate_delta(prof, sim$data, method = "sensor")
#' plot(cs, main = "Sensor calibration")
#' @export
plot.delta_calibration <- function(x, ...) {
  cv <- x$curve[order(x$curve$delta), ]
  zq <- stats::qnorm(1 - (1 - x$level) / 2)
  yl <- range(c(cv$predicted, x$observed + c(-1, 1) * zq * x$se), na.rm = TRUE)
  graphics::plot(cv$delta, cv$predicted, type = "b", pch = 16, ylim = yl,
                 xlab = expression(paste("sensitivity value ", delta)), ylab = paste("model-implied statistic (", x$method, ")"), ...)
  graphics::abline(h = x$observed, lwd = 2); graphics::abline(h = x$observed + c(-1, 1) * zq * x$se, lty = 2)
  if (!is.na(x$delta_hat)) graphics::abline(v = x$delta_hat, lty = 3)
  invisible(cv)
}

#' Simulate data from a fitted tilt model
#'
#' Generates a data set from the fitted VAR(1) with the fitted person means
#' and (when the fit has person-specific intercepts) response intercepts,
#' both resampled jointly with replacement so that the association between a
#' person's level and propensity is preserved, and the probit self-censoring
#' model at the fit's sensitivity value and calibrated intercept. An optional
#' burden term lowers the propensity after a skipped prompt.
#' @param fit A \code{tilt_fit} (or a plain \code{pairs_fit}, which is treated
#'   as a fit at \code{delta = 0}).
#' @param N,T Persons and prompts.
#' @param days Optional number of prompts per day; adds a \code{day} column and
#'   restarts the burden term at the first prompt of each day.
#' @param kappa_R Burden term: probit drop in the response propensity at a
#'   prompt whose predecessor (within the day) was skipped.
#' @param seed Optional seed; the caller's random-number state is restored.
#' @return A long data frame with \code{id}, \code{time}, optional \code{day},
#'   \code{R} (1 = answered) and the state columns (\code{NA} at skipped
#'   prompts), in the format accepted by every estimator of the package.
#' @examples
#' sim <- simulate_ema(N = 30, T = 20, motifs = "M2", seed = 1)
#' f <- fit_tilt(sim$data, sim$vars, delta = -1)
#' d <- simulate_from_fit(f, N = 30, T = 20, seed = 2)
#' mean(d$R)
#' head(d)
#' # with a day structure and a burden term
#' d2 <- simulate_from_fit(f, N = 30, T = 20, days = 5, kappa_R = 1, seed = 2)
#' fatigue_check(d2, day = "day")
#' @export
simulate_from_fit <- function(fit, N, T, days = NULL, kappa_R = 0, seed = NULL) {
  with_seed(seed, {
    vars <- fit$vars; p <- length(vars)
    Phi <- fit$Phi; Psi <- fit$Psi
    if (max(Mod(eigen(Phi, only.values = TRUE)$values)) >= 1)
      stop("the fitted lagged-coefficient matrix is not stable (spectral radius >= 1); no stationary process exists at this sensitivity value")
    delta <- if (is.null(fit$delta)) rep(0, p) else fit$delta
    idx <- sample(nrow(fit$mu_dyn), N, replace = TRUE)     # resample persons jointly: means and propensities stay paired
    mu_i <- fit$mu_dyn[idx, , drop = FALSE]
    # at delta = 0 (an untilted fit) the response model is a constant probit whose intercept reproduces the response rate
    alpha <- if (is.null(fit$alpha) || all(is.na(fit$alpha))) rep(stats::qnorm(fit$pairs$response_rate), N) else if (length(fit$alpha) == 1L) rep(fit$alpha, N) else fit$alpha[idx]
    if (anyNA(alpha)) alpha[is.na(alpha)] <- mean(alpha, na.rm = TRUE)
    Sinf <- stationary_cov(Phi, Psi); Ls <- t(chol(Sinf)); Lp <- t(chol(Psi))
    x <- mu_i + t(Ls %*% matrix(stats::rnorm(N * p), p, N))
    out <- matrix(NA_real_, N * T, p); Rv <- integer(N * T)
    r_prev <- rep(1L, N)
    for (tt in seq_len(T)) {
      x <- mu_i + (x - mu_i) %*% t(Phi) + t(Lp %*% matrix(stats::rnorm(N * p), p, N))
      first_of_day <- !is.null(days) && ((tt - 1) %% days == 0)
      bt <- if (kappa_R != 0 && !first_of_day) -kappa_R * (1L - r_prev) else 0
      r <- as.integer(stats::runif(N) < stats::pnorm(alpha + as.vector(x %*% delta) + bt))
      rows <- (seq_len(N) - 1) * T + tt
      out[rows, ] <- x; Rv[rows] <- r; r_prev <- r
    }
    out[Rv == 0L, ] <- NA
    colnames(out) <- vars
    d <- data.frame(id = rep(seq_len(N), each = T), time = rep(seq_len(T), N))
    if (!is.null(days)) d$day <- (d$time - 1) %/% days + 1
    d$R <- Rv
    cbind(d, as.data.frame(out))
  })
}
