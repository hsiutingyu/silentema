#' Full-information maximum likelihood under missing at random (state-space EM)
#'
#' Fits the two-level VAR(1) with random person means by maximum likelihood,
#' treating skipped prompts as missing at random. The likelihood is computed
#' by the Kalman filter over the augmented state \eqn{(X_t, \mu_i)} and
#' maximized by the EM algorithm of Shumway and Stoffer (1982). This is the
#' likelihood that dynamic structural equation modeling and related
#' software maximize (with a diffuse or informative prior instead of the
#' flat likelihood), and it is the "default pipeline" against which Yu (2026)
#' measures the damage of state-dependent nonresponse.
#'
#' @inheritParams make_pairs
#' @param start Optional list with \code{Phi}, \code{Psi}, \code{mu},
#'   \code{Sigma_mu}; by default the answered-pairs estimates are used.
#' @param max_iter,tol EM control.
#' @param verbose Print the EM trace.
#' @details The prompts of a person are treated as one equally spaced series
#'   (the overnight transition is a lag-1 transition, as in the default
#'   specification of dynamic structural equation models with equally spaced
#'   occasions); there is no \code{day} argument. Standard errors are not
#'   computed (the estimator serves as the MAR benchmark in the simulation
#'   studies of Yu, 2026). \code{loglik} holds the log-likelihood evaluated at the
#'   parameters entering each EM iteration, so its last element belongs to the
#'   penultimate iterate; \code{\link{loglik_fiml}} evaluates the returned
#'   estimates.
#' @return An object of class \code{"fiml_fit"}: a list with \code{Phi},
#'   \code{Psi}, \code{pcor}, \code{mu}, \code{Sigma_mu}, \code{Sigma0}
#'   (the stationary covariance of the within-person deviations), \code{loglik}
#'   (the trace of log-likelihood values over the EM iterations),
#'   \code{iterations}, \code{converged}, \code{N}, \code{n_prompts} (the
#'   number of prompts spanned by the prompt index) and \code{vars}. It has a \code{print} method.
#' @examples
#' sim <- simulate_ema(N = 30, n_prompts = 20, motifs = "M1", seed = 1)
#' f <- fit_fiml(sim$data, sim$vars)
#' f
#' f$converged; length(f$loglik)
#' # under self-censoring the MAR likelihood is biased like the answered-pairs estimator
#' sim2 <- simulate_ema(N = 30, n_prompts = 20, motifs = "M2", delta = -1, seed = 1)
#' c(fiml = fit_fiml(sim2$data, sim2$vars)$Phi[1, 1],
#'   pairs = fit_pairs(sim2$data, sim2$vars)$Phi[1, 1],
#'   truth = default_params()$Phi[1, 1])
#' @references Shumway, R. H., & Stoffer, D. S. (1982). An approach to time
#'   series smoothing and forecasting using the EM algorithm. \emph{Journal of
#'   Time Series Analysis, 3}, 253-264. \doi{10.1111/j.1467-9892.1982.tb00349.x}
#' @export
fit_fiml <- function(data, vars, id = "id", time = "time", R = "R", start = NULL,
                     max_iter = 200L, tol = 1e-5, verbose = FALSE) {
  p <- length(vars)
  check_columns(data, c(vars, id, time))
  check_time(data, time)
  Rv <- response_indicator(data, vars, R)
  ord <- order(data[[id]], data[[time]]); data <- data[ord, , drop = FALSE]; Rv <- Rv[ord]
  ids <- unique(data[[id]]); N <- length(ids)
  Tn <- max(data[[time]]) - min(data[[time]]) + 1L
  if (Tn < 2L) stop("the prompt index must span at least two prompts")
  if (nrow(data) < N * Tn) warning("some persons have no row for some prompt indices between the first and the last; those prompts are treated as skipped")
  t0 <- min(data[[time]])
  Y <- array(NaN, c(p, Tn, N)); Rm <- matrix(0L, Tn, N)
  Xm <- as.matrix(data[, vars, drop = FALSE])
  ii <- match(data[[id]], ids); tt <- data[[time]] - t0 + 1L
  for (r in seq_len(nrow(data))) if (Rv[r] == 1L) { Y[, tt[r], ii[r]] <- Xm[r, ]; Rm[tt[r], ii[r]] <- 1L }
  if (is.null(start)) {
    f0 <- tryCatch(suppressWarnings(fit_pairs(data, vars, id, time, NULL, R)), error = function(e) NULL)   # starting values only
    if (is.null(f0)) start <- list(Phi = diag(0.3, p), Psi = diag(stats::var(Xm, na.rm = TRUE)) * diag(p) * 0.7,
                                   mu = colMeans(Xm, na.rm = TRUE), Sigma_mu = diag(0.3, p))
    else {
      # the covariance of the observed person means is always positive definite; the jackknifed
      # between-person covariance need not be, and a singular start is an absorbing state of EM
      Smu0 <- if (!is.null(f0$Sigma_mu_obs) && all(is.finite(f0$Sigma_mu_obs))) f0$Sigma_mu_obs else
        if (all(is.finite(f0$Sigma_mu_ls))) f0$Sigma_mu_ls else diag(0.3 * mean(diag(f0$Psi)), p)   # few persons: no between-person covariance yet
      mu0 <- if (all(is.finite(f0$mu_obs_person))) f0$mu_obs_person else colMeans(Xm, na.rm = TRUE)
      start <- list(Phi = f0$Phi, Psi = f0$Psi, mu = mu0, Sigma_mu = Smu0)
    }
  }
  # regularise the start: eigenvalues of Sigma_mu bounded below by 5% of its average variance
  ev <- eigen(start$Sigma_mu, symmetric = TRUE)
  floor_ev <- 0.05 * mean(diag(start$Sigma_mu))
  if (min(ev$values) < floor_ev) start$Sigma_mu <- ev$vectors %*% diag(pmax(ev$values, floor_ev), p) %*% t(ev$vectors)
  ev <- eigen(start$Psi, symmetric = TRUE)
  if (min(ev$values) < 1e-3) start$Psi <- ev$vectors %*% diag(pmax(ev$values, 1e-3), p) %*% t(ev$vectors)
  if (max(Mod(eigen(start$Phi)$values)) >= 0.95) start$Phi <- start$Phi * 0.9 / max(Mod(eigen(start$Phi)$values))
  Sigma0 <- stationary_cov(start$Phi, start$Psi)
  res <- kalman_em_cpp(Y, Rm, start$Phi, start$Psi, start$mu, start$Sigma_mu, Sigma0,
                       max_iter = as.integer(max_iter), tol = tol, verbose = verbose)
  Phi <- res$Phi; Psi <- res$Psi
  dimnames(Phi) <- dimnames(Psi) <- dimnames(res$Sigma_mu) <- list(vars, vars)
  structure(list(Phi = Phi, Psi = Psi, pcor = partial_cors(Psi), mu = stats::setNames(as.vector(res$mu), vars),
                 Sigma_mu = res$Sigma_mu, Sigma0 = res$Sigma0, loglik = res$loglik,
                 iterations = res$iterations, converged = res$converged, N = N, n_prompts = Tn, vars = vars),
            class = "fiml_fit")
}

#' @export
print.fiml_fit <- function(x, digits = 3, ...) {
  cat("FIML (MAR) two-level VAR(1) by state-space EM:", x$N, "persons,", x$n_prompts, "prompts;",
      x$iterations, "iterations,", if (x$converged) "converged" else "NOT converged", "\n")
  cat("Lagged coefficients Phi:\n"); print(round(x$Phi, digits))
  cat("Contemporaneous partial correlations:\n"); print(round(x$pcor, digits))
  cat("Means:", round(x$mu, digits), "\n")
  invisible(x)
}

#' Log-likelihood of the two-level VAR(1) at given parameters (MAR)
#'
#' Evaluates, without fitting, the likelihood that \code{\link{fit_fiml}}
#' maximizes: the Gaussian likelihood of the answered states under the
#' two-level VAR(1) with random person means, skipped prompts integrated out
#' as missing at random. Useful for comparing fits, for profile likelihoods,
#' and for checking a fit against the trace of the EM iterations.
#' @inheritParams fit_fiml
#' @param Phi,Psi,mu,Sigma_mu,Sigma0 Parameter values: the \eqn{p \times p}
#'   lagged-coefficient matrix, the \eqn{p \times p} innovation covariance,
#'   the length-\eqn{p} population mean, the \eqn{p \times p} between-person
#'   covariance of the person means, and the \eqn{p \times p} covariance of
#'   the within-person deviation at the first prompt (by default the stationary
#'   covariance implied by \code{Phi} and \code{Psi}).
#' @return A single number, the log-likelihood of the observed states under
#'   the two-level VAR(1) with the given parameters, skipped prompts treated as
#'   missing at random.
#' @examples
#' sim <- simulate_ema(N = 20, n_prompts = 15, seed = 1)
#' f <- fit_fiml(sim$data, sim$vars)
#' loglik_fiml(sim$data, sim$vars, f$Phi, f$Psi, f$mu, f$Sigma_mu)
#' # the likelihood at the population parameters is lower than at the ML estimates
#' p <- default_params()
#' loglik_fiml(sim$data, sim$vars, p$Phi, p$Psi, p$mu, p$Sigma_mu)
#' @export
loglik_fiml <- function(data, vars, Phi, Psi, mu, Sigma_mu, Sigma0 = stationary_cov(Phi, Psi),
                        id = "id", time = "time", R = "R") {
  p <- length(vars)
  check_columns(data, c(vars, id, time))
  check_time(data, time)
  Rv <- response_indicator(data, vars, R)
  ord <- order(data[[id]], data[[time]]); data <- data[ord, , drop = FALSE]; Rv <- Rv[ord]
  ids <- unique(data[[id]]); N <- length(ids)
  Tn <- max(data[[time]]) - min(data[[time]]) + 1L; t0 <- min(data[[time]])
  Y <- array(NaN, c(p, Tn, N)); Rm <- matrix(0L, Tn, N)
  Xm <- as.matrix(data[, vars, drop = FALSE]); ii <- match(data[[id]], ids); tt <- data[[time]] - t0 + 1L
  for (r in seq_len(nrow(data))) if (Rv[r] == 1L) { Y[, tt[r], ii[r]] <- Xm[r, ]; Rm[tt[r], ii[r]] <- 1L }
  kalman_loglik_cpp(Y, Rm, Phi, Psi, mu, Sigma_mu, Sigma0)
}
