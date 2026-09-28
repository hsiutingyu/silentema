#' Stationary covariance of a VAR(1) process
#'
#' Solves \eqn{\Sigma = \Phi \Sigma \Phi' + \Psi} by vectorization.
#' @param Phi Lagged-coefficient matrix.
#' @param Psi Innovation covariance matrix.
#' @return The stationary covariance matrix (same dimension as \code{Phi}).
#'   An error is raised when \code{Phi} is not stable.
#' @examples
#' p <- default_params()
#' round(diag(stationary_cov(p$Phi, p$Psi)), 3)
#' S <- stationary_cov(p$Phi, p$Psi)
#' all.equal(S, p$Phi %*% S %*% t(p$Phi) + p$Psi)
#' @export
stationary_cov <- function(Phi, Psi) {
  p <- nrow(Phi)
  if (max(Mod(eigen(Phi, only.values = TRUE)$values)) >= 1)
    stop("'Phi' is not stable (spectral radius >= 1): no stationary covariance exists")
  A <- diag(p * p) - kronecker(Phi, Phi)
  matrix(solve(A, as.vector(Psi)), p, p)
}

#' Partial correlations from a covariance matrix
#'
#' Inverts the covariance matrix (or takes the precision matrix as given) and
#' standardizes its negative off-diagonal elements. Applied to the innovation
#' covariance of the VAR(1), the result is the contemporaneous network: the
#' partial correlations of the states at the same prompt given the previous
#' prompt and the other states.
#' @param S A covariance (or precision, with \code{precision = TRUE}) matrix.
#' @param precision Logical; is \code{S} already a precision matrix?
#' @return The matrix of partial correlations (the contemporaneous network
#'   when \code{S} is an innovation covariance) with unit diagonal.
#' @examples
#' round(partial_cors(default_params()$Psi), 2)
#' # from a precision matrix directly
#' K <- solve(default_params()$Psi)
#' all.equal(partial_cors(K, precision = TRUE), partial_cors(default_params()$Psi))
#' @export
partial_cors <- function(S, precision = FALSE) {
  K <- if (precision) S else solve(S)
  d <- 1 / sqrt(diag(K))
  P <- -K * outer(d, d)
  diag(P) <- 1
  P
}

#' Default population parameters used in the simulation studies of Yu (2026)
#'
#' A four-variable within-person system with autoregressive coefficients
#' between .30 and .45, a few cross-lagged effects, an innovation covariance
#' with a sparse contemporaneous network, and stationary variances equal to
#' one.
#' @return A list with \code{Phi} (4 x 4 lagged coefficients), \code{Psi}
#'   (innovation covariance, scaled so that the stationary variances equal
#'   one), \code{mu} (population means), \code{Sigma_mu} (between-person
#'   covariance of the person means) and \code{vars}, the variable names
#'   \code{c("NegA", "PosA", "Stress", "Fatigue")}.
#' @examples
#' p <- default_params()
#' p$Phi
#' # contemporaneous network: NegA--PosA -.30, NegA--Stress .30, Stress--Fatigue .20
#' round(partial_cors(p$Psi), 2)
#' @export
default_params <- function() {
  Phi <- matrix(c(0.40, 0.00, 0.15, 0.00,
                  -0.10, 0.35, 0.00, 0.00,
                  0.15, 0.00, 0.45, 0.10,
                  0.00, -0.10, 0.00, 0.30), 4, 4, byrow = TRUE)
  # innovation precision with a sparse contemporaneous network
  K <- matrix(c(1.0, 0.30, -0.30, 0.00,
                0.30, 1.0, 0.00, 0.00,
                -0.30, 0.00, 1.0, -0.20,
                0.00, 0.00, -0.20, 1.0), 4, 4, byrow = TRUE)
  Psi_raw <- solve(K)
  # scale so that stationary variances equal one
  S_raw <- stationary_cov(Phi, Psi_raw)
  D <- diag(1 / sqrt(diag(S_raw)))
  Psi <- D %*% Psi_raw %*% D
  # note: rescaling Psi by D does not preserve stationary unit variance exactly
  # because Phi is not diagonal; iterate a few times
  for (k in 1:25) {
    S <- stationary_cov(Phi, Psi)
    D <- diag(1 / sqrt(diag(S)))
    Psi <- D %*% Psi %*% D
  }
  Sigma_mu <- matrix(c(0.5, 0.15, 0.20, 0.10,
                       0.15, 0.5, 0.00, 0.00,
                       0.20, 0.00, 0.5, 0.10,
                       0.10, 0.00, 0.10, 0.5), 4, 4, byrow = TRUE)
  vn <- c("NegA", "PosA", "Stress", "Fatigue")
  dimnames(Phi) <- dimnames(Psi) <- dimnames(Sigma_mu) <- list(vn, vn)
  list(Phi = Phi, Psi = Psi, mu = stats::setNames(c(2.5, 4.0, 2.8, 3.0), vn), Sigma_mu = Sigma_mu, vars = vn)
}

#' Simulate experience-sampling data with a declared missingness mechanism
#'
#' Generates a two-level VAR(1) process for \code{N} persons over \code{n_prompts}
#' prompts and then deletes whole prompts according to one or more
#' missingness motifs. Response probabilities use a probit link. Intercepts
#' are calibrated numerically so that the realized response rate matches
#' \code{compliance}.
#'
#' @param N Number of persons.
#' @param n_prompts Number of prompts per person (called \code{T} in the
#'   archived versions 0.2.x of the package).
#' @param params List with \code{Phi}, \code{Psi}, \code{mu}, \code{Sigma_mu}
#'   (see \code{\link{default_params}}).
#' @param motifs Character vector of motif codes (see \code{\link{dm_graph}}).
#' @param compliance Target overall response rate.
#' @param gamma Coefficient vector of \code{X_{t-1}} in the response model
#'   (motif M1); a scalar is applied to the first variable.
#' @param delta Coefficient vector of \code{X_t} in the response model
#'   (motif M2); a scalar is applied to the first variable. Negative values
#'   mean that high states are skipped. Both \code{gamma} and \code{delta}
#'   act on the raw (uncentered) state, so a self-censoring person with a high
#'   typical state also responds less often overall.
#' @param kappa_R Coefficient of \code{R_{t-1}} in the response model (M3);
#'   the term enters as \code{kappa_R * (R_{t-1} - compliance)}, so the
#'   propensity after an answered prompt exceeds that after a skipped prompt
#'   by \code{kappa_R} probit units.
#' @param burden Coefficient of a linear time trend in the response model
#'   (M3); the term enters as \code{burden * trend}, where the trend runs from
#'   -0.5 to 0.5 over the burn-in and the recorded prompts together, so a
#'   negative value gives declining compliance over the study and the realized
#'   change over the recorded prompts is about
#'   \code{burden * n_prompts / (n_prompts + burnin)}
#'   probit units.
#' @param rho_propensity Correlation between the person-level response
#'   propensity and the person's mean on the first variable (M4).
#' @param sd_propensity Standard deviation of the person-level propensity
#'   (M4), in probit units.
#' @param p_context Probability that the (binary) context \code{C_t = 1}
#'   (M5).
#' @param rho_context Persistence of the context: with this probability
#'   \code{C_t} copies \code{C_{t-1}}, otherwise it is drawn afresh with
#'   probability \code{p_context}; \code{0} gives serially independent
#'   contexts (the marginal rate is \code{p_context} either way).
#' @param gamma_C Shift of the state vector when \code{C_t = 1} (M5); a scalar
#'   is applied to the first variable.
#' @param kappa_C Coefficient of \code{C_t} in the response model (M5).
#' @param rho_react Shift of the state vector at prompt \code{t} when the
#'   previous prompt was answered (M6, reactivity); the term enters as
#'   \code{rho_react * (R_{t-1} - compliance)}, so the contrast between an
#'   answered and a skipped predecessor is \code{rho_react}; a scalar applies
#'   to the first variable.
#' @param sensor_cor Correlation between an always-observed sensor channel
#'   and the within-person fluctuation of the first state variable (its
#'   deviation from the person mean, standardized); \code{NULL} for no sensor.
#' @param p_probe Probability that a prompt is a randomized probe that forces a
#'   response; \code{0} for no probes.
#' @param burnin Number of burn-in prompts discarded before recording.
#' @param seed Optional integer seed (the caller's random-number state is
#'   restored afterwards).
#' @param days Optional number of prompts per day; when given, a \code{day}
#'   column is added (the response model itself does not use the day).
#'
#' @return A list of class \code{"ema_sim"} with \code{data} (long data frame:
#'   \code{id}, \code{time}, optional \code{day}, \code{R}, the state columns
#'   with \code{NA} at skipped prompts, and the optional columns \code{C}
#'   (context, M5), \code{S} (sensor) and \code{Z} (probe)), \code{full}
#'   (the state columns without deletion, with \code{id} and \code{time}),
#'   \code{mu_i} (the true person means), \code{alpha0} (the calibrated
#'   response intercept), \code{a_i} (the person propensities, M4),
#'   \code{params}, \code{motifs}, \code{settings} (the design and the
#'   response-model arguments) and \code{vars}. It has a \code{print} method.
#' @seealso \code{\link{dm_graph}} for the same motifs as a graph,
#'   \code{\link{simulate_from_fit}} to simulate from a fitted model,
#'   \code{\link{default_params}}.
#' @examples
#' sim <- simulate_ema(N = 20, n_prompts = 15, motifs = c("M2", "M4"), compliance = 0.7, seed = 1)
#' sim
#' head(sim$data)
#' # the deleted states are kept for checking estimators
#' head(sim$full)
#' # a sensor, probes, a persistent observed context and a day structure
#' sim2 <- simulate_ema(N = 20, n_prompts = 20, motifs = c("M2", "M5"), sensor_cor = 0.6,
#'                      p_probe = 0.1, rho_context = 0.5, days = 5, seed = 2)
#' names(sim2$data)
#' mean(sim2$data$R[sim2$data$Z == 1])   # probes are always answered
#' @export
simulate_ema <- function(N = 100, n_prompts = 56, params = default_params(), motifs = "M0",
                         compliance = 0.75, gamma = -0.5, delta = -1.0, kappa_R = 1.0,
                         burden = 0, rho_propensity = 0.4, sd_propensity = 0.5,
                         p_context = 0.3, rho_context = 0, gamma_C = 0.6, kappa_C = 1.0, rho_react = 0.2,
                         sensor_cor = NULL, p_probe = 0, burnin = 50, seed = NULL, days = NULL) {
  motifs <- toupper(motifs)
  bad <- setdiff(motifs, c("M0", "M1", "M2", "M3", "M4", "M5", "M6"))
  if (length(bad)) stop("Unknown motif code(s): ", paste(bad, collapse = ", "))
  if (!is.numeric(N) || !is.numeric(n_prompts) || is.logical(N) || is.logical(n_prompts) || N < 1 || n_prompts < 3)
    stop("'N' must be at least 1 and 'n_prompts' at least 3 (numeric)")
  if (compliance <= 0 || compliance > 1) stop("'compliance' must lie in (0, 1]")
  with_seed(seed, {
  Phi <- params$Phi; Psi <- params$Psi; mu <- params$mu; Sigma_mu <- params$Sigma_mu
  p <- nrow(Phi); vn <- if (!is.null(params$vars)) params$vars else paste0("X", seq_len(p))
  expand <- function(v) { if (length(v) == 1L) { out <- rep(0, p); out[1] <- v; out } else v }
  gamma <- expand(gamma); delta <- expand(delta); gamma_C <- expand(gamma_C); rho_react <- expand(rho_react)
  Ltmu <- t(chol(Sigma_mu)); Lpsi <- t(chol(Psi))
  Tt <- n_prompts + burnin
  # person means
  mu_i <- t(mu + Ltmu %*% matrix(stats::rnorm(N * p), p, N))   # N x p
  # context
  C <- NULL
  if ("M5" %in% motifs) {
    C <- matrix(stats::rbinom(N * Tt, 1, p_context), N, Tt)
    if (rho_context > 0) for (tt in 2:Tt) {
      copy <- stats::runif(N) < rho_context
      C[copy, tt] <- C[copy, tt - 1]
    }
  }
  # person propensity (M4): probit intercept deviation correlated with mu_i[,1]
  a_i <- rep(0, N)
  if ("M4" %in% motifs) {
    z1 <- (mu_i[, 1] - mu[1]) / sqrt(Sigma_mu[1, 1])
    a_i <- sd_propensity * (rho_propensity * z1 + sqrt(1 - rho_propensity^2) * stats::rnorm(N))
  }
  # generate states and responses jointly (M6 needs R_{t-1})
  Sinf <- stationary_cov(Phi, Psi); Ls <- t(chol(Sinf))
  x_prev <- mu_i + t(Ls %*% matrix(stats::rnorm(N * p), p, N))
  r_prev <- stats::rbinom(N, 1, compliance)
  eps <- array(stats::rnorm(N * Tt * p), c(N, Tt, p))
  trend <- seq(-0.5, 0.5, length.out = Tt)
  # first pass: generate the latent process; responses need intercept calibration,
  # which depends on the realized states. Because M6 feeds R back into X, we
  # calibrate the intercept on a pilot pass and then generate for real.
  gen <- function(alpha0) {   # returns the states and responses generated with response intercept alpha0
    Xg <- array(NA_real_, c(N, Tt, p)); Rg <- matrix(NA_integer_, N, Tt)
    xp <- x_prev; rp <- r_prev
    for (tt in seq_len(Tt)) {
      e <- t(Lpsi %*% t(matrix(eps[, tt, ], N, p)))
      xt <- mu_i + (xp - mu_i) %*% t(Phi) + e
      if (!is.null(C)) xt <- xt + outer(C[, tt], gamma_C)
      if ("M6" %in% motifs) xt <- xt + outer(rp - compliance, rho_react)
      lin <- rep(alpha0, N) + a_i
      if ("M1" %in% motifs) lin <- lin + as.vector(xp %*% gamma)
      if ("M2" %in% motifs) lin <- lin + as.vector(xt %*% delta)
      if ("M3" %in% motifs) lin <- lin + kappa_R * (rp - compliance) + burden * trend[tt]
      if (!is.null(C)) lin <- lin - kappa_C * (C[, tt] - p_context)
      u <- U[, tt]
      rt <- as.integer(u < stats::pnorm(lin))
      if (p_probe > 0) rt[Zm[, tt] == 1L] <- 1L
      Xg[, tt, ] <- xt; Rg[, tt] <- rt
      xp <- xt; rp <- rt
    }
    list(X = Xg, Rm = Rg)
  }
  U <- matrix(stats::runif(N * Tt), N, Tt)
  Zm <- if (p_probe > 0) matrix(stats::rbinom(N * Tt, 1, p_probe), N, Tt) else NULL
  # calibrate alpha0 so that the realized response rate over the recorded prompts equals compliance
  rate_at <- function(alpha0) mean(gen(alpha0)$Rm[, (burnin + 1):Tt]) - compliance
  alpha0 <- tryCatch(stats::uniroot(rate_at, c(-30, 30), tol = 1e-4)$root,
                     error = function(e) stop("Could not calibrate the response intercept to the target compliance: ", conditionMessage(e)))
  g <- gen(alpha0); X <- g$X; Rm <- g$Rm
  keep <- (burnin + 1):Tt
  Xk <- X[, keep, , drop = FALSE]; Rk <- Rm[, keep, drop = FALSE]
  Ck <- if (!is.null(C)) C[, keep, drop = FALSE] else NULL
  Zk <- if (!is.null(Zm)) Zm[, keep, drop = FALSE] else NULL
  # sensor
  Sk <- NULL
  if (!is.null(sensor_cor)) {
    sdx <- sqrt(Sinf[1, 1])
    Sk <- matrix(NA_real_, N, n_prompts)
    for (tt in seq_len(n_prompts)) Sk[, tt] <- sensor_cor * (Xk[, tt, 1] - mu_i[, 1]) / sdx + sqrt(1 - sensor_cor^2) * stats::rnorm(N)
  }
  # long data frame
  id <- rep(seq_len(N), each = n_prompts); time <- rep(seq_len(n_prompts), times = N)
  full <- matrix(NA_real_, N * n_prompts, p, dimnames = list(NULL, vn))
  for (j in seq_len(p)) full[, j] <- as.vector(t(Xk[, , j]))
  Rv <- as.vector(t(Rk))
  obs <- full; obs[Rv == 0L, ] <- NA
  d <- data.frame(id = id, time = time)
  if (!is.null(days)) d$day <- (time - 1) %/% days + 1
  d$R <- Rv
  d <- cbind(d, as.data.frame(obs))
  if (!is.null(Ck)) d$C <- as.vector(t(Ck))
  if (!is.null(Sk)) d$S <- as.vector(t(Sk))
  if (!is.null(Zk)) d$Z <- as.vector(t(Zk))
  structure(list(data = d, full = cbind(id = id, time = time, as.data.frame(full)), mu_i = mu_i,
                 alpha0 = alpha0, a_i = a_i, params = params, motifs = motifs,
                 settings = list(N = N, n_prompts = n_prompts, compliance = compliance, gamma = gamma, delta = delta,
                                 kappa_R = kappa_R, burden = burden, rho_propensity = rho_propensity,
                                 sd_propensity = sd_propensity, p_context = p_context, rho_context = rho_context, gamma_C = gamma_C,
                                 kappa_C = kappa_C, rho_react = rho_react, sensor_cor = sensor_cor,
                                 p_probe = p_probe),
                 vars = vn), class = "ema_sim")
  })
}

#' @export
print.ema_sim <- function(x, ...) {
  cat("Simulated EMA data:", x$settings$N, "persons x", x$settings$n_prompts, "prompts,", length(x$vars), "variables\n")
  cat("  motifs:", paste(x$motifs, collapse = " + "), "  realized response rate:", round(mean(x$data$R), 3), "\n")
  invisible(x)
}
