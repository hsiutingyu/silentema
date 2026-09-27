#' Inverse-probability-weighted within-person VAR(1) for observed context confounding
#'
#' For motif M5 with an observed context variable (\code{C_t -> X_t} and
#' \code{C_t -> R_t}), complete adjacent pairs are reweighted by the
#' stabilized weight \eqn{w_t = P(R_t = 1 \mid x_{t-1}) / P(R_t = 1 \mid c_t, x_{t-1})}
#' so that the context-marginal transition kernel is recovered (the observed-
#' context case of the recoverability results in Yu, 2026). Both response models are probit regressions fitted on all
#' prompts whose predecessor was answered (the response indicator is always
#' observed), with person-specific intercepts absorbed by including the
#' person's observed response rate as an offset-like covariate when
#' \code{propensity = "person"}.
#'
#' @inheritParams fit_pairs
#' @param context Name(s) of the observed context column(s).
#' @param propensity \code{"common"} or \code{"person"} (adds the person's
#'   response rate on other prompts as a covariate of the response model).
#' @return An object of class \code{c("ipw_fit", "pairs_fit")}: a
#'   \code{\link{fit_pairs}} object (weighted) with the additional elements
#'   \code{response_model} (the fitted denominator probit model, a
#'   \code{glm}) and \code{ess} (the effective number of pairs); the weights
#'   are in \code{weights}. The weights are exact only when the context is serially
#'   independent; under a persistent context the pair selection through
#'   \eqn{R_{t-1}} depends on \eqn{C_{t-1}}, and covariate adjustment
#'   (\code{fit_pairs(covariates = )}) is preferable.
#' @seealso \code{\link{fit_pairs}} with \code{covariates} for covariate
#'   adjustment, which is preferable under a persistent context.
#' @examples
#' sim <- simulate_ema(N = 40, T = 30, motifs = "M5", seed = 1)   # the context is recorded in column C
#' f <- fit_ipw(sim$data, sim$vars, context = "C")
#' f
#' f$ess / f$n_pairs
#' summary(f$response_model)$coefficients
#' @export
fit_ipw <- function(data, vars, context = "C", id = "id", time = "time", day = NULL, R = "R",
                    propensity = c("common", "person"), pairs = NULL) {
  propensity <- match.arg(propensity)
  check_columns(data, context)
  pr <- if (is.null(pairs)) make_pairs(data, vars, id, time, day, R) else pairs
  lo <- pr$lag_obs
  Cl <- as.matrix(data[lo$row_t, context, drop = FALSE])
  Xl <- lo$x
  df_lo <- data.frame(Rr = lo$R, Xl, Cl)
  names(df_lo) <- c("Rr", paste0("x", seq_len(ncol(Xl))), paste0("c", seq_len(ncol(Cl))))
  if (propensity == "person") {
    # leave-one-out person response rate over all prompts
    Rv <- pr$Rv; ids <- data[[id]]
    tot <- tapply(Rv, ids, sum); n_i <- tapply(Rv, ids, length)
    loo <- (tot[as.character(lo$id)] - lo$R) / pmax(n_i[as.character(lo$id)] - 1, 1)
    df_lo$prop <- as.numeric(loo)
  }
  xn <- paste0("x", seq_len(ncol(Xl))); cn <- paste0("c", seq_len(ncol(Cl)))
  f_num <- stats::as.formula(paste("Rr ~", paste(c(xn, if (propensity == "person") "prop"), collapse = " + ")))
  f_den <- stats::as.formula(paste("Rr ~", paste(c(xn, cn, if (propensity == "person") "prop"), collapse = " + ")))
  m_num <- stats::glm(f_num, family = stats::binomial("probit"), data = df_lo)
  m_den <- stats::glm(f_den, family = stats::binomial("probit"), data = df_lo)
  # evaluate at the complete pairs
  Cp <- as.matrix(data[pr$row_t, context, drop = FALSE])
  df_p <- data.frame(pr$x, Cp); names(df_p) <- c(xn, cn)
  if (propensity == "person") df_p$prop <- as.numeric((tot[as.character(pr$id)] - 1) / pmax(n_i[as.character(pr$id)] - 1, 1))
  num <- stats::predict(m_num, newdata = df_p, type = "response")
  den <- stats::predict(m_den, newdata = df_p, type = "response")
  w <- num / pmax(den, 1e-6)
  fit <- fit_pairs(data, vars, id, time, day, R, weights = w, pairs = pr)
  fit$response_model <- m_den
  fit$ess <- sum(w)^2 / sum(w^2)
  class(fit) <- c("ipw_fit", "pairs_fit")
  fit
}
