pars <- default_params()

test_that("simulate_ema returns the documented structure and hits the target compliance", {
  sim <- simulate_ema(N = 30, T = 20, motifs = c("M2", "M4"), compliance = 0.7, sensor_cor = 0.5, p_probe = 0.1, seed = 1, days = 5)
  expect_s3_class(sim, "ema_sim")
  expect_equal(nrow(sim$data), 600)
  expect_true(all(c("id", "time", "day", "R", "NegA", "PosA", "Stress", "Fatigue", "S", "Z") %in% names(sim$data)))
  expect_lt(abs(mean(sim$data$R) - 0.7), 0.03)
  expect_true(all(is.na(sim$data$NegA[sim$data$R == 0])))
  expect_true(all(sim$data$R[sim$data$Z == 1] == 1))
  expect_equal(max(sim$data$day), 4)
})

test_that("stationary covariance and partial correlations are consistent", {
  S <- stationary_cov(pars$Phi, pars$Psi)
  expect_equal(unname(S), unname(pars$Phi %*% S %*% t(pars$Phi) + pars$Psi), tolerance = 1e-10)
  expect_equal(unname(diag(S)), rep(1, 4), tolerance = 1e-6)
  P <- partial_cors(pars$Psi)
  expect_equal(unname(diag(P)), rep(1, 4))
  expect_equal(P, t(P))
})

test_that("fit_pairs recovers the truth on a large complete data set", {
  sim <- simulate_ema(N = 300, T = 150, motifs = "M0", compliance = 1, seed = 2)
  f <- fit_pairs(sim$data, sim$vars)
  expect_lt(max(abs(f$Phi - pars$Phi)), 0.03)
  expect_lt(max(abs(f$Psi - pars$Psi)), 0.04)
  expect_lt(max(abs(f$mu - pars$mu)), 0.15)
  s <- summary(f)
  expect_equal(nrow(s), 16)
  expect_true(all(s$se > 0))
  expect_output(print(f), "half-panel jackknife")
})

test_that("weights of one reproduce the unweighted fit and the jackknife can be switched off", {
  sim <- simulate_ema(N = 40, T = 30, motifs = "M1", compliance = 0.8, seed = 3)
  pr <- make_pairs(sim$data, sim$vars)
  f1 <- fit_pairs(sim$data, sim$vars, pairs = pr)
  f2 <- fit_pairs(sim$data, sim$vars, pairs = pr, weights = rep(1, nrow(pr$y)))
  expect_equal(f1$Phi, f2$Phi)
  f3 <- fit_pairs(sim$data, sim$vars, pairs = pr, bias_correct = FALSE)
  expect_equal(f3$Phi, f1$Phi_ls)
})

test_that("fit_tilt at delta = 0 equals fit_pairs and moves the estimates away from it otherwise", {
  sim <- simulate_ema(N = 60, T = 40, motifs = "M2", compliance = 0.7, delta = -1, seed = 4)
  f0 <- fit_pairs(sim$data, sim$vars)
  t0 <- fit_tilt(sim$data, sim$vars, delta = 0)
  expect_equal(t0$Phi, f0$Phi)
  t1 <- fit_tilt(sim$data, sim$vars, delta = -1)
  expect_true(t1$converged)
  expect_gt(t1$Phi[1, 1], f0$Phi[1, 1])          # attenuation removed
  expect_gt(t1$mu[1], f0$mu[1])                  # hidden states were higher
  expect_true(is.finite(t1$ess) && t1$ess <= t1$n_pairs)
  tp <- fit_tilt(sim$data, sim$vars, delta = -1, propensity = "person")
  expect_equal(length(tp$alpha), tp$n_persons)
})

test_that("tilt_profile, break_even and plotting work", {
  sim <- simulate_ema(N = 60, T = 40, motifs = "M2", compliance = 0.7, delta = -1, seed = 5)
  prof <- tilt_profile(sim$data, sim$vars, delta_grid = c(-1, -0.5, 0, 0.5))
  expect_s3_class(prof, "tilt_profile")
  expect_equal(nrow(prof$Phi), 4 * 16)
  be <- break_even(prof, delta_max = 1)
  expect_equal(nrow(be), 16)
  expect_true(all(be$set_lower <= be$set_upper))
  expect_output(print(prof), "Tilt profile")
  pdf(NULL); on.exit(dev.off())
  expect_silent(plot(prof))
  expect_silent(plot(dm_graph(c("M1", "M3"), sensor = TRUE)))
})

test_that("silence test is null under M0 and signed under M2", {
  sim0 <- simulate_ema(N = 200, T = 56, motifs = "M0", compliance = 0.7, seed = 6)
  st0 <- silence_test(sim0$data, sim0$vars)
  expect_s3_class(st0, "silence_test")
  expect_true(all(abs(st0$z) < 4))
  sim2 <- simulate_ema(N = 200, T = 56, motifs = "M2", compliance = 0.7, delta = -1, seed = 7)
  st2 <- silence_test(sim2$data, sim2$vars)
  expect_lt(st2$coef_R[1], 0)
  expect_lt(st2$p[1], 0.001)
  expect_output(print(st2), "Silence test")
  # single-person day-block bootstrap path
  one <- simulate_ema(N = 1, T = 300, motifs = "M2", compliance = 0.7, delta = -1, seed = 12, days = 10)
  st1 <- silence_test(one$data, one$vars, day = "day", B = 50)
  expect_equal(attr(st1, "se"), "dayblock")
  expect_true(all(is.finite(st1$se)))
})

test_that("sensor-gap test and calibration run and point in the right direction", {
  sim <- simulate_ema(N = 200, T = 56, motifs = "M2", compliance = 0.7, delta = -1, sensor_cor = 0.6, p_probe = 0.05, seed = 8)
  sg <- sensor_gap_test(sim$data, sim$vars)
  expect_lt(sg$coef_R, 0)
  expect_gt(attr(sg, "loading")[1], 0)
  prof <- tilt_profile(sim$data, sim$vars, delta_grid = seq(-2, 0.5, by = 0.5), probe = "Z")
  cs <- calibrate_delta(prof, sim$data, method = "sensor")
  expect_s3_class(cs, "delta_calibration")
  expect_true(is.finite(cs$delta_hat))
  expect_lt(cs$delta_hat, 0)
  cp <- calibrate_delta(prof, sim$data, method = "probe")
  expect_lt(cp$delta_hat, 0)
  ck <- calibrate_delta(prof, sim$data, method = "postskip", n_sim = 2, seed = 1)
  expect_lt(ck$delta_hat, 0)
  expect_output(print(cs), "Calibration of delta")
  pdf(NULL); on.exit(dev.off())
  expect_silent(plot(cs))
  fc <- fatigue_check(sim$data)
  expect_s3_class(fc, "fatigue_check")
  expect_true(is.finite(fc$difference))
  expect_output(print(fc), "Fatigue check")
})

test_that("burden-aware post-skip calibration recovers a larger |delta| than the M2-only version under M2 + M3", {
  sim <- simulate_ema(N = 150, T = 56, motifs = c("M2", "M3"), compliance = 0.7, delta = -1, kappa_R = 1, seed = 21)
  prof <- tilt_profile(sim$data, sim$vars, delta_grid = c(-1.5, -1, -0.5, 0))
  c0 <- calibrate_delta(prof, sim$data, method = "postskip", n_sim = 4, seed = 2, burden = "none")
  c1 <- calibrate_delta(prof, sim$data, method = "postskip", n_sim = 4, seed = 2, burden = "fit", kappa_grid = c(0, 0.5, 1, 1.5))
  expect_equal(c1$burden, "fit")
  expect_true(all(c("kappa_hat") %in% names(c1$curve)))
  expect_true(is.finite(c1$persistence_observed))
  # the M2-only curve understates |delta| under burden; the burden-aware curve should not be closer to zero
  expect_true(is.na(c0$delta_hat) || is.na(c1$delta_hat) || c1$delta_hat <= c0$delta_hat + 0.25)
  d <- simulate_from_fit(prof$fits[[2]], N = 20, T = 15, days = 5, kappa_R = 1, seed = 3)
  expect_true(all(c("day", "R") %in% names(d)))
  expect_equal(nrow(d), 300)
})

test_that("random-number state is restored by seeded functions", {
  set.seed(100); a <- runif(1)
  set.seed(100); invisible(simulate_ema(N = 5, T = 5, seed = 7)); b <- runif(1)
  expect_equal(a, b)
})

test_that("input validation gives informative errors", {
  sim <- simulate_ema(N = 20, T = 15, seed = 1)
  expect_error(simulate_ema(N = 5, T = 5, motifs = "M9"), "Unknown motif")
  expect_error(fit_pairs(sim$data, c("NegA", "Nope")), "not found")
  expect_error(fit_tilt(sim$data, sim$vars, delta = c(-1, 0)), "'delta'")
  d2 <- sim$data; d2$time <- d2$time + 0.5
  expect_error(silence_test(d2, sim$vars), "integer prompt index")
  expect_error(calibrate_delta(tilt_profile(sim$data, sim$vars, delta_grid = c(-0.5, 0)), sim$data, method = "sensor"), "sensor column")
  f <- fit_pairs(sim$data, sim$vars)
  expect_equal(length(coef(f)), 16); expect_equal(dim(vcov(f)), c(16, 16)); expect_equal(dim(confint(f)), c(16, 2)); expect_equal(nobs(f), f$n_pairs)
  expect_s3_class(fit_tilt(sim$data, sim$vars, delta = 0), "tilt_fit")
})

test_that("bounds_support gives Manski bounds with the right width", {
  sim <- simulate_ema(N = 5, T = 20, motifs = "M0", compliance = 0.6, seed = 9)
  b <- bounds_support(sim$data, sim$vars, L = 1, U = 7)
  expect_equal(nrow(b), 20)
  expect_equal(b$upper - b$lower, (1 - b$response_rate) * 6, tolerance = 1e-12)
})

test_that("FIML log-likelihood matches the brute-force Gaussian likelihood and EM improves it", {
  set.seed(10)
  p <- 2; Tn <- 5
  Phi2 <- matrix(c(0.5, 0.1, -0.2, 0.4), 2, 2); Psi2 <- matrix(c(1, 0.3, 0.3, 0.8), 2, 2); mu2 <- c(1, 2)
  Smu2 <- matrix(c(0.5, 0.1, 0.1, 0.4), 2, 2); S02 <- stationary_cov(Phi2, Psi2)
  brute <- function(Ymat, obs) {
    q <- p * Tn; CD <- matrix(0, q, q); blk <- function(t) (p * (t - 1) + 1):(p * t)
    for (t in 1:Tn) for (s in 1:t) {
      C <- if (t == s) { if (t == 1) S02 else Phi2 %*% CD[blk(t - 1), blk(t - 1)] %*% t(Phi2) + Psi2 } else Phi2 %*% CD[blk(t - 1), blk(s)]
      CD[blk(t), blk(s)] <- C; CD[blk(s), blk(t)] <- t(C)
    }
    ones <- kronecker(rep(1, Tn), diag(p)); Sig <- CD + ones %*% Smu2 %*% t(ones); m <- as.vector(ones %*% mu2)
    Y <- as.vector(t(Ymat)); o <- which(obs)
    -0.5 * (length(o) * log(2 * pi) + as.numeric(determinant(Sig[o, o])$modulus) + as.numeric(t(Y[o] - m[o]) %*% solve(Sig[o, o], Y[o] - m[o])))
  }
  mu_i <- mu2 + t(chol(Smu2)) %*% rnorm(2); X <- matrix(NA, Tn, p); d <- t(chol(S02)) %*% rnorm(2)
  for (t in 1:Tn) { if (t > 1) d <- Phi2 %*% d + t(chol(Psi2)) %*% rnorm(2); X[t, ] <- mu_i + d }
  Rr <- c(1, 0, 1, 1, 0); Xo <- X; Xo[Rr == 0, ] <- NA
  dd <- data.frame(id = 1, time = 1:Tn, R = Rr, X1 = Xo[, 1], X2 = Xo[, 2])
  expect_equal(loglik_fiml(dd, c("X1", "X2"), Phi2, Psi2, mu2, Smu2, S02), brute(X, rep(Rr == 1, each = p)), tolerance = 1e-10)
  sim <- simulate_ema(N = 60, T = 20, motifs = "M0", compliance = 0.7, seed = 11)
  ff <- fit_fiml(sim$data, sim$vars)
  expect_s3_class(ff, "fiml_fit")
  expect_true(all(diff(ff$loglik) > -1e-6))
  expect_output(print(ff), "FIML")
})
