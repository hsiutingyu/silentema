// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// Full-information maximum likelihood for the two-level VAR(1) with random
// person means under a missing-at-random assumption, by the EM algorithm of
// Shumway and Stoffer (1982) applied to the augmented state (X_t, mu_i).
//
// State:       a_t = (X_t, mu_i),  a_t = F a_{t-1} + e_t,  e_t ~ N(0, Q)
//              F = [Phi, I - Phi; 0, I],  Q = blockdiag(Psi, 0)
// Observation: y_t = H a_t exactly (H = [I, 0]) when R_t = 1, absent otherwise
// Initial:     a_1 ~ N((mu, mu), [Sigma0 + Sigma_mu, Sigma_mu; Sigma_mu, Sigma_mu])
//
// Y: cube p x T x N (NaN where missing), R: T x N integer matrix.

static mat sym(const mat& A) { return 0.5 * (A + A.t()); }

// [[Rcpp::export]]
List kalman_em_cpp(const arma::cube& Y, const arma::imat& R, arma::mat Phi, arma::mat Psi,
                   arma::vec mu, arma::mat Sigma_mu, arma::mat Sigma0,
                   int max_iter = 200, double tol = 1e-6, bool verbose = false) {
  const int p = Y.n_rows, T = Y.n_cols, N = Y.n_slices, q = 2 * p;
  mat I_p = eye(p, p);
  mat H = join_rows(I_p, zeros(p, p));         // p x q
  mat L = join_rows(I_p, -I_p);                // D_t = X_t - mu_i
  std::vector<double> ll_trace;
  bool converged = false;
  int iter = 0;
  double ll_old = -datum::inf;
  for (iter = 1; iter <= max_iter; ++iter) {
    Rcpp::checkUserInterrupt();
    mat F = join_cols(join_rows(Phi, I_p - Phi), join_rows(zeros(p, p), I_p));
    mat Q = zeros(q, q); Q.submat(0, 0, p - 1, p - 1) = Psi;
    vec m1 = join_cols(mu, mu);
    mat P1 = join_cols(join_rows(Sigma0 + Sigma_mu, Sigma_mu), join_rows(Sigma_mu, Sigma_mu));
    // accumulators
    mat S11 = zeros(p, p), S10 = zeros(p, p), S00 = zeros(p, p);
    vec mu_acc = zeros(p); mat Smu_acc = zeros(p, p); mat S0_acc = zeros(p, p);
    double ll = 0.0; int n_trans = 0;
    // storage
    std::vector<vec> mp(T), mf(T), ms(T);
    std::vector<mat> Pp(T), Pf(T), Ps(T), J(T);
    for (int i = 0; i < N; ++i) {
      // ---- filter
      for (int t = 0; t < T; ++t) {
        if (t == 0) { mp[t] = m1; Pp[t] = P1; }
        else { mp[t] = F * mf[t - 1]; Pp[t] = sym(F * Pf[t - 1] * F.t() + Q); }
        if (R(t, i) == 1) {
          vec y = Y.slice(i).col(t);
          vec v = y - H * mp[t];
          mat S = sym(Pp[t].submat(0, 0, p - 1, p - 1));
          mat Sinv; bool ok = inv_sympd(Sinv, S);
          if (!ok) Sinv = pinv(S);
          mat K = Pp[t] * H.t() * Sinv;                 // q x p
          mf[t] = mp[t] + K * v;
          mat IKH = eye(q, q) - K * H;
          Pf[t] = sym(IKH * Pp[t] * IKH.t());           // Joseph form (no measurement noise)
          double ldet, sgn; log_det(ldet, sgn, S);
          ll += -0.5 * (p * std::log(2.0 * datum::pi) + ldet + as_scalar(v.t() * Sinv * v));
        } else { mf[t] = mp[t]; Pf[t] = Pp[t]; }
      }
      // ---- smoother
      ms[T - 1] = mf[T - 1]; Ps[T - 1] = Pf[T - 1];
      for (int t = T - 2; t >= 0; --t) {
        mat Ppinv; bool ok = inv_sympd(Ppinv, Pp[t + 1]);
        if (!ok) Ppinv = pinv(Pp[t + 1]);
        J[t] = Pf[t] * F.t() * Ppinv;
        ms[t] = mf[t] + J[t] * (ms[t + 1] - mp[t + 1]);
        Ps[t] = sym(Pf[t] + J[t] * (Ps[t + 1] - Pp[t + 1]) * J[t].t());
      }
      // ---- accumulate sufficient statistics
      for (int t = 1; t < T; ++t) {
        mat Plag = Ps[t] * J[t - 1].t();                // Cov(a_t, a_{t-1} | Y)
        mat E_tt   = L * (Ps[t] + ms[t] * ms[t].t()) * L.t();
        mat E_t1t1 = L * (Ps[t - 1] + ms[t - 1] * ms[t - 1].t()) * L.t();
        mat E_tt1  = L * (Plag + ms[t] * ms[t - 1].t()) * L.t();
        S00 += E_tt; S11 += E_t1t1; S10 += E_tt1; ++n_trans;
      }
      vec mu_i = ms[0].subvec(p, q - 1);
      mat Pmu = Ps[0].submat(p, p, q - 1, q - 1);
      mu_acc += mu_i; Smu_acc += Pmu + mu_i * mu_i.t();
      vec d0 = L * ms[0]; mat Pd0 = L * Ps[0] * L.t();
      S0_acc += Pd0 + d0 * d0.t();
    }
    ll_trace.push_back(ll);
    // ---- M-step
    mat Phi_new = S10 * inv_sympd(sym(S11));
    mat Psi_new = sym((S00 - Phi_new * S10.t()) / n_trans);
    vec mu_new = mu_acc / N;
    mat Smu_new = sym(Smu_acc / N - mu_new * mu_new.t());
    mat S0_new = sym(S0_acc / N);
    // guard positive definiteness
    vec ev; mat evec;
    // keep the covariance matrices positive definite (a zero variance component would be an absorbing state)
    eig_sym(ev, evec, Smu_new); if (ev.min() < 1e-6) { ev = clamp(ev, 1e-6, datum::inf); Smu_new = evec * diagmat(ev) * evec.t(); }
    eig_sym(ev, evec, Psi_new); if (ev.min() < 1e-6) { ev = clamp(ev, 1e-6, datum::inf); Psi_new = evec * diagmat(ev) * evec.t(); }
    eig_sym(ev, evec, S0_new); if (ev.min() < 1e-6) { ev = clamp(ev, 1e-6, datum::inf); S0_new = evec * diagmat(ev) * evec.t(); }
    double change = std::max({abs(Phi_new - Phi).max(), abs(Psi_new - Psi).max(), abs(mu_new - mu).max(), abs(Smu_new - Sigma_mu).max()});
    Phi = Phi_new; Psi = Psi_new; mu = mu_new; Sigma_mu = Smu_new; Sigma0 = S0_new;
    if (verbose) Rcpp::Rcout << "iter " << iter << " loglik " << ll << " change " << change << std::endl;
    if (change < tol || (iter > 5 && std::abs(ll - ll_old) < 1e-7 * std::max(1.0, std::abs(ll)))) { converged = true; break; }
    ll_old = ll;
  }
  return List::create(_["Phi"] = Phi, _["Psi"] = Psi, _["mu"] = mu, _["Sigma_mu"] = Sigma_mu,
                      _["Sigma0"] = Sigma0, _["loglik"] = ll_trace, _["iterations"] = std::min(iter, max_iter),
                      _["converged"] = converged);
}

// Log-likelihood only (for checks): same filter, no EM.
// [[Rcpp::export]]
double kalman_loglik_cpp(const arma::cube& Y, const arma::imat& R, const arma::mat& Phi, const arma::mat& Psi,
                         const arma::vec& mu, const arma::mat& Sigma_mu, const arma::mat& Sigma0) {
  const int p = Y.n_rows, T = Y.n_cols, N = Y.n_slices, q = 2 * p;
  mat I_p = eye(p, p);
  mat H = join_rows(I_p, zeros(p, p));
  mat F = join_cols(join_rows(Phi, I_p - Phi), join_rows(zeros(p, p), I_p));
  mat Q = zeros(q, q); Q.submat(0, 0, p - 1, p - 1) = Psi;
  vec m1 = join_cols(mu, mu);
  mat P1 = join_cols(join_rows(Sigma0 + Sigma_mu, Sigma_mu), join_rows(Sigma_mu, Sigma_mu));
  double ll = 0.0;
  for (int i = 0; i < N; ++i) {
    vec m = m1; mat P = P1;
    for (int t = 0; t < T; ++t) {
      if (t > 0) { m = F * m; P = sym(F * P * F.t() + Q); }
      if (R(t, i) == 1) {
        vec y = Y.slice(i).col(t);
        vec v = y - H * m;
        mat S = sym(P.submat(0, 0, p - 1, p - 1));
        mat Sinv; bool ok = inv_sympd(Sinv, S); if (!ok) Sinv = pinv(S);
        mat K = P * H.t() * Sinv;
        m = m + K * v;
        mat IKH = eye(q, q) - K * H;
        P = sym(IKH * P * IKH.t());
        double ldet, sgn; log_det(ldet, sgn, S);
        ll += -0.5 * (p * std::log(2.0 * datum::pi) + ldet + as_scalar(v.t() * Sinv * v));
      }
    }
  }
  return ll;
}
