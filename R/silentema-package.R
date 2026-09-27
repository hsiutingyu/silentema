#' silentema: dynamic missingness graphs and sensitivity analysis for EMA data
#'
#' Participants in ecological momentary assessment (EMA) and other
#' experience-sampling studies skip prompts, and the reasons are rarely
#' unrelated to the momentary states being measured. \pkg{silentema}
#' provides a graphical language for stating which nonresponse mechanism is
#' assumed, a checker that says which quantities of a two-level vector
#' autoregressive model of order one (VAR(1)) remain estimable under that
#' mechanism and by which estimator, two tests that ask the data whether
#' skipped prompts were informative, estimators for the recoverable case, a
#' sensitivity analysis for self-censoring (the mechanism that estimation
#' alone cannot repair), three calibration designs for the sensitivity
#' parameter, worst-case bounds, a reporting template, and a simulator.
#'
#' @section The workflow:
#' The recommended order of analysis is
#' \enumerate{
#'   \item \strong{Declare} the assumed mechanism as a dynamic missingness
#'     graph with \code{\link{dm_graph}}, built from the motifs M0 (completely
#'     random), M1 (lagged-state dependence), M2 (self-censoring), M3 (burden
#'     or fatigue), M4 (person propensity), M5 (context confounding) and M6
#'     (reactivity), with optional sensor, probe and context nodes; inspect it
#'     with \code{plot} and query it with \code{\link{dsep}}.
#'   \item \strong{Check} recoverability with \code{\link{recoverability}}:
#'     which estimands (transition kernel, person means, between-person law)
#'     are structurally recoverable, and which estimator recovers them.
#'   \item \strong{Test} whether silence was informative with
#'     \code{\link{silence_test}} (the state after a skipped prompt) and, when
#'     an always-observed channel exists, \code{\link{sensor_gap_test}};
#'     describe response persistence with \code{\link{fatigue_check}}.
#'   \item \strong{Estimate} the within-person VAR(1) from answered adjacent
#'     prompts with \code{\link{fit_pairs}} (person intercepts, half-panel
#'     jackknife, cluster-robust inference through \code{\link{summary.pairs_fit}}),
#'     with covariate adjustment or inverse-probability weighting for an
#'     observed context (\code{\link{fit_ipw}}), or by full-information
#'     maximum likelihood under missing at random (\code{\link{fit_fiml}}).
#'   \item \strong{Profile} the estimates over a self-censoring sensitivity
#'     value with \code{\link{fit_tilt}} and \code{\link{tilt_profile}}, and
#'     summarize the profile with \code{\link{break_even}} (sign changes,
#'     significance changes, identified sets and bands).
#'   \item \strong{Calibrate} the sensitivity value from a passive sensor,
#'     randomized probes or the post-skip contrast with
#'     \code{\link{calibrate_delta}}; compare with the worst-case bounds of
#'     \code{\link{bounds_support}}.
#'   \item \strong{Report} with \code{\link{missingness_declaration}}, which
#'     writes a Markdown declaration for preregistrations and papers.
#' }
#'
#' @section Simulation and design:
#' \code{\link{simulate_ema}} generates two-level VAR(1) data under any
#' combination of motifs with a target response rate, optional sensor, probe,
#' context and day structure; \code{\link{simulate_from_fit}} simulates from a
#' fitted (tilted) model; \code{\link{default_params}},
#' \code{\link{stationary_cov}} and \code{\link{partial_cors}} give the
#' population parameters used in the accompanying article and the derived
#' quantities.
#'
#' @section Data format:
#' All estimators take a long data frame with one row per \emph{scheduled}
#' prompt, answered or not: a person identifier (\code{id}), an integer prompt
#' index that increases by one from one scheduled prompt to the next
#' (\code{time}), a response indicator (\code{R}, 1 = answered), the state
#' columns (\code{NA} at skipped prompts), and optionally a \code{day} column
#' (pairs are then formed within days only), a sensor column, a probe column
#' and context columns. See \code{\link{make_pairs}}.
#'
#' @references
#' Yu, H.-T. (2026). What skipped prompts hide: Detecting, diagnosing, and
#' correcting informative nonresponse in ecological momentary assessment.
#' Manuscript under review. Materials: \url{https://osf.io/x6d2t/}.
#'
#' Mohan, K., & Pearl, J. (2021). Graphical models for processing missing
#' data. \emph{Journal of the American Statistical Association, 116},
#' 1023-1037. \doi{10.1080/01621459.2021.1874961}
#'
#' Dhaene, G., & Jochmans, K. (2015). Split-panel jackknife estimation of
#' fixed-effect models. \emph{The Review of Economic Studies, 82}, 991-1030.
#' \doi{10.1093/restud/rdv007}
#'
#' Shumway, R. H., & Stoffer, D. S. (1982). An approach to time series
#' smoothing and forecasting using the EM algorithm. \emph{Journal of Time
#' Series Analysis, 3}, 253-264. \doi{10.1111/j.1467-9892.1982.tb00349.x}
#'
#' @seealso The vignettes: \code{vignette("silentema-workflow")} for a
#'   start-to-finish analysis, \code{vignette("dm-graphs")} for the motif
#'   taxonomy and recoverability, \code{vignette("testing-informativeness")}
#'   for the tests, \code{vignette("sensitivity-analysis")} for tilting and
#'   calibration, and \code{vignette("simulation-and-design")} for the
#'   simulator.
#' @aliases silentema
#' @useDynLib silentema, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom stats pnorm dnorm qnorm rnorm runif rbinom uniroot cov complete.cases setNames coef vcov confint nobs
"_PACKAGE"
