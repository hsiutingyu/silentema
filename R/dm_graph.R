#' Declare a dynamic missingness graph (dm-graph)
#'
#' A dm-graph is a directed acyclic graph over the unrolled within-person
#' process of an experience-sampling study: person-level latent nodes
#' (\code{eta}, the person's typical state; \code{zeta}, a person-level
#' response propensity), momentary states \code{X_t}, response indicators
#' \code{R_t}, and optional context (\code{C_t}), sensor (\code{S_t}) and
#' probe (\code{Z_t}) nodes. The structural edges \code{X_{t-1} -> X_t} and
#' \code{eta -> X_t} are always present. The missingness mechanism is
#' declared as a set of motifs, each of which adds edges into \code{R_t}.
#'
#' @param motifs Character vector of motif codes. \code{"M0"} (no edges into
#'   R), \code{"M1"} (\code{X_{t-1} -> R_t}, lagged-state dependence),
#'   \code{"M2"} (\code{X_t -> R_t}, self-censoring), \code{"M3"}
#'   (\code{R_{t-1} -> R_t}, burden or fatigue), \code{"M4"}
#'   (\code{zeta -> R_t} with \code{zeta} associated with \code{eta},
#'   person propensity), \code{"M5"} (\code{C_t -> X_t} and \code{C_t -> R_t},
#'   context confounding), \code{"M6"} (\code{R_{t-1} -> X_t}, reactivity).
#'   Motifs may be combined.
#' @param context One of \code{"none"}, \code{"latent"}, \code{"observed"}.
#'   Relevant for \code{"M5"}; \code{"observed"} means the context variable is
#'   recorded at every prompt (for example by a phone sensor).
#' @param context_persistent Logical; if \code{TRUE} the context process has
#'   edges \code{C_{t-1} -> C_t}.
#' @param sensor Logical; add an always-observed passive-sensor node
#'   \code{S_t} with \code{X_t -> S_t}.
#' @param probe Logical; add a randomized probe indicator \code{Z_t} with
#'   \code{Z_t -> R_t} (a probe forces a response).
#' @param window Integer number of prompts in the unrolled graph (at least 4).
#'
#' @return An object of class \code{"dm_graph"}: a list with elements
#'   \code{nodes} (character vector of node names), \code{edges} (two-column
#'   character matrix \code{from}, \code{to}), \code{motifs}, \code{context},
#'   \code{context_persistent}, \code{sensor}, \code{probe}, \code{window},
#'   \code{latent} (names of the latent nodes) and \code{observed_always}
#'   (names of the nodes observed at every prompt). It has \code{print} and
#'   \code{\link[=plot.dm_graph]{plot}} methods.
#' @examples
#' g <- dm_graph(c("M1", "M4"))
#' g
#' recoverability(g)
#' g2 <- dm_graph(c("M2", "M3"), sensor = TRUE)
#' recoverability(g2)
#' # an observed context adds M5 and observed context nodes
#' g5 <- dm_graph("M1", context = "observed", context_persistent = TRUE)
#' g5$edges[g5$edges[, "from"] == "C2", ]
#' @seealso \code{\link{recoverability}}, \code{\link{dsep}},
#'   \code{\link{plot.dm_graph}}, \code{\link{simulate_ema}} (which generates
#'   data under the same motifs)
#' @export
dm_graph <- function(motifs = "M0", context = c("none", "latent", "observed"),
                     context_persistent = FALSE, sensor = FALSE, probe = FALSE,
                     window = 5L) {
  context <- match.arg(context)
  motifs <- toupper(unique(motifs))
  valid <- c("M0", "M1", "M2", "M3", "M4", "M5", "M6")
  bad <- setdiff(motifs, valid)
  if (length(bad)) stop("Unknown motif code(s): ", paste(bad, collapse = ", "))
  if ("M5" %in% motifs && context == "none") { context <- "latent"; message("M5 declared without 'context': a latent context is assumed") }
  if (context != "none" && !("M5" %in% motifs)) { motifs <- c(motifs, "M5"); message("'context' given: motif M5 added to the declaration") }
  if (length(motifs) > 1L) motifs <- setdiff(motifs, "M0")
  if (window < 4L) stop("window must be at least 4")
  W <- as.integer(window)
  Xn <- paste0("X", seq_len(W)); Rn <- paste0("R", seq_len(W))
  edges <- list()
  add <- function(from, to) edges[[length(edges) + 1L]] <<- c(from, to)
  # structural edges
  for (t in seq_len(W)) add("eta", Xn[t])
  for (t in 2:W) add(Xn[t - 1], Xn[t])
  nodes <- c("eta", Xn, Rn)
  latent <- "eta"
  observed_always <- Rn
  # M4: zeta associated with eta through a common latent U
  if ("M4" %in% motifs) {
    nodes <- c(nodes, "U", "zeta"); latent <- c(latent, "U", "zeta")
    add("U", "eta"); add("U", "zeta")
    for (t in seq_len(W)) add("zeta", Rn[t])
  }
  if ("M1" %in% motifs) for (t in 2:W) add(Xn[t - 1], Rn[t])
  if ("M2" %in% motifs) for (t in seq_len(W)) add(Xn[t], Rn[t])
  if ("M3" %in% motifs) for (t in 2:W) add(Rn[t - 1], Rn[t])
  if ("M6" %in% motifs) for (t in 2:W) add(Rn[t - 1], Xn[t])
  if (context != "none") {
    Cn <- paste0("C", seq_len(W)); nodes <- c(nodes, Cn)
    if (context == "latent") latent <- c(latent, Cn) else observed_always <- c(observed_always, Cn)
    for (t in seq_len(W)) { add(Cn[t], Xn[t]); add(Cn[t], Rn[t]) }
    if (context_persistent) for (t in 2:W) add(Cn[t - 1], Cn[t])
  }
  if (sensor) {
    Sn <- paste0("S", seq_len(W)); nodes <- c(nodes, Sn); observed_always <- c(observed_always, Sn)
    for (t in seq_len(W)) add(Xn[t], Sn[t])
  }
  if (probe) {
    Zn <- paste0("Z", seq_len(W)); nodes <- c(nodes, Zn); observed_always <- c(observed_always, Zn)
    for (t in seq_len(W)) add(Zn[t], Rn[t])
  }
  em <- do.call(rbind, edges)
  colnames(em) <- c("from", "to")
  structure(list(nodes = nodes, edges = em, motifs = motifs, context = context,
                 context_persistent = context_persistent, sensor = sensor,
                 probe = probe, window = W, latent = latent,
                 observed_always = observed_always),
            class = "dm_graph")
}

#' @export
print.dm_graph <- function(x, ...) {
  cat("Dynamic missingness graph (window of", x$window, "prompts)\n")
  cat("  motifs:", paste(x$motifs, collapse = " + "), "\n")
  cat("  context:", x$context, if (x$context_persistent) "(persistent)" else "", "\n")
  cat("  sensor:", x$sensor, "  probe:", x$probe, "\n")
  cat("  nodes:", length(x$nodes), "  edges:", nrow(x$edges), "\n")
  invisible(x)
}

#' d-separation in a dm-graph
#'
#' Tests whether two node sets are d-separated given a conditioning set,
#' using the reachability ("Bayes-ball") algorithm (Koller & Friedman, 2009,
#' Algorithm 3.1).
#'
#' @param g A \code{dm_graph} (or any list with an \code{edges} matrix with
#'   columns \code{from}, \code{to}).
#' @param x,y Character vectors of node names (window indices: \code{"X3"} is
#'   the state at the third prompt of the window; \code{recoverability} uses
#'   the third prompt as \eqn{t}).
#' @param z Character vector of conditioning nodes (may be empty).
#' @return A single logical value: \code{TRUE} if every node in \code{x} is
#'   d-separated from every node in \code{y} given \code{z}, \code{FALSE}
#'   otherwise.
#' @references Koller, D., & Friedman, N. (2009). \emph{Probabilistic graphical
#'   models: Principles and techniques}. MIT Press.
#'
#'   Shachter, R. D. (1998). Bayes-ball: The rational pastime (for determining
#'   irrelevance and requisite information in belief networks and influence
#'   diagrams). In \emph{Proceedings of the Fourteenth Conference on
#'   Uncertainty in Artificial Intelligence} (pp. 480-487). Morgan Kaufmann.
#' @examples
#' g <- dm_graph("M1")
#' dsep(g, "R3", "X3", c("X2", "eta"))      # TRUE: the kernel is recoverable
#' g2 <- dm_graph("M2")
#' dsep(g2, "R3", "X3", c("X2", "eta"))     # FALSE: self-censoring
#' # any edge list works, not only dm-graphs: a collider A -> B <- C
#' coll <- list(nodes = c("A", "B", "C"),
#'              edges = cbind(from = c("A", "C"), to = c("B", "B")))
#' dsep(coll, "A", "C")          # TRUE
#' dsep(coll, "A", "C", "B")     # FALSE: conditioning on the collider opens the path
#' @export
dsep <- function(g, x, y, z = character(0)) {
  em <- g$edges
  nodes <- unique(c(em[, 1], em[, 2], g$nodes))
  miss <- setdiff(c(x, y, z), nodes)
  if (length(miss)) stop("Unknown node(s): ", paste(miss, collapse = ", "))
  parents <- split(em[, "from"], factor(em[, "to"], levels = nodes))
  children <- split(em[, "to"], factor(em[, "from"], levels = nodes))
  # ancestors of z (including z)
  anc <- character(0); stack <- z
  while (length(stack)) {
    v <- stack[1]; stack <- stack[-1]
    if (!(v %in% anc)) { anc <- c(anc, v); stack <- c(stack, parents[[v]]) }
  }
  # traversal: (node, direction) direction "up" = arrived from a child, "down" = arrived from a parent
  visited <- character(0)
  queue <- lapply(x, function(v) list(node = v, dir = "up"))
  reached <- character(0)
  while (length(queue)) {
    cur <- queue[[1]]; queue <- queue[-1]
    key <- paste(cur$node, cur$dir)
    if (key %in% visited) next
    visited <- c(visited, key)
    v <- cur$node
    if (!(v %in% z)) reached <- c(reached, v)
    if (cur$dir == "up" && !(v %in% z)) {
      for (p in parents[[v]]) queue[[length(queue) + 1L]] <- list(node = p, dir = "up")
      for (ch in children[[v]]) queue[[length(queue) + 1L]] <- list(node = ch, dir = "down")
    } else if (cur$dir == "down") {
      if (!(v %in% z)) for (ch in children[[v]]) queue[[length(queue) + 1L]] <- list(node = ch, dir = "down")
      if (v %in% anc) for (p in parents[[v]]) queue[[length(queue) + 1L]] <- list(node = p, dir = "up")
    }
  }
  !any(y %in% reached)
}

#' Recoverability report for a dm-graph
#'
#' Applies the d-separation conditions derived in Yu (2026) to a declared
#' dm-graph and reports, estimand by estimand, whether the estimand is
#' structurally recoverable and by which estimator, and whether the silence
#' test and the sensor-gap test are valid tests of the recoverable class under
#' the declared mechanism. Recoverability is used in the sense of Mohan and
#' Pearl (2021): a consistent estimator exists that uses only the observed
#' part of the data.
#'
#' @param g A \code{dm_graph}.
#' @param t Index of the focal prompt inside the window (default: the third).
#' @return A data frame of class \code{"recoverability"} with one row per
#'   estimand and the columns \code{estimand}, \code{recoverable} (logical;
#'   \code{NA} for the two test rows, which report expected behavior rather
#'   than recoverability), \code{estimator} (the estimator that recovers the
#'   estimand, or the reason it is not recoverable), \code{condition} (the
#'   d-separation statement that was checked), and the attributes
#'   \code{motifs} and \code{silence_null} (logical: whether the null
#'   hypothesis of the silence test is implied by the graph).
#' @examples
#' recoverability(dm_graph("M1"))
#' recoverability(dm_graph(c("M1", "M3")))
#' recoverability(dm_graph("M6"))
#' r <- recoverability(dm_graph("M2", sensor = TRUE))
#' r$recoverable
#' attr(r, "silence_null")
#' @details The within-person conditioning set always contains every
#'   person-level latent node (\code{eta}, and \code{zeta} when declared),
#'   which is what person-specific intercepts (fixed effects) realize in
#'   estimation. Observed context nodes are added to the conditioning set for
#'   the context-conditional transition kernel. The probe rule
#'   (\code{Z_t = 1} forces \code{R_t = 1}) is a deterministic relation that
#'   d-separation cannot express; it is applied as a separate rule.
#' @references Mohan, K., & Pearl, J. (2021). Graphical models for processing
#'   missing data. \emph{Journal of the American Statistical Association, 116},
#'   1023-1037. \doi{10.1080/01621459.2021.1874961}
#'
#'   Yu, H.-T. (2026). What skipped prompts hide: Detecting, diagnosing, and
#'   correcting informative nonresponse in ecological momentary assessment.
#'   Manuscript under review.
#' @seealso \code{\link{dm_graph}}, \code{\link{dsep}},
#'   \code{\link{missingness_declaration}}
#' @export
recoverability <- function(g, t = 3L) {
  W <- g$window
  if (t < 3L || t > W - 1L) stop("t must satisfy 3 <= t <= window - 1")
  X <- function(k) paste0("X", k); R <- function(k) paste0("R", k)
  pl <- intersect(c("eta", "zeta"), g$nodes)         # person-level latents
  Cobs <- if (g$context == "observed") paste0("C", t) else character(0)
  rows <- list()
  push <- function(estimand, ok, estimator, cond) rows[[length(rows) + 1L]] <<-
    data.frame(estimand = estimand, recoverable = ok, estimator = estimator, condition = cond,
               stringsAsFactors = FALSE)
  # --- transition kernel (Phi, Psi) from complete adjacent pairs
  za <- c(X(t - 1), pl, Cobs)
  a <- dsep(g, R(t), X(t), za)
  b <- dsep(g, R(t - 1), X(t), c(za, R(t)))
  reactive <- any(g$edges[, "from"] == R(t - 1) & g$edges[, "to"] == X(t))
  kern_ok <- a && b
  if (reactive) kern_ok <- dsep(g, R(t), X(t), c(za, R(t - 1)))   # assessment-conditioned kernel p(X_t | X_{t-1}, eta, R_{t-1} = 1)
  if (kern_ok && !reactive) {
    est <- if (length(Cobs)) "complete pairs, context-conditional kernel (covariate adjustment); context-marginal kernel by stabilized IPW on C" else "complete adjacent pairs, within-person (person intercepts)"
    push("transition kernel (Phi, Psi, contemporaneous network)", TRUE, est,
         sprintf("R_t _||_ X_t | {%s} and R_{t-1} _||_ X_t | {%s, R_t}", paste(za, collapse = ","), paste(za, collapse = ",")))
  } else if (reactive && kern_ok) {
    push("transition kernel (Phi, Psi, contemporaneous network)", TRUE,
         "complete adjacent pairs recover the assessment-conditioned kernel p(X_t | X_{t-1}, eta, R_{t-1} = 1): Phi and Psi are unaffected when reactivity shifts the level (additive), but the person intercept absorbs the shift",
         "reactivity edge R_{t-1} -> X_t declared; all complete pairs have R_{t-1} = 1")
  } else if (reactive) {
    push("transition kernel (Phi, Psi, contemporaneous network)", FALSE,
         "reactivity combined with a motif that makes R_t depend on X_t: sensitivity analysis required",
         "reactivity edge R_{t-1} -> X_t declared and R_t, X_t d-connected")
  } else {
    why <- if (!a) sprintf("R_t and X_t are d-connected given {%s}", paste(za, collapse = ",")) else
      sprintf("R_{t-1} and X_t are d-connected given {%s, R_t}", paste(za, collapse = ","))
    push("transition kernel (Phi, Psi, contemporaneous network)", FALSE,
         "sensitivity analysis (tilt profile) or a calibration design is required", why)
  }
  # --- probe rule
  if (g$probe) {
    zp <- c(X(t - 1), pl, Cobs)
    okp <- dsep(g, paste0("Z", t), X(t), zp)
    push("transition kernel from probe prompts", okp,
         if (okp) "complete pairs restricted to probe prompts (Z_t = 1); delta calibrated by the probe contrast" else "probe not independent of X_t",
         sprintf("Z_t _||_ X_t | {%s}; R_t = 1 whenever Z_t = 1", paste(zp, collapse = ",")))
  }
  # --- person mean
  om <- dsep(g, R(t), X(t), c(pl, Cobs))
  push("person mean via observed within-person mean", om,
       if (om) "mean of answered prompts, per person" else "biased",
       sprintf("R_t _||_ X_t | {%s}", paste(c(pl, Cobs), collapse = ",")))
  push("person mean via recovered dynamics", kern_ok && !reactive,
       if (kern_ok && !reactive) "mu_i = (I - Phi)^{-1} c_i from the recovered kernel" else
         if (reactive) "biased: the answered-pair intercept is the intercept of the assessed process (it absorbs the reactivity shift)" else "not available",
       if (reactive) "transition kernel recoverable and no reactivity" else "transition kernel recoverable")
  # --- between-person law
  any_mean <- om || (kern_ok && !reactive)
  push("between-person law (mu, Sigma_mu), person-weighted", any_mean,
       if (any_mean) paste0("one recovered mean per person (", if (om) "observed means" else "dynamics-recovered means", "), persons weighted equally; positivity assumed") else "not recoverable",
       "a person-mean estimator is available and P(R_t = 1 | eta) > 0")
  pw <- dsep(g, R(t), X(t), Cobs)
  push("between-person law, prompt-weighted (pooling answered prompts)", pw,
       if (pw) "pooled answered prompts" else "biased: response rate depends on the person's states",
       sprintf("R_t _||_ X_t | {%s}", if (length(Cobs)) paste(Cobs, collapse = ",") else "empty set"))
  # --- silence test validity (null holds under recoverable class)
  zs <- c(X(t - 1), pl, Cobs, R(t - 1), R(t + 1))
  st <- dsep(g, R(t), X(t + 1), zs)
  burden <- any(g$edges[, "from"] == R(t) & g$edges[, "to"] == R(t + 1))
  push("silence test (coefficient of R_t in X_{t+1} ~ X_{t-1} + R_t, within person)",
       NA, if (st) "null expected (test valid as a test of the recoverable class)" else
         if (kern_ok && !reactive) "null violated although the kernel is recoverable (collider at R_{t+1}); use the sensor-gap test or check for fatigue first" else
           if (burden) "non-null expected, but the burden collider at R_{t+1} works against the state-dependent signal: low power; prefer the sensor-gap test" else
           "non-null expected: silence is informative",
       sprintf("R_t _||_ X_{t+1} | {%s}", paste(zs, collapse = ",")))
  attr(rows, "silence_null") <- st
  if (g$sensor) {
    zg <- c(X(t - 1), pl, Cobs, R(t - 1))
    sg <- dsep(g, R(t), paste0("S", t), zg)
    push("sensor-gap test (coefficient of R_t in S_t ~ X_{t-1} + R_t, within person)", NA,
         if (sg) "null expected (test valid as a test of the recoverable class)" else
           if (kern_ok && !reactive) "null violated although the kernel is recoverable" else
             "non-null expected: delta can be calibrated from the sensor gap",
         sprintf("R_t _||_ S_t | {%s}", paste(zg, collapse = ",")))
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  class(out) <- c("recoverability", "data.frame")
  attr(out, "motifs") <- g$motifs
  attr(out, "silence_null") <- st
  out
}

#' @export
print.recoverability <- function(x, ...) {
  if (!all(c("estimand", "recoverable", "estimator", "condition") %in% names(x))) return(print.data.frame(x, ...))
  cat("Recoverability report for motifs:", paste(attr(x, "motifs"), collapse = " + "), "\n\n")
  tn <- function(v) {   # window indices -> t-notation (the focal prompt is the third of the window)
    v <- gsub("([XRCSZ])(\\d+)", "\\1_{\\2}", v)
    for (k in 1:8) v <- gsub(paste0("_\\{", k, "\\}"), if (k == 3) "_t" else if (k < 3) paste0("_{t-", 3 - k, "}") else paste0("_{t+", k - 3, "}"), v)
    v
  }
  for (i in seq_len(nrow(x))) {
    verdict <- if (is.na(x$recoverable[i])) "" else if (x$recoverable[i]) "[recoverable] " else "[NOT recoverable] "
    cat(sprintf("* %s\n    %s%s\n    condition: %s\n", x$estimand[i], verdict, x$estimator[i], tn(x$condition[i])))
  }
  invisible(x)
}

# Subsetting a report returns a plain data frame (the attributes that the print method needs are lost anyway)
#' @export
`[.recoverability` <- function(x, ...) {
  out <- NextMethod()
  if (is.data.frame(out)) class(out) <- "data.frame"
  out
}

#' Plot a dm-graph
#'
#' Draws the unrolled graph with prompts on the horizontal axis and node
#' types in rows (person-level latents on top, then context, states, sensor,
#' response indicators and probes). Structural edges are gray, edges added by
#' the declared motifs are black and dashed, so the plot remains readable in
#' grayscale. Latent nodes are drawn as open circles.
#'
#' @param x A \code{dm_graph}.
#' @param ... Further arguments passed to \code{plot}.
#' @return Invisibly, a numeric matrix with one row per node and the columns
#'   \code{x} and \code{y} (the plotting coordinates). Called for its side
#'   effect, the plot.
#' @examples
#' plot(dm_graph(c("M2", "M4"), sensor = TRUE))
#' plot(dm_graph("M5", context = "observed", context_persistent = TRUE))
#' @export
plot.dm_graph <- function(x, ...) {
  W <- x$window
  nodes <- x$nodes
  pos <- matrix(NA_real_, nrow = length(nodes), ncol = 2, dimnames = list(nodes, c("x", "y")))
  rowy <- c(C = 5.2, X = 4.2, S = 3.3, R = 2.2, Z = 1.2)
  for (v in nodes) {
    if (v %in% c("eta", "zeta", "U")) {
      pos[v, ] <- switch(v, eta = c(1.5, 6.6), U = c((W + 1) / 2, 7.4), zeta = c(W - 0.5, 6.6))
    } else {
      typ <- substr(v, 1, 1); k <- as.integer(substr(v, 2, nchar(v)))
      pos[v, ] <- c(k, rowy[[typ]])
    }
  }
  graphics::plot(NA, xlim = c(0.5, W + 0.5), ylim = c(0.6, 7.9), axes = FALSE, xlab = "prompt", ylab = "", ...)
  graphics::axis(1, at = seq_len(W), labels = ifelse(seq_len(W) == 3, "t", ifelse(seq_len(W) < 3, paste0("t-", 3 - seq_len(W)), paste0("t+", seq_len(W) - 3))))
  structural <- (x$edges[, "from"] == "eta" & grepl("^X", x$edges[, "to"])) |
    (grepl("^X", x$edges[, "from"]) & grepl("^X", x$edges[, "to"])) |
    (grepl("^C", x$edges[, "from"]) & grepl("^C", x$edges[, "to"])) |
    (grepl("^X", x$edges[, "from"]) & grepl("^S", x$edges[, "to"])) |
    (grepl("^Z", x$edges[, "from"]) & grepl("^R", x$edges[, "to"])) | x$edges[, "from"] == "U"
  rad <- 0.22   # node radius in user units (approximate; arrows stop at the node border)
  for (i in seq_len(nrow(x$edges))) {
    f <- x$edges[i, 1]; tt <- x$edges[i, 2]
    dx <- pos[tt, 1] - pos[f, 1]; dy <- pos[tt, 2] - pos[f, 2]; len <- sqrt(dx^2 + dy^2)
    if (len < 1e-8) next
    ux <- dx / len; uy <- dy / len
    graphics::arrows(pos[f, 1] + rad * ux, pos[f, 2] + rad * uy, pos[tt, 1] - rad * ux, pos[tt, 2] - rad * uy,
                     length = 0.07, col = if (structural[i]) "grey55" else "black", lty = if (structural[i]) 1 else 2,
                     lwd = if (structural[i]) 1 else 1.6)
  }
  lat <- nodes %in% x$latent
  graphics::symbols(pos[, 1], pos[, 2], circles = rep(rad, nrow(pos)), inches = FALSE, add = TRUE,
                    bg = ifelse(lat, "white", "grey90"), fg = "black", lwd = ifelse(lat, 1.2, 1))
  lab <- as.list(nodes)
  lab[nodes == "eta"] <- list(expression(eta)); lab[nodes == "zeta"] <- list(expression(zeta))
  for (i in seq_along(nodes)) graphics::text(pos[i, 1], pos[i, 2], labels = lab[[i]], cex = 0.75)
  graphics::title(main = paste("dm-graph:", paste(x$motifs, collapse = " + ")))
  graphics::legend("bottomleft", legend = c("structural edge", "motif edge", "latent node"), lty = c(1, 2, NA), pch = c(NA, NA, 1),
                   col = c("grey55", "black", "black"), lwd = c(1, 1.6, NA), bty = "n", cex = 0.8, horiz = TRUE)
  invisible(pos)
}
