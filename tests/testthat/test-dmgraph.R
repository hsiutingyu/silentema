test_that("dm_graph builds the expected edges", {
  g <- dm_graph("M0")
  expect_s3_class(g, "dm_graph")
  expect_true(all(c("eta", "X1", "R1") %in% g$nodes))
  expect_equal(sum(g$edges[, "to"] == "R3"), 0)
  g1 <- dm_graph("M1")
  expect_true(any(g1$edges[, "from"] == "X2" & g1$edges[, "to"] == "R3"))
  g2 <- dm_graph("M2")
  expect_true(any(g2$edges[, "from"] == "X3" & g2$edges[, "to"] == "R3"))
  g3 <- dm_graph("M3")
  expect_true(any(g3$edges[, "from"] == "R2" & g3$edges[, "to"] == "R3"))
  g4 <- dm_graph("M4")
  expect_true(all(c("U", "zeta") %in% g4$nodes))
  expect_message(g5 <- dm_graph("M5"), "latent context")
  expect_equal(g5$context, "latent")
  g6 <- dm_graph("M6")
  expect_true(any(g6$edges[, "from"] == "R2" & g6$edges[, "to"] == "X3"))
  expect_error(dm_graph("M9"))
  expect_message(g7 <- dm_graph("M0", context = "observed"), "M5 added")
  expect_true("M5" %in% g7$motifs)
  r6 <- recoverability(g6)
  expect_true(r6$recoverable[grepl("transition kernel", r6$estimand)][1])
  expect_true(r6$recoverable[grepl("observed within-person mean", r6$estimand)])
  expect_false(r6$recoverable[grepl("recovered dynamics", r6$estimand)])
  r23 <- recoverability(dm_graph(c("M2", "M3")))
  expect_match(r23$estimator[grepl("silence test", r23$estimand)], "low power")
  expect_output(print(r23), "X_\\{t-1\\}")
  pdf(NULL); on.exit(dev.off())
  expect_silent(plot(dm_graph(c("M2", "M4", "M5"), sensor = TRUE, probe = TRUE, context = "observed")))
})

test_that("d-separation gives textbook answers on small graphs", {
  # chain A -> B -> C
  g <- list(nodes = c("A", "B", "C"), edges = cbind(from = c("A", "B"), to = c("B", "C")))
  expect_false(dsep(g, "A", "C"))
  expect_true(dsep(g, "A", "C", "B"))
  # fork A <- B -> C
  g <- list(nodes = c("A", "B", "C"), edges = cbind(from = c("B", "B"), to = c("A", "C")))
  expect_false(dsep(g, "A", "C"))
  expect_true(dsep(g, "A", "C", "B"))
  # collider A -> B <- C
  g <- list(nodes = c("A", "B", "C"), edges = cbind(from = c("A", "C"), to = c("B", "B")))
  expect_true(dsep(g, "A", "C"))
  expect_false(dsep(g, "A", "C", "B"))
  # descendant of a collider opens the path
  g <- list(nodes = c("A", "B", "C", "D"), edges = cbind(from = c("A", "C", "B"), to = c("B", "B", "D")))
  expect_true(dsep(g, "A", "C"))
  expect_false(dsep(g, "A", "C", "D"))
})

test_that("d-separation agrees with linear-Gaussian partial correlations on dm-graphs", {
  set.seed(1)
  for (ms in list("M0", "M1", "M2", "M3", "M4", "M5", "M6", c("M1", "M3"), c("M1", "M4"), c("M2", "M4"))) {
    g <- suppressMessages(dm_graph(ms, sensor = TRUE))
    nodes <- g$nodes; k <- length(nodes)
    B <- matrix(0, k, k, dimnames = list(nodes, nodes))
    for (i in seq_len(nrow(g$edges))) B[g$edges[i, "to"], g$edges[i, "from"]] <- runif(1, 0.3, 0.9) * sample(c(-1, 1), 1)
    S <- solve(diag(k) - B) %*% diag(runif(k, 0.5, 1.5)) %*% t(solve(diag(k) - B))
    dimnames(S) <- list(nodes, nodes)
    pl <- intersect(c("eta", "zeta"), nodes)
    queries <- list(list("R3", "X3", c("X2", pl)), list("R2", "X3", c("X2", pl, "R3")), list("R3", "X3", pl),
                    list("R3", "X3", character(0)), list("R3", "X4", c("X2", pl, "R2", "R4")), list("R3", "S3", c("X2", pl, "R2")))
    for (q in queries) {
      idx <- c(q[[1]], q[[2]], q[[3]]); K <- solve(S[idx, idx]); pc <- -K[1, 2] / sqrt(K[1, 1] * K[2, 2])
      expect_equal(dsep(g, q[[1]], q[[2]], q[[3]]), abs(pc) < 1e-10, info = paste(paste(ms, collapse = "+"), q[[1]], q[[2]]))
    }
  }
})

test_that("recoverability verdicts match the hand-derived table (Study 4)", {
  expected <- list(M0 = c(TRUE, TRUE, TRUE, TRUE, TRUE), M1 = c(TRUE, FALSE, TRUE, TRUE, FALSE),
                   M2 = c(FALSE, FALSE, FALSE, FALSE, FALSE), M3 = c(TRUE, TRUE, TRUE, TRUE, TRUE),
                   M4 = c(TRUE, TRUE, TRUE, TRUE, FALSE), M5 = c(FALSE, FALSE, FALSE, FALSE, FALSE),
                   M6 = c(TRUE, TRUE, FALSE, TRUE, TRUE), "M1+M3" = c(TRUE, FALSE, TRUE, TRUE, FALSE),
                   "M1+M4" = c(TRUE, FALSE, TRUE, TRUE, FALSE), "M3+M4" = c(TRUE, TRUE, TRUE, TRUE, FALSE),
                   "M2+M4" = c(FALSE, FALSE, FALSE, FALSE, FALSE), "M1+M3+M4" = c(TRUE, FALSE, TRUE, TRUE, FALSE))
  silence_expected <- c(M0 = TRUE, M1 = TRUE, M2 = FALSE, M3 = TRUE, M4 = TRUE, M5 = FALSE, M6 = FALSE,
                        "M1+M3" = FALSE, "M1+M4" = TRUE, "M3+M4" = TRUE, "M2+M4" = FALSE, "M1+M3+M4" = FALSE)
  for (nm in names(expected)) {
    r <- recoverability(suppressMessages(dm_graph(strsplit(nm, "\\+")[[1]])))
    expect_equal(r$recoverable[1:5], expected[[nm]], info = nm)
    expect_equal(attr(r, "silence_null"), unname(silence_expected[nm]), info = nm)
  }
  r5 <- recoverability(dm_graph("M5", context = "observed"))
  expect_true(r5$recoverable[1])
  r5p <- recoverability(dm_graph("M5", context = "observed", context_persistent = TRUE))
  expect_true(r5p$recoverable[1])
  rp <- recoverability(dm_graph("M2", probe = TRUE))
  expect_true(rp$recoverable[rp$estimand == "transition kernel from probe prompts"])
  rs <- recoverability(dm_graph(c("M1", "M3"), sensor = TRUE))
  expect_true(grepl("null expected", rs$estimator[grepl("sensor-gap", rs$estimand)]))
  expect_output(print(r5), "Recoverability report")
})

test_that("missingness declaration writes a template", {
  g <- dm_graph(c("M2", "M4"), sensor = TRUE)
  tf <- tempfile(fileext = ".md")
  missingness_declaration(g, file = tf)
  txt <- readLines(tf)
  expect_true(any(grepl("Motifs: M2 \\+ M4", txt)))
  expect_true(any(grepl("NOT recoverable", txt)))
})

test_that("subsetting a report returns a plain data frame and plot honors a user title", {
  r <- recoverability(dm_graph("M2", sensor = TRUE))
  sub <- r[, c("estimand", "recoverable")]
  expect_identical(class(sub), "data.frame")
  expect_equal(nrow(sub), nrow(r))
  expect_output(print(r[1:2, ]), "estimand")
  pdf(NULL); on.exit(dev.off())
  expect_silent(plot(dm_graph("M1"), main = "custom title"))
  expect_error(recoverability(dm_graph("M1"), t = 2), "t must satisfy")
})
