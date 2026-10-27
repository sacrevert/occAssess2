make_env_bias_data <- function() {
  data.frame(
    species = rep("species-a", 6),
    x = seq_len(6),
    y = seq_len(6) + 10,
    spatialUncertainty = rep(1, 6),
    year = c(2000, 2001, 2003, 2000, 2001, 2003),
    identifier = rep(c("source-a", "source-b"), each = 3)
  )
}

run_env_bias <- function(dat, envDat, periods, ...) {
  assessEnvBias(
    dat = dat,
    species = "species",
    x = "x",
    y = "y",
    spatialUncertainty = "spatialUncertainty",
    year = "year",
    identifier = "identifier",
    periods = periods,
    envDat = envDat,
    ...
  )
}

test_that("uncertainty filtering preserves occurrence-environment pairing", {
  dat <- make_env_bias_data()[1:4, ]
  dat$identifier <- "source-a"
  dat$year <- 2000:2003
  dat$spatialUncertainty <- c(1, 100, 2, NA)
  envDat <- data.frame(
    temperature = c(1, 100, 4, 1000),
    precipitation = c(2, 200, 9, 2000)
  )

  result <- run_env_bias(
    dat,
    envDat,
    periods = list(2000:2001, 2002:2003),
    maxSpatUncertainty = 10
  )$data

  expected_pca <- stats::prcomp(envDat[c(1, 3), , drop = FALSE])

  expect_equal(nrow(result), 2)
  expect_equal(result$Period, c("p1", "p2"))
  expect_equal(result$identifier, rep("source-a", 2))
  expect_equal(
    unname(as.matrix(result[, c("scores.PC1", "scores.PC2")])),
    unname(expected_pca$x)
  )
  expect_equal(result$xVar, rep(summary(expected_pca)$importance[2, 1] * 100, 2))
  expect_equal(result$yVar, rep(summary(expected_pca)$importance[2, 2] * 100, 2))
})

test_that("unfiltered records retain identifiers and periods in PCA output", {
  dat <- make_env_bias_data()
  envDat <- data.frame(
    temperature = c(1, 2, 5, 10, 14, 13),
    precipitation = c(8, 3, 6, 2, 9, 15)
  )

  result <- run_env_bias(
    dat,
    envDat,
    periods = list(2000:2001, 2002:2003)
  )$data

  expected_a <- stats::prcomp(envDat[1:3, , drop = FALSE])
  expected_b <- stats::prcomp(envDat[4:6, , drop = FALSE])

  expect_equal(result$Period, rep(c("p1", "p1", "p2"), 2))
  expect_equal(result$identifier, rep(c("source-a", "source-b"), each = 3))
  expect_equal(
    unname(as.matrix(result[, c("scores.PC1", "scores.PC2")])),
    unname(rbind(expected_a$x, expected_b$x))
  )
})

test_that("missing and out-of-period records drop their paired environment rows", {
  dat <- make_env_bias_data()
  dat$identifier <- "source-a"
  dat$year <- c(2000, NA, 2001, 1999, 2001, 2003)
  envDat <- data.frame(
    temperature = c(1, 20, 3, 40, 5, 6),
    precipitation = c(8, 20, 2, 40, NA, 9)
  )

  expect_warning(
    result <- run_env_bias(
      dat,
      envDat,
      periods = list(2000:2001, 2003:2004)
    )$data,
    "Removing data without a specified year"
  )

  expected_pca <- stats::prcomp(envDat[c(1, 3, 6), , drop = FALSE])

  expect_equal(result$Period, c("p1", "p1", "p2"))
  expect_equal(
    unname(as.matrix(result[, c("scores.PC1", "scores.PC2")])),
    unname(expected_pca$x)
  )
})

test_that("background PCA projects paired records for every identifier", {
  dat <- make_env_bias_data()[c(1, 2, 4, 5), ]
  envDat <- data.frame(
    temperature = c(2, 5, 11, 13),
    precipitation = c(9, 3, 4, 12)
  )
  backgroundEnvDat <- data.frame(
    temperature = c(0, 3, 7, 12, 15),
    precipitation = c(1, 8, 4, 10, 14)
  )

  result <- run_env_bias(
    dat,
    envDat,
    periods = list(2000:2001),
    backgroundEnvDat = backgroundEnvDat
  )$data

  expected_pca <- stats::prcomp(backgroundEnvDat)
  expected_a <- stats::predict(expected_pca, envDat[1:2, , drop = FALSE])
  expected_b <- stats::predict(expected_pca, envDat[3:4, , drop = FALSE])
  expected_scores <- rbind(
    expected_a,
    expected_pca$x,
    expected_b,
    expected_pca$x
  )

  expect_equal(nrow(result), 14)
  expect_equal(result$identifier, rep(c("source-a", "source-b"), each = 7))
  expect_equal(
    result$Period,
    rep(c("p1", "p1", rep("background", 5)), 2)
  )
  expect_equal(
    unname(as.matrix(result[, c("scores.PC1", "scores.PC2")])),
    unname(expected_scores)
  )
  expect_equal(
    unique(result$xVar),
    summary(expected_pca)$importance[2, 1] * 100
  )
  expect_equal(
    unique(result$yVar),
    summary(expected_pca)$importance[2, 2] * 100
  )
})

test_that("environmental data must have one row per occurrence", {
  dat <- make_env_bias_data()
  envDat <- data.frame(
    temperature = seq_len(5),
    precipitation = seq_len(5)
  )

  expect_error(
    run_env_bias(dat, envDat, periods = list(2000:2003)),
    "nrow of environmental data does not equal nrow of species occurrence data",
    fixed = TRUE
  )
})
