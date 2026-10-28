coverage_background <- function() {
  data.frame(
    temperature = c(-2.2, -1.7, -1.2, -0.8, -0.2, 0.1, 0.5, 0.9,
                    1.2, 1.6, 2.0, 2.4),
    rainfall = c(1.7, 0.9, 1.4, 0.2, 0.8, -0.5, 0.1, -1.1,
                 -0.3, -1.6, -0.8, -2.0),
    elevation = c(-1.4, -0.7, -1.0, 0.2, -0.4, 0.9, 0.3, 1.5,
                  0.6, 1.9, 1.2, 2.3)
  )
}

coverage_occurrences <- function(n, years = rep(2000, n),
                                 identifiers = rep("survey", n)) {
  data.frame(
    species = rep("example", n),
    x = seq_len(n), y = seq_len(n), year = years,
    spatialUncertainty = 0, identifier = identifiers,
    stringsAsFactors = FALSE
  )
}

coverage_call <- function(dat, env, background = coverage_background(), ...,
                          periods = list(2000), bandwidth = diag(c(0.25, 0.20)),
                          gridSize = 41) {
  assessEnvCoverage(
    dat = dat, species = "species", x = "x", y = "y", year = "year",
    spatialUncertainty = "spatialUncertainty", identifier = "identifier",
    periods = periods, envDat = env, backgroundEnvDat = background,
    bandwidth = bandwidth, gridSize = gridSize, ...
  )
}

manual_kde_density <- function(evaluation, sample, bandwidth) {
  inverse <- solve(bandwidth)
  normalizer <- 2 * pi * sqrt(det(bandwidth))
  vapply(seq_len(nrow(evaluation)), function(index) {
    differences <- sweep(sample, 2, evaluation[index, ], "-")
    exponent <- rowSums((differences %*% inverse) * differences)
    mean(exp(-0.5 * exponent) / normalizer)
  }, numeric(1))
}

manual_hdr <- function(sample, background, bandwidth, bounds, gridSize, probs) {
  x <- seq(bounds$min[[1]], bounds$max[[1]], length.out = gridSize)
  y <- seq(bounds$min[[2]], bounds$max[[2]], length.out = gridSize)
  grid <- as.matrix(expand.grid(x, y))
  density <- manual_kde_density(grid, sample, bandwidth)
  area <- (x[[2]] - x[[1]]) * (y[[2]] - y[[1]])
  order_high <- order(density, decreasing = TRUE)
  cumulative <- cumsum(density[order_high] * area)
  cumulative <- cumulative / tail(cumulative, 1)
  thresholds <- vapply(probs, function(probability) {
    density[order_high[[which(cumulative >= probability)[[1]]]]]
  }, numeric(1))
  background_density <- manual_kde_density(background, sample, bandwidth)
  inside <- vapply(thresholds, function(threshold) {
    sum(background_density >= threshold)
  }, integer(1))
  list(thresholds = thresholds, inside = inside,
       coverage = inside / nrow(background), density = background_density)
}

test_that("background PCA is scaled and occurrence scores use that PCA", {
  background <- coverage_background()
  env <- background[c(2, 4, 5, 7, 9, 11), ]
  dat <- coverage_occurrences(nrow(env))
  result <- coverage_call(dat, env, background)
  expected <- stats::prcomp(background, center = TRUE, scale. = TRUE,
                            retx = FALSE)
  expected_scores <- stats::predict(expected, env)[, 1:2]

  expect_equal(result$metadata$pca$center, expected$center)
  expect_equal(result$metadata$pca$scale, expected$scale)
  expect_equal(result$metadata$pca$rotation, expected$rotation)
  expect_equal(
    unname(as.matrix(result$metadata$sampleScores[, c("xScore", "yScore")])),
    unname(expected_scores)
  )
  expect_equal(result$metadata$varianceExplained,
               expected$sdev^2 / sum(expected$sdev^2))
})

test_that("fixed KDE thresholds and membership match independent calculations", {
  background <- coverage_background()
  env <- background[c(1, 3, 4, 6, 8, 10, 12), ]
  dat <- coverage_occurrences(nrow(env))
  probabilities <- c(0.5, 0.8, 0.95)
  fixed <- matrix(c(0.28, 0.04, 0.04, 0.20), 2, 2)
  result <- coverage_call(
    dat, env, background, hdrProbs = probabilities,
    bandwidth = fixed, gridSize = 51
  )
  pca <- stats::prcomp(background, center = TRUE, scale. = TRUE, retx = FALSE)
  sample_scores <- stats::predict(pca, env)[, 1:2]
  background_scores <- stats::predict(pca, background)[, 1:2]
  expected <- manual_hdr(
    sample_scores, background_scores, fixed,
    result$metadata$gridBounds, 51, probabilities
  )

  expect_equal(result$data$hdrThreshold, expected$thresholds, tolerance = 2e-7)
  expect_equal(result$data$nBackgroundInside, expected$inside)
  expect_equal(result$data$backgroundCoverage, expected$coverage)
  expect_true(all(result$data$status == "ok"))
})

test_that("filtering preserves occurrence-environment row alignment", {
  background <- coverage_background()
  env <- background[1:8, ]
  env$temperature[[3]] <- NA
  dat <- coverage_occurrences(
    8, years = c(2000, 2000, 2000, 1990, 2000, 2000, 2000, 2000)
  )
  dat$spatialUncertainty[[2]] <- 100
  result <- suppressWarnings(coverage_call(
    dat, env, background, maxSpatUncertainty = 10
  ))
  kept <- c(1, 5, 6, 7, 8)
  pca <- stats::prcomp(background, center = TRUE, scale. = TRUE, retx = FALSE)
  expected <- stats::predict(pca, env[kept, ])[, 1:2]

  expect_equal(result$metadata$sampleScores$row, kept)
  expect_equal(
    unname(as.matrix(result$metadata$sampleScores[, c("xScore", "yScore")])),
    unname(expected)
  )
  expect_equal(unique(result$data$nSample), 5L)
})

test_that("coverage is monotone and detects restricted environmental sampling", {
  set.seed(107)
  background <- data.frame(
    v1 = rnorm(250),
    v2 = 0.55 * rnorm(250) + rnorm(250, sd = 0.75)
  )
  broad <- data.frame(v1 = rnorm(100), v2 = rnorm(100))
  restricted <- data.frame(
    v1 = rnorm(100, -1.4, 0.25),
    v2 = rnorm(100, -1.2, 0.25)
  )
  env <- rbind(broad, restricted)
  dat <- coverage_occurrences(200, years = c(rep(2000, 100), rep(2010, 100)))
  result <- coverage_call(
    dat, env, background, periods = list(2000, 2010),
    hdrProbs = c(0.5, 0.8, 0.95), bandwidth = diag(c(0.18, 0.18)),
    gridSize = 61
  )

  coverages <- split(result$data$backgroundCoverage, result$data$period)
  expect_true(all(vapply(coverages, function(values) {
    all(diff(values) >= -1e-12)
  }, logical(1))))
  broad95 <- subset(result$data, period == 1 & hdrProbability == 0.95)
  restricted95 <- subset(result$data, period == 2 & hdrProbability == 0.95)
  expect_gt(broad95$backgroundCoverage, restricted95$backgroundCoverage)
  expect_true(all(result$data$backgroundCoverage >= 0 &
                  result$data$backgroundCoverage <= 1))
})

test_that("periods and identifiers are assessed separately", {
  background <- coverage_background()
  env <- rbind(background[1:6, ], background[4:9, ],
               background[2:7, ], background[7:12, ])
  dat <- coverage_occurrences(
    24,
    years = rep(c(2000, 2010), each = 6, times = 2),
    identifiers = rep(c("alpha", "beta"), each = 12)
  )
  result <- coverage_call(
    dat, env, background, periods = list(2000, 2010), hdrProbs = c(0.5, 0.9)
  )

  combinations <- unique(result$data[, c("identifier", "period")])
  expect_equal(
    paste(combinations$identifier, combinations$period),
    c("alpha 1", "alpha 2", "beta 1", "beta 2")
  )
  expect_equal(result$data$nSample, rep(6L, 8))
})

test_that("fixed and pooled bandwidths are reproducible and common", {
  background <- coverage_background()
  env <- background[c(1:12, 2:9), ]
  dat <- coverage_occurrences(
    20, years = rep(c(2000, 2010), each = 10)
  )
  fixed <- matrix(c(0.3, 0.02, 0.02, 0.25), 2, 2)
  one <- coverage_call(dat, env, background, periods = list(2000, 2010),
                       bandwidth = fixed, bandwidthAdjust = 1.5)
  two <- coverage_call(dat, env, background, periods = list(2000, 2010),
                       bandwidth = fixed, bandwidthAdjust = 1.5)
  expect_equal(one$metadata$bandwidth, fixed * 1.5^2)
  expect_equal(one$data, two$data)

  pooled_one <- coverage_call(dat, env, background, periods = list(2000, 2010),
                              bandwidth = NULL, gridSize = 31)
  pooled_two <- coverage_call(dat, env, background, periods = list(2000, 2010),
                              bandwidth = NULL, gridSize = 31)
  expect_equal(pooled_one$metadata$bandwidth, pooled_two$metadata$bandwidth)
  expect_equal(pooled_one$data, pooled_two$data)
  expect_equal(pooled_one$metadata$bandwidthSource, "pooled_Hpi")
})

test_that("insufficient, degenerate, and empty samples have statuses", {
  background <- coverage_background()
  varied <- background[1:4, ]
  identical <- background[rep(6, 5), ]
  env <- rbind(varied, identical)
  dat <- coverage_occurrences(
    9,
    years = c(rep(2000, 4), rep(2000, 5)),
    identifiers = c(rep("few", 4), rep("same", 5))
  )
  result <- coverage_call(dat, env, background, periods = list(2000, 2010))
  statuses <- unique(result$data[, c("identifier", "period", "status")])

  expect_equal(
    statuses$status,
    c("insufficient_observations", "no_observations",
      "degenerate_sample", "no_observations")
  )
  expect_true(all(is.na(result$data$backgroundCoverage)))
  expect_s3_class(result$plot, "ggplot")
})

test_that("missing environmental values are reported and removed", {
  background <- coverage_background()
  background$rainfall[[1]] <- NA
  env <- coverage_background()[2:8, ]
  env$elevation[[2]] <- Inf
  dat <- coverage_occurrences(7)
  warnings <- character()
  result <- withCallingHandlers(
    coverage_call(dat, env, background),
    warning = function(condition) {
      warnings <<- c(warnings, conditionMessage(condition))
      invokeRestart("muffleWarning")
    }
  )
  expect_match(paste(warnings, collapse = " "), "background row")
  expect_match(paste(warnings, collapse = " "), "occurrence row")
  expect_equal(result$metadata$nBackgroundOriginal, 12L)
  expect_equal(result$metadata$nBackground, 11L)
  expect_equal(unique(result$data$nSample), 6L)
})

test_that("uncertain dates allocate each observation at most once", {
  background <- coverage_background()
  env <- rbind(background[1:5, ], background[6:10, ])
  dat <- coverage_occurrences(10, years = NA_real_)
  dat$start <- c(rep("1995", 5), rep("1999", 5))
  dat$end <- c(rep("1995", 5), rep("2003", 5))
  result <- suppressWarnings(coverage_call(
    dat, env, background, periods = list(1990:1999, 2000:2009),
    dateHandling = list(
      start = "start", end = "end", ambiguous = "maxOverlap"
    )
  ))

  expect_equal(unique(result$data$nSample[result$data$period == 1]), 5L)
  expect_equal(unique(result$data$nSample[result$data$period == 2]), 5L)
  expect_equal(sum(attr(result, "temporalDiagnostics")$periodCounts), 10L)
  expect_equal(nrow(result$metadata$sampleScores), 10L)
})

test_that("plot contours use numerical HDR levels and correct facets", {
  background <- coverage_background()
  env <- rbind(background[1:6, ], background[7:12, ])
  dat <- coverage_occurrences(12, years = rep(c(2000, 2010), each = 6))
  fixed <- diag(c(0.22, 0.22))
  result <- coverage_call(
    dat, env, background, periods = list(2000, 2010),
    hdrProbs = c(0.5, 0.9), bandwidth = fixed, gridSize = 51
  )
  contour_data <- result$plot$layers[[3]]$data
  expect_setequal(unique(contour_data$hdrProbability), c(0.5, 0.9))
  expect_setequal(as.character(unique(contour_data$period)), c("p1", "p2"))

  first <- contour_data[
    contour_data$period == "p1" & contour_data$hdrProbability == 0.5, , drop = FALSE
  ]
  sample <- as.matrix(subset(
    result$metadata$sampleScores, period == 1,
    select = c("xScore", "yScore")
  ))
  contour_density <- manual_kde_density(
    as.matrix(first[, c("xScore", "yScore")]), sample, fixed
  )
  threshold <- subset(
    result$data, period == 1 & hdrProbability == 0.5
  )$hdrThreshold
  expect_equal(stats::median(contour_density), threshold, tolerance = 3e-3)

  layout <- ggplot2::ggplot_build(result$plot)$layout$layout
  expect_setequal(as.character(layout$period), c("p1", "p2"))
  expect_equal(unique(as.character(layout$identifier)), "survey")
  original_data <- result$data
  result$plot <- result$plot + ggplot2::theme_void()
  expect_equal(result$data, original_data)
})

test_that("inputs, PCA rank, bandwidth, and grid settings are validated", {
  background <- coverage_background()
  env <- background[1:6, ]
  dat <- coverage_occurrences(6)
  expect_error(coverage_call(dat, env, hdrProbs = c(0, 0.5)), "hdrProbs")
  expect_error(coverage_call(dat, env, hdrProbs = c(0.5, 0.5)), "hdrProbs")
  expect_error(coverage_call(dat, env, bandwidthAdjust = 0), "bandwidthAdjust")
  expect_error(coverage_call(dat, env, gridSize = 10), "gridSize")
  expect_error(coverage_call(dat, env, bandwidth = diag(c(1, 0))), "bandwidth")
  expect_error(
    coverage_call(dat, env[, 2:1], background),
    "identical column names and order"
  )
  singular_background <- data.frame(
    temperature = 1:6, rainfall = 2 * (1:6)
  )
  expect_error(
    coverage_call(dat, env[, 1:2], singular_background),
    "not estimable"
  )
  expect_error(coverage_call(dat, env, xPC = 1, yPC = 1), "distinct")
})
