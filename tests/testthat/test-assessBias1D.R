bias1d_reference_data <- function() {
  population <- data.frame(
    environment = rep(1:4, c(20, 10, 40, 30)),
    sampled = 0L
  )
  sampled_rows <- c(1, 21:22, 31:33, 71:74)
  population$sampled[sampled_rows] <- 1L
  list(population = population, sampled_rows = sampled_rows)
}

bias1d_call <- function(population, ..., breaks = c(0.5, 1.5, 2.5, 3.5, 4.5)) {
  assessBias1D(
    pop = population, R = "sampled", x = "environment",
    breaks = breaks, RNames = "survey", returnMetrics = TRUE, ...
  )
}

test_that("identical binned distributions have zero discrepancy", {
  population <- data.frame(
    environment = rep(1:4, each = 2),
    sampled = rep(c(1, 0), 4)
  )
  result <- bias1d_call(population)

  expect_equal(result$metrics$totalVariation, 0)
  expect_equal(result$metrics$meanAbsoluteError, 0)
  expect_equal(result$metrics$nBins, 4L)
})

test_that("disjoint distributions have total variation one", {
  metric <- occAssess:::.bias1DMetricRow(
    id = "survey", comparison = "sample", weightType = NA_character_,
    sampleDistribution = c(1, 0, 0, 0),
    populationDistribution = c(0, 0, 0, 1),
    nBins = 4, nPopulation = 10, nSample = 10, status = "ok"
  )

  expect_equal(metric$totalVariation, 1)
  expect_equal(metric$meanAbsoluteError, 0.5)
})

test_that("explicit four-bin reference example gives TV 0.2 and MAE 0.1", {
  reference <- bias1d_reference_data()
  result <- bias1d_call(reference$population)

  expect_equal(result$metrics$totalVariation, 0.2)
  expect_equal(result$metrics$meanAbsoluteError, 0.1)
  expect_equal(result$metrics$meanAbsoluteError,
               2 * result$metrics$totalVariation / result$metrics$nBins)

  frequencies <- reshape2::dcast(
    result$data, variable ~ bin, value.var = "value"
  )
  expect_equal(as.numeric(frequencies[frequencies$variable == "population", -1]),
               c(0.2, 0.1, 0.4, 0.3))
  expect_equal(as.numeric(frequencies[frequencies$variable == "sample", -1]),
               c(0.1, 0.2, 0.3, 0.4))
})

test_that("all common bins are retained including empty bins", {
  population <- data.frame(
    environment = c(0.5, 0.5, 3.5, 3.5),
    sampled = c(1, 0, 1, 0)
  )
  result <- bias1d_call(population, breaks = 0:4)

  expect_setequal(unique(result$data$bin), 1:4)
  population_values <- subset(result$data, variable == "population")$value
  expect_equal(population_values, c(0.5, 0, 0, 0.5))
  expect_equal(sum(population_values), 1)
})

test_that("multiple indicators share population bins and metrics", {
  population <- data.frame(
    environment = rep(1:4, each = 2),
    early = c(1, 0, 1, 0, 1, 0, 1, 0),
    late = c(0, 0, 1, 1, 0, 0, 1, 1)
  )
  result <- assessBias1D(
    population, R = c("early", "late"), x = "environment", breaks = 4,
    RNames = c("early survey", "late survey"), returnMetrics = TRUE
  )

  expect_equal(result$metrics$id, c("early survey", "late survey"))
  expect_equal(result$metrics$nBins, c(4L, 4L))
  population_rows <- subset(result$data, variable == "population")
  by_id <- split(population_rows, population_rows$id)
  expect_equal(by_id[[1]]$value, by_id[[2]]$value)
  expect_equal(unname(vapply(
    split(result$data, interaction(result$data$id, result$data$variable)),
    function(rows) sum(rows$value), numeric(1)
  )),
               rep(1, 4))
})

test_that("metrics use proportions rather than population counts", {
  reference <- bias1d_reference_data()
  doubled <- reference$population[rep(seq_len(nrow(reference$population)), 2), ]
  doubled$sampled <- 0L
  sampled_by_bin <- c(1, 2, 3, 4)
  for (bin in 1:4) {
    rows <- which(doubled$environment == bin)
    doubled$sampled[rows[seq_len(sampled_by_bin[[bin]])]] <- 1L
  }

  result <- bias1d_call(doubled)
  expect_equal(result$metrics$totalVariation, 0.2)
  expect_equal(result$metrics$meanAbsoluteError, 0.1)
})

test_that("population-row weights can exactly reproduce the population", {
  reference <- bias1d_reference_data()
  calibrated <- numeric(nrow(reference$population))
  calibrated[reference$sampled_rows] <- rep(
    c(0.2, 0.05, 4 / 30, 0.075), c(1, 2, 3, 4)
  )
  result <- bias1d_call(
    reference$population, weights = calibrated, WNames = "calibrated"
  )

  unweighted <- subset(result$metrics, comparison == "sample")
  weighted <- subset(result$metrics, comparison == "weighted_sample")
  expect_equal(unweighted$totalVariation, 0.2)
  expect_equal(unweighted$meanAbsoluteError, 0.1)
  expect_equal(weighted$totalVariation, 0, tolerance = 1e-14)
  expect_equal(weighted$meanAbsoluteError, 0, tolerance = 1e-14)
  expect_equal(sum(result$data$value[result$data$variable == "weighted_sample"]),
               1)
})

test_that("multiple weighting schemes add one weighted metric each", {
  reference <- bias1d_reference_data()
  calibrated <- numeric(nrow(reference$population))
  calibrated[reference$sampled_rows] <- rep(
    c(0.2, 0.05, 4 / 30, 0.075), c(1, 2, 3, 4)
  )
  equal <- as.numeric(reference$population$sampled == 1)
  result <- bias1d_call(
    reference$population,
    weights = list(calibrated, equal),
    WNames = c("calibrated", "equal")
  )

  expect_equal(result$metrics$comparison,
               c("sample", "weighted_sample", "weighted_sample"))
  expect_equal(result$metrics$weightType,
               c(NA, "calibrated", "equal"))
  expect_equal(result$metrics$totalVariation, c(0.2, 0, 0.2),
               tolerance = 1e-14)
  expect_equal(sum(result$metrics$comparison == "sample"), 1)
})

test_that("legacy nested sample-aligned weights remain supported", {
  population <- data.frame(
    environment = rep(1:4, each = 2),
    early = c(1, 0, 1, 0, 1, 0, 1, 0),
    late = c(1, 1, 0, 0, 1, 1, 0, 0)
  )
  nested <- list(list(c(1, 2, 3, 4), c(4, 3, 2, 1)))
  result <- assessBias1D(
    population, R = c("early", "late"), x = "environment", breaks = 4,
    RNames = c("early", "late"), weights = nested, WNames = "legacy",
    returnMetrics = TRUE
  )

  expect_equal(nrow(result$metrics), 4L)
  expect_equal(subset(result$metrics, comparison == "weighted_sample")$weightType,
               c("legacy", "legacy"))
  expect_true(all(subset(result$metrics,
                         comparison == "weighted_sample")$status == "ok"))
})

test_that("weights stay aligned when missing x values are removed", {
  population <- data.frame(
    environment = c(1, NA, 2, 3, 4, Inf),
    sampled = c(1, 1, 1, 0, 1, 0)
  )
  full_weights <- c(1, 100, 2, 0, 4, 0)
  full <- bias1d_call(
    population, weights = full_weights, WNames = "full",
    breaks = c(0.5, 1.5, 2.5, 3.5, 4.5)
  )
  nested <- bias1d_call(
    population, weights = list(list(c(1, 100, 2, 4))),
    WNames = "nested", breaks = c(0.5, 1.5, 2.5, 3.5, 4.5)
  )

  full_distribution <- subset(full$data, variable == "weighted_sample")$value
  nested_distribution <- subset(nested$data,
                                variable == "weighted_sample")$value
  expect_equal(full_distribution, c(1, 2, 0, 4) / 7)
  expect_equal(nested_distribution, full_distribution)
  expect_equal(unique(full$metrics$nPopulation), 4L)
  expect_equal(unique(full$metrics$nSample), 3L)
})

test_that("invalid and zero-total weights are rejected", {
  reference <- bias1d_reference_data()
  expect_error(
    bias1d_call(reference$population, weights = rep(1, 3), WNames = "bad"),
    "nrow\\(pop\\)"
  )
  expect_error(
    bias1d_call(reference$population,
                weights = c(rep(1, 99), -1), WNames = "bad"),
    "non-negative"
  )
  expect_error(
    bias1d_call(reference$population,
                weights = c(rep(1, 99), Inf), WNames = "bad"),
    "finite"
  )
  expect_error(
    bias1d_call(reference$population,
                weights = numeric(100), WNames = "zero"),
    "positive total weight"
  )
  expect_error(
    bias1d_call(reference$population, weights = rep(1, 100)),
    "WNames"
  )
})

test_that("empty and unusable samples return explicit unavailable metrics", {
  population <- data.frame(
    environment = 1:4,
    sampled = c(0, 0, NA, 0)
  )
  result <- bias1d_call(population)

  expect_equal(result$metrics$status, "no_eligible_sample")
  expect_true(is.na(result$metrics$totalVariation))
  expect_true(is.na(result$metrics$meanAbsoluteError))
  expect_true(all(is.na(subset(result$data, variable == "sample")$value)))
})

test_that("constant x, explicit boundaries, and boundary values are handled", {
  constant <- data.frame(environment = rep(5, 6), sampled = c(1, 0, 1, 0, 1, 0))
  constant_result <- bias1d_call(constant, breaks = 4)
  expect_equal(constant_result$metrics$totalVariation, 0)
  expect_equal(length(unique(constant_result$data$bin)), 4L)

  boundary <- data.frame(
    environment = c(0, 1, 2),
    sampled = c(1, 1, 1)
  )
  boundary_result <- bias1d_call(boundary, breaks = c(0, 1, 2))
  expect_equal(subset(boundary_result$data,
                      variable == "population")$value, c(2 / 3, 1 / 3))
})

test_that("invalid inputs fail clearly", {
  reference <- bias1d_reference_data()$population
  expect_error(bias1d_call(reference[FALSE, ]), "at least one row")
  no_x <- transform(reference, environment = NA_real_)
  expect_error(bias1d_call(no_x), "no finite values")
  invalid_indicator <- reference
  invalid_indicator$sampled[[1]] <- 2
  expect_error(bias1d_call(invalid_indicator), "only 0, 1, or NA")
  expect_error(bias1d_call(reference, breaks = 0), "positive integer")
  expect_error(bias1d_call(reference, breaks = c(0, 2, 1, 5)),
               "strictly increasing")
  expect_error(bias1d_call(reference, breaks = c(2, 5)), "cover every")
  expect_error(
    assessBias1D(transform(reference, sampled2 = sampled),
                 c("sampled", "sampled2"), "environment", breaks = 4,
                 RNames = c("duplicate", "duplicate")),
    "one unique"
  )
  expect_error(
    bias1d_call(reference, weights = list(rep(1, 100), rep(2, 100)),
                WNames = c("duplicate", "duplicate")),
    "uniquely name"
  )
})

test_that("default interface and optional metrics structure are preserved", {
  reference <- bias1d_reference_data()$population
  default <- assessBias1D(
    reference, R = "sampled", x = "environment", breaks = 4,
    RNames = "survey"
  )
  extended <- bias1d_call(reference, breaks = 4)

  expect_named(default, c("plot", "data"), ignore.order = FALSE)
  expect_named(extended, c("plot", "data", "metrics"), ignore.order = FALSE)
  expect_s3_class(default$plot, "ggplot")
  expect_equal(default$data, extended$data)
  expect_named(
    extended$metrics,
    c("id", "comparison", "weightType", "totalVariation",
      "meanAbsoluteError", "nBins", "nPopulation", "nSample", "status")
  )
  expect_true(all(extended$metrics$totalVariation >= 0 &
                  extended$metrics$totalVariation <= 1))
})
