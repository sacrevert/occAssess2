rarity_records <- function(identifier, year, occupied, counts,
                           species = paste0("sp", seq_along(occupied)),
                           x_offset = 0, y_offset = 0, res = 0.5) {
  rows <- lapply(seq_along(species), function(i) {
    unique_x <- x_offset + (seq_len(occupied[[i]]) - 1) * res + 0.1
    data.frame(
      species = species[[i]],
      x = c(unique_x, rep(unique_x[[1]], counts[[i]] - occupied[[i]])),
      y = y_offset + 0.1,
      year = year,
      spatialUncertainty = 0,
      identifier = identifier,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

rarity_call <- function(dat, ..., periods = list(2000),
                        prevPerPeriod = TRUE, res = 0.5) {
  assessRarityBias(
    dat = dat, species = "species", x = "x", y = "y", year = "year",
    spatialUncertainty = "spatialUncertainty", identifier = "identifier",
    periods = periods, res = res, prevPerPeriod = prevPerPeriod, ...
  )
}

test_that("species counts, OLS diagnostics, and influence values are numerical", {
  occupied <- c(1, 2, 3, 4, 5, 6)
  counts <- c(2, 4, 6, 8, 30, 7)
  dat <- rarity_records("survey", 2000, occupied, counts)
  result <- rarity_call(dat, returnDiagnostics = TRUE)

  observed <- result$speciesDiagnostics
  expected_model <- stats::lm(counts ~ occupied)

  expect_equal(observed$occupiedCells, occupied)
  expect_equal(observed$recordCount, counts)
  expect_equal(observed$fittedCount, unname(stats::fitted(expected_model)))
  expect_equal(observed$residual, counts - unname(stats::fitted(expected_model)))
  expect_equal(observed$standardizedResidual,
               unname(stats::rstandard(expected_model)))
  expect_equal(observed$leverage, unname(stats::hatvalues(expected_model)))
  expect_equal(observed$cooksDistance,
               unname(stats::cooks.distance(expected_model)))
  expect_gt(observed$residual[observed$species == "sp5"], 0)
  expect_lt(observed$residual[observed$species == "sp6"], 0)

  regression <- result$regressionDiagnostics
  expect_equal(regression$intercept, unname(stats::coef(expected_model)[[1]]))
  expect_equal(regression$slope, unname(stats::coef(expected_model)[[2]]))
  expect_equal(regression$rSquared, summary(expected_model)$r.squared)
  expect_equal(regression$correlation, stats::cor(counts, occupied))
  expect_equal(regression$residualRMSE,
               sqrt(mean(stats::residuals(expected_model)^2)))
  expect_equal(regression$status, "ok")
})

test_that("metric selects the established index without changing the model", {
  occupied <- c(1, 2, 3, 4, 5, 6)
  counts <- c(2, 4, 6, 8, 30, 7)
  dat <- rarity_records("survey", 2000, occupied, counts)

  default <- rarity_call(dat)
  r2 <- rarity_call(dat, metric = "r2", returnDiagnostics = TRUE)
  correlation <- rarity_call(dat, metric = "cor", returnDiagnostics = TRUE)

  expect_named(default, c("plot", "data"), ignore.order = FALSE)
  expect_s3_class(default$plot, "ggplot")
  expect_named(default$data, c("period", "id", "index"))
  expect_equal(r2$data$index, summary(stats::lm(counts ~ occupied))$r.squared)
  expect_equal(correlation$data$index, stats::cor(counts, occupied))
  expect_equal(r2$speciesDiagnostics, correlation$speciesDiagnostics)
})

test_that("identifiers and periods remain separate", {
  first <- rbind(
    rarity_records("alpha", 2000, 1:6, c(2, 4, 6, 8, 10, 12)),
    rarity_records("alpha", 2010, 1:6, c(3, 5, 7, 9, 11, 13), x_offset = 10)
  )
  second <- rbind(
    rarity_records("beta", 2000, 1:6, c(7, 7, 8, 9, 10, 12), y_offset = 10),
    rarity_records("beta", 2010, 1:6, c(2, 3, 5, 8, 13, 21),
                   x_offset = 10, y_offset = 10)
  )
  result <- rarity_call(
    rbind(first, second), periods = list(2000, 2010),
    plotRegressions = TRUE
  )

  expect_equal(
    result$regressionDiagnostics[, c("identifier", "period", "nSpecies")],
    data.frame(
      identifier = c("alpha", "alpha", "beta", "beta"),
      period = c(1L, 2L, 1L, 2L), nSpecies = rep(6L, 4)
    )
  )
  expect_equal(nrow(result$speciesDiagnostics), 24)
  expect_named(result$regressionPlots, c("alpha", "beta"))
})

test_that("prevPerPeriod changes only the occupied-cell window", {
  period1 <- rarity_records("survey", 2000, 1:6, c(2, 3, 4, 5, 6, 7))
  period2 <- rarity_records(
    "survey", 2010, rep(1, 6), rep(1, 6), x_offset = 10
  )
  dat <- rbind(period1, period2)

  per_period <- rarity_call(
    dat, periods = list(2000, 2010), prevPerPeriod = TRUE,
    returnDiagnostics = TRUE
  )
  all_periods <- rarity_call(
    dat, periods = list(2000, 2010), prevPerPeriod = FALSE,
    returnDiagnostics = TRUE
  )

  p1 <- per_period$speciesDiagnostics$period == 1
  a1 <- all_periods$speciesDiagnostics$period == 1
  expect_equal(per_period$speciesDiagnostics$recordCount[p1], 2:7)
  expect_equal(all_periods$speciesDiagnostics$recordCount[a1], 2:7)
  expect_equal(per_period$speciesDiagnostics$occupiedCells[p1], 1:6)
  expect_equal(all_periods$speciesDiagnostics$occupiedCells[a1], 2:7)
})

test_that("missing identities and coordinates have explicit treatment", {
  dat <- rarity_records("survey", 2000, 1:6, c(2, 3, 4, 5, 6, 7))
  dat <- rbind(
    dat,
    data.frame(
      species = c(NA, "sp1"), x = c(1, NA), y = c(1, Inf), year = 2000,
      spatialUncertainty = 0, identifier = "survey"
    )
  )

  warnings <- character()
  result <- withCallingHandlers(
    rarity_call(dat, returnDiagnostics = TRUE),
    warning = function(condition) {
      warnings <<- c(warnings, conditionMessage(condition))
      invokeRestart("muffleWarning")
    }
  )
  expect_match(paste(warnings, collapse = " "), "without a species identity")
  expect_match(paste(warnings, collapse = " "),
               "cannot contribute to occupied-cell counts")
  sp1 <- result$speciesDiagnostics[result$speciesDiagnostics$species == "sp1", ]
  expect_equal(sp1$recordCount, 3L)
  expect_equal(sp1$occupiedCells, 1L)
})

test_that("insufficient and degenerate combinations retain observed counts", {
  insufficient <- rarity_records("few", 2000, 1:4, c(2, 3, 4, 5))
  expect_warning(
    few <- rarity_call(insufficient, returnDiagnostics = TRUE),
    "Fewer than five species"
  )
  expect_equal(few$regressionDiagnostics$status, "insufficient_species")
  expect_true(all(is.na(few$speciesDiagnostics$fittedCount)))
  expect_equal(few$speciesDiagnostics$recordCount, c(2L, 3L, 4L, 5L))

  constant_x <- rarity_records("constant-x", 2000, rep(1, 5), 2:6)
  singular <- rarity_call(constant_x, returnDiagnostics = TRUE)
  expect_equal(singular$regressionDiagnostics$status,
               "constant_occupied_cells")
  expect_true(is.na(singular$data$index))
  expect_true(all(is.na(singular$speciesDiagnostics$residual)))

  constant_y <- rarity_records("constant-y", 2000, 1:5, rep(6, 5))
  flat <- rarity_call(constant_y, returnDiagnostics = TRUE)
  expect_equal(flat$regressionDiagnostics$status, "constant_record_count")
  expect_true(is.na(flat$data$index))
  expect_equal(flat$speciesDiagnostics$fittedCount, rep(6, 5), tolerance = 1e-12)
  expect_true(all(is.na(flat$speciesDiagnostics$standardizedResidual)))
  expect_true(all(is.na(flat$speciesDiagnostics$cooksDistance)))
})

test_that("empty periods and all-missing species receive statuses", {
  dat <- rarity_records("survey", 2000, 1:5, 2:6)
  result <- suppressWarnings(rarity_call(
    dat, periods = list(2000, 2010), returnDiagnostics = TRUE
  ))
  expect_equal(result$regressionDiagnostics$status,
               c("ok", "no_species"))
  expect_equal(result$regressionDiagnostics$nSpecies, c(5L, 0L))

  dat$species <- NA_character_
  all_missing <- suppressWarnings(rarity_call(
    dat, returnDiagnostics = TRUE
  ))
  expect_equal(all_missing$regressionDiagnostics$status, "no_species")
  expect_equal(nrow(all_missing$speciesDiagnostics), 0)
})

test_that("fractional zero-aligned cells and boundaries are deterministic", {
  dat <- rarity_records("survey", 2000, 1:5, rep(5, 5))
  dat$x[dat$species == "sp1"] <- c(-0.5, -0.1, 0, 0.49, 0.5)
  result <- rarity_call(dat, returnDiagnostics = TRUE)
  sp1 <- result$speciesDiagnostics[result$speciesDiagnostics$species == "sp1", ]

  # At res 0.5 these coordinates independently map to bins -1, -1, 0, 0, 1.
  expect_equal(sp1$occupiedCells, 3L)
})

test_that("uncertain observations contribute to no more than one period", {
  dat <- rarity_records("survey", 2000, 1:5, rep(5, 5))
  dat$start <- as.character(dat$year)
  dat$end <- as.character(dat$year)
  crossing <- dat[1, ]
  crossing$year <- NA
  crossing$start <- "1999"
  crossing$end <- "2003"
  dat <- rbind(dat, crossing)

  result <- suppressWarnings(rarity_call(
    dat, periods = list(1990:1999, 2000:2009),
    dateHandling = list(
      start = "start", end = "end", ambiguous = "maxOverlap"
    ),
    returnDiagnostics = TRUE
  ))
  expect_equal(sum(result$speciesDiagnostics$recordCount), nrow(dat))
  expect_equal(attr(result, "temporalDiagnostics")$periodCounts,
               c(p1 = 0L, p2 = nrow(dat)))
})

test_that("regression plots use supplied fits, facets, and restrained labels", {
  occupied <- 1:6
  counts <- c(2, 4, 6, 8, 30, 7)
  dat <- rbind(
    rarity_records("survey", 2000, occupied, counts),
    rarity_records("survey", 2010, occupied, counts + c(1, 0, 2, 1, 0, 3),
                   x_offset = 10)
  )
  result <- rarity_call(
    dat, periods = list(2000, 2010),
    plotRegressions = TRUE, labelTopResiduals = 2
  )

  expect_named(
    result,
    c("plot", "data", "speciesDiagnostics", "regressionDiagnostics",
      "regressionPlots")
  )
  expect_named(result$regressionPlots, "survey")
  expect_s3_class(result$regressionPlots$survey, "ggplot")
  built <- ggplot2::ggplot_build(result$regressionPlots$survey)
  expect_equal(nrow(built$data[[1]]), 12)
  expect_equal(built$data[[2]]$intercept,
               result$regressionDiagnostics$intercept)
  expect_equal(built$data[[2]]$slope,
               result$regressionDiagnostics$slope)
  expect_equal(as.character(built$layout$layout$period), c("1", "2"))
  expect_equal(nrow(built$data[[4]]), 4)
})

test_that("new arguments and resolution are validated", {
  dat <- rarity_records("survey", 2000, 1:5, 2:6)
  expect_error(rarity_call(dat, res = 0), "res must")
  expect_error(rarity_call(dat, prevPerPeriod = NA), "prevPerPeriod")
  expect_error(rarity_call(dat, returnDiagnostics = 1), "returnDiagnostics")
  expect_error(rarity_call(dat, plotRegressions = NA), "plotRegressions")
  expect_error(rarity_call(dat, labelTopResiduals = 1.5), "non-negative integer")
  dat$x <- as.character(dat$x)
  expect_error(rarity_call(dat), "coordinate columns must be numeric")
})

test_that("default identifier type is preserved", {
  dat <- rarity_records("survey", 2000, 1:5, 2:6)
  dat$identifier <- factor(dat$identifier)
  result <- rarity_call(dat)
  expect_s3_class(result$data$id, "factor")
  expect_equal(as.character(result$data$id), "survey")
})
