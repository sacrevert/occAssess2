make_spatial_uncertainty_data <- function(year, uncertainty, identifier = "source-a") {
  n <- length(year)
  data.frame(
    species = rep("species-a", n),
    x = seq_len(n),
    y = seq_len(n) + 10,
    year = year,
    spatialUncertainty = uncertainty,
    identifier = identifier
  )
}

test_that("all specified uncertainties produce a proportion of one", {
  periods <- list(2000:2002)
  dat <- make_spatial_uncertainty_data(
    year = 2000:2002,
    uncertainty = c(0, 100, 20000)
  )

  result <- assessSpatialUncertainty(dat, periods)

  expected <- matrix(
    1,
    ncol = 1,
    dimnames = list("source-a", "proportion_specified")
  )
  attr(expected, "periods") <- periods

  expect_equal(result$data, expected)
  expect_equal(attr(result$data, "periods"), periods)
  expect_s3_class(result$plot, "ggplot")

  built <- ggplot2::ggplot_build(result$plot)
  expect_gte(max(built$data[[1]]$xmax), 20000)
  expect_equal(
    built$data[[2]]$label,
    "Proportion specified (pooled) = 1.00"
  )
  expect_false(result$plot$layers[[2]]$inherit.aes)
})

test_that("all missing uncertainties produce a proportion of zero", {
  periods <- list(2000:2002)
  dat <- make_spatial_uncertainty_data(
    year = 2000:2002,
    uncertainty = rep(NA_real_, 3)
  )

  result <- assessSpatialUncertainty(dat, periods)

  expect_equal(unname(result$data[, 1]), 0)
  expect_warning(
    built <- ggplot2::ggplot_build(result$plot),
    "Removed 3 rows containing non-finite"
  )
  expect_equal(nrow(built$data[[1]]), 0)
  expect_equal(
    built$data[[2]]$label,
    "Proportion specified (pooled) = 0.00"
  )
})

test_that("pooled proportions retain identifiers across different period patterns", {
  periods <- list(2000:2001, 2002:2003)
  dat <- make_spatial_uncertainty_data(
    year = c(2000, 2001, 2002, 2003, 2000, 2001, 2002, 2003),
    uncertainty = c(1, 2, NA, NA, NA, 5, 6, 7),
    identifier = rep(c("source-a", "source-b"), each = 4)
  )

  result <- assessSpatialUncertainty(dat, periods)

  expect_equal(
    unname(result$data[, "proportion_specified"]),
    c(1 / 2, 3 / 4)
  )
  expect_equal(rownames(result$data), c("source-a", "source-b"))
  expect_equal(attr(result$data, "periods"), periods)
  expect_equal(
    result$plot$data$Period,
    c("p1", "p1", "p1", "p1", "p2", "p2", "p2", "p2")
  )

  annotation_data <- result$plot$layers[[2]]$data
  expect_equal(annotation_data$identifier, c("source-a", "source-b"))
  expect_false("Period" %in% names(annotation_data))
})

test_that("missing years and out-of-period observations are excluded", {
  periods <- list(2000:2001)
  dat <- make_spatial_uncertainty_data(
    year = c(2000, NA, 1999, 2001),
    uncertainty = c(10, 20, 30, NA)
  )

  expect_warning(
    result <- assessSpatialUncertainty(dat, periods),
    "Removing 1 record because it does not have a year associated",
    fixed = TRUE
  )

  expect_equal(unname(result$data[, 1]), 1 / 2)
  expect_equal(result$plot$data$year, c(2000, 2001))
  expect_equal(result$plot$data$spatialUncertainty, c(10, NA))
  expect_equal(result$plot$data$Period, c("p1", "p1"))
})

test_that("empty selections fail with an informative error", {
  empty <- make_spatial_uncertainty_data(
    year = numeric(),
    uncertainty = numeric(),
    identifier = character()
  )
  outside <- make_spatial_uncertainty_data(
    year = c(1998, 1999),
    uncertainty = c(1, 2)
  )

  expect_error(
    assessSpatialUncertainty(empty, list(2000:2001)),
    "No records fall within the requested periods.",
    fixed = TRUE
  )
  expect_error(
    assessSpatialUncertainty(outside, list(2000:2001)),
    "No records fall within the requested periods.",
    fixed = TRUE
  )
})
