make_spatial_cov_data <- function(x, y, year, identifier = "source-a",
                                  uncertainty = 0) {
  data.frame(
    species = rep("species-a", length(x)),
    x = x,
    y = y,
    year = year,
    spatialUncertainty = rep(uncertainty, length.out = length(x)),
    identifier = rep(identifier, length.out = length(x))
  )
}

run_spatial_cov <- function(dat, periods, res, output = "density", ...) {
  assessSpatialCov(
    dat = dat,
    species = "species",
    x = "x",
    y = "y",
    year = "year",
    spatialUncertainty = "spatialUncertainty",
    identifier = "identifier",
    periods = periods,
    res = res,
    output = output,
    countries = NULL,
    returnRaster = TRUE,
    ...
  )
}

raster_values <- function(result, layer = 1) {
  unname(raster::getValues(result$rasters[[layer]]))
}

make_mode_test_data <- function() {
  p1 <- make_spatial_cov_data(
    x = c(0.2, rep(0.2, 3), rep(1.2, 2)),
    y = c(1.2, rep(0.2, 3), rep(0.2, 2)),
    year = rep(2000, 6)
  )
  p2 <- make_spatial_cov_data(
    x = c(rep(1.2, 4), 0.2),
    y = c(rep(1.2, 4), 0.2),
    year = rep(2001, 5)
  )
  rbind(p1, p2)
}

test_that("density rasters contain explicit occurrence-record counts", {
  dat <- make_spatial_cov_data(
    x = c(0.2, rep(0.2, 3), rep(1.2, 2)),
    y = c(1.2, rep(0.2, 3), rep(0.2, 2)),
    year = rep(2000, 6)
  )

  result <- run_spatial_cov(dat, list(2000), res = 1)

  # Raster order is upper-left, upper-right, lower-left, lower-right.
  expect_equal(raster_values(result), c(1, NA, 3, 2))
  expect_equal(raster::res(result$rasters), c(1, 1))
  expect_equal(c(raster::nrow(result$rasters), raster::ncol(result$rasters)), c(2, 2))
  expect_equal(names(result$rasters), "id.source.a.period.p1")
})

test_that("requested periods drive density layer order including empty periods", {
  dat <- make_mode_test_data()
  result <- run_spatial_cov(dat, list(2000, 2001, 2002), res = 1)

  expect_equal(raster_values(result, 1), c(1, NA, 3, 2))
  expect_equal(raster_values(result, 2), c(NA, 4, 1, NA))
  expect_true(all(is.na(raster_values(result, 3))))
  expect_equal(
    names(result$rasters),
    paste0("id.source.a.period.p", 1:3)
  )
  expect_equal(levels(result$plots[[1]]$data$variable), paste0("p", 1:3))
})

test_that("nPeriods counts raw period presence independently of logCount", {
  dat <- make_mode_test_data()
  periods <- list(2000, 2001, 2002)

  plain <- run_spatial_cov(dat, periods, res = 1, output = "nPeriods")
  logged <- run_spatial_cov(
    dat, periods, res = 1, output = "nPeriods", logCount = TRUE
  )

  expect_equal(raster_values(plain), c(1, 1, 2, 1))
  expect_equal(raster_values(logged), raster_values(plain))
  expect_equal(names(plain$rasters), "id.source.a")
  expect_equal(plain$plots[[1]]$data$value, raster_values(plain))
  expect_equal(
    plain$plots[[1]]$scales$get_scales("fill")$name,
    "Number of periods sampled"
  )
})

test_that("overlap applies valid thresholds to raw sampled-period counts", {
  dat <- make_mode_test_data()
  periods <- list(2000, 2001, 2002)

  at_least_one <- run_spatial_cov(
    dat, periods, 1, output = "overlap", minPeriods = 1
  )
  at_least_two <- run_spatial_cov(
    dat, periods, 1, output = "overlap", minPeriods = 2
  )
  at_least_two_logged <- run_spatial_cov(
    dat, periods, 1, output = "overlap", minPeriods = 2,
    logCount = TRUE
  )
  all_requested <- run_spatial_cov(dat, periods, 1, output = "overlap")

  expect_equal(raster_values(at_least_one), rep(TRUE, 4))
  expect_equal(raster_values(at_least_two), c(NA, NA, TRUE, NA))
  expect_equal(raster_values(at_least_two_logged), raster_values(at_least_two))
  expect_true(all(is.na(raster_values(all_requested))))
  expect_no_error(ggplot2::ggplot_build(all_requested$plots[[1]]))
})

test_that("logCount transforms density values only after raw counting", {
  dat <- make_mode_test_data()
  result <- run_spatial_cov(
    dat, list(2000, 2001), res = 1, logCount = TRUE
  )

  expect_equal(raster_values(result, 1), c(0, NA, log10(3), log10(2)))
  expect_equal(raster_values(result, 2), c(NA, log10(4), 0, NA))
  expect_equal(
    result$plots[[1]]$scales$get_scales("fill")$name,
    "log10(record count)"
  )
  expect_equal(
    result$plots[[1]]$data$value,
    unlist(lapply(1:2, function(i) raster_values(result, i)))
  )
})

test_that("multiple identifiers use one grid with ordered unambiguous layers", {
  dat <- make_spatial_cov_data(
    x = c(0.2, 0.2, 0.2, 1.2, 1.2, 1.2),
    y = rep(0.2, 6),
    year = c(2000, 2000, 2000, 2000, 2001, 2001),
    identifier = c(rep("source-a", 3), rep("source-b", 3))
  )

  density <- run_spatial_cov(dat, list(2000, 2001), res = 1)
  n_periods <- run_spatial_cov(
    dat, list(2000, 2001), res = 1, output = "nPeriods"
  )

  expect_equal(
    names(density$rasters),
    c(
      "id.source.a.period.p1", "id.source.a.period.p2",
      "id.source.b.period.p1", "id.source.b.period.p2"
    )
  )
  expect_equal(raster_values(density, 1), c(3, NA))
  expect_true(all(is.na(raster_values(density, 2))))
  expect_equal(raster_values(density, 3), c(NA, 1))
  expect_equal(raster_values(density, 4), c(NA, 2))
  expect_equal(names(density$plots), c("source-a", "source-b"))

  expect_equal(names(n_periods$rasters), c("id.source.a", "id.source.b"))
  expect_equal(raster_values(n_periods, 1), c(1, NA))
  expect_equal(raster_values(n_periods, 2), c(NA, 2))
})

test_that("geographic and projected coordinates use the requested resolution", {
  geographic <- make_spatial_cov_data(
    x = c(-3.8, -3.2, -2.8),
    y = rep(55.2, 3),
    year = rep(2000, 3)
  )
  projected <- make_spatial_cov_data(
    x = c(430100, 430900, 431100),
    y = rep(350100, 3),
    year = rep(2000, 3)
  )

  geographic_result <- run_spatial_cov(geographic, list(2000), res = 1)
  projected_result <- run_spatial_cov(projected, list(2000), res = 1000)

  expect_equal(raster_values(geographic_result), c(2, 1))
  expect_equal(raster::res(geographic_result$rasters), c(1, 1))
  expect_equal(raster_values(projected_result), c(2, 1))
  expect_equal(raster::res(projected_result$rasters), c(1000, 1000))
})

test_that("non-integer resolutions and exact boundaries are deterministic", {
  positive <- make_spatial_cov_data(
    x = c(0, 0.299, 0.3, 0.6),
    y = rep(0.1, 4),
    year = rep(2000, 4)
  )
  negative <- make_spatial_cov_data(
    x = c(-0.3, -0.001, 0),
    y = rep(-0.1, 3),
    year = rep(2000, 3)
  )

  positive_result <- run_spatial_cov(positive, list(2000), res = 0.3)
  negative_result <- run_spatial_cov(negative, list(2000), res = 0.3)

  expect_equal(raster_values(positive_result), c(2, 1, 1))
  expect_equal(raster::res(positive_result$rasters), c(0.3, 0.3))
  expect_equal(raster::ncol(positive_result$rasters), 3)
  expect_equal(raster_values(negative_result), c(2, 1))
  expect_equal(raster::ncol(negative_result$rasters), 2)
})

test_that("missing data are explicit and outside-period rows do not affect the grid", {
  missing_year <- make_spatial_cov_data(
    x = c(0.2, 0.2), y = c(0.2, 0.2), year = c(2000, NA)
  )
  missing_coordinates <- make_spatial_cov_data(
    x = c(0.2, NA, 0.2), y = c(0.2, 0.2, Inf), year = rep(2000, 3)
  )
  outside <- make_spatial_cov_data(
    x = c(0.2, 1e9), y = c(0.2, 1e9), year = c(2000, 1999)
  )

  expect_warning(
    year_result <- run_spatial_cov(missing_year, list(2000), 1),
    "Removing 1 record with a missing or non-finite year",
    fixed = TRUE
  )
  expect_warning(
    coordinate_result <- run_spatial_cov(missing_coordinates, list(2000), 1),
    "Removing 2 records with missing or non-finite coordinates",
    fixed = TRUE
  )
  outside_result <- run_spatial_cov(outside, list(2000), 1)

  expect_equal(raster_values(year_result), 1)
  expect_equal(raster_values(coordinate_result), 1)
  expect_equal(raster_values(outside_result), 1)
  expect_equal(raster::ncell(outside_result$rasters), 1)
})

test_that("single locations and equal counts construct valid plots", {
  dat <- make_spatial_cov_data(x = 12, y = 48, year = 2000)
  result <- run_spatial_cov(dat, list(2000), res = 10)

  expect_equal(raster_values(result), 1)
  expect_equal(raster::res(result$rasters), c(10, 10))
  expect_s3_class(result$plots[[1]], "ggplot")
  expect_no_error(ggplot2::ggplot_build(result$plots[[1]]))
  expect_equal(result$plots[[1]]$data$value, 1)
})

test_that("period labels remain in numeric request order beyond p9", {
  periods <- lapply(2001:2011, function(value) value)
  dat <- make_spatial_cov_data(
    x = rep(0.2, 11), y = rep(0.2, 11), year = 2001:2011
  )
  result <- run_spatial_cov(dat, periods, res = 1)

  expected_periods <- paste0("p", 1:11)
  expect_equal(
    names(result$rasters),
    paste0("id.source.a.period.", expected_periods)
  )
  expect_equal(levels(result$plots[[1]]$data$variable), expected_periods)
  expect_equal(vapply(1:11, function(i) raster_values(result, i), numeric(1)), rep(1, 11))
})

test_that("invalid resolutions, outputs, periods, and thresholds are rejected", {
  dat <- make_spatial_cov_data(x = 0, y = 0, year = 2000)

  for (bad_resolution in list(0, -1, Inf, NA_real_, c(1, 2))) {
    expect_error(
      run_spatial_cov(dat, list(2000), bad_resolution),
      "res must be a single positive, finite number",
      fixed = TRUE
    )
  }
  expect_error(
    run_spatial_cov(dat, list(2000), 1, output = "unsupported"),
    "output must be one of",
    fixed = TRUE
  )
  expect_error(
    run_spatial_cov(dat, list(2000:2001, 2001:2002), 1),
    "Periods must not overlap",
    fixed = TRUE
  )
  expect_error(
    run_spatial_cov(dat, list(numeric()), 1),
    "Each period must be a non-empty finite numeric vector",
    fixed = TRUE
  )

  for (bad_threshold in list(0, 3, 1.5, Inf, NA_real_, c(1, 2))) {
    expect_error(
      run_spatial_cov(
        dat, list(2000, 2001), 1, output = "overlap",
        minPeriods = bad_threshold
      ),
      "minPeriods must be a single integer between 1 and length(periods)",
      fixed = TRUE
    )
  }
})

test_that("empty selections, missing identifiers, and huge grids fail clearly", {
  dat <- make_spatial_cov_data(x = 0, y = 0, year = 2000)
  missing_identifier <- dat
  missing_identifier$identifier <- NA_character_
  huge <- make_spatial_cov_data(
    x = c(0, 1e6), y = c(0, 1e6), year = c(2000, 2000)
  )

  expect_error(
    run_spatial_cov(dat[0, ], list(2000), 1),
    "dat contains no records",
    fixed = TRUE
  )
  expect_error(
    run_spatial_cov(dat, list(1999), 1),
    "No records fall within the requested periods",
    fixed = TRUE
  )
  expect_error(
    run_spatial_cov(missing_identifier, list(2000), 1),
    "Missing values are not permitted in identifier",
    fixed = TRUE
  )
  expect_error(
    run_spatial_cov(huge, list(2000), 1),
    "exceed ten million cell values",
    fixed = TRUE
  )
})

test_that("return interfaces and plotted density values are unambiguous", {
  dat <- make_mode_test_data()
  with_rasters <- run_spatial_cov(dat, list(2000, 2001), res = 1)
  plots_only <- assessSpatialCov(
    dat = dat,
    species = "species",
    x = "x",
    y = "y",
    year = "year",
    spatialUncertainty = "spatialUncertainty",
    identifier = "identifier",
    periods = list(2000, 2001),
    res = 1,
    output = "density",
    countries = NULL,
    returnRaster = FALSE
  )

  expect_named(with_rasters, c("plots", "rasters"))
  expect_s4_class(with_rasters$rasters, "RasterStack")
  expect_named(plots_only, "source-a")
  expect_s3_class(plots_only[[1]], "ggplot")
  expect_equal(
    with_rasters$plots[[1]]$data$value,
    c(raster_values(with_rasters, 1), raster_values(with_rasters, 2))
  )
  expect_equal(
    with_rasters$plots[[1]]$scales$get_scales("fill")$name,
    "Record count"
  )
})
