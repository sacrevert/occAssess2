make_repeat_visit_data <- function(lon, lat, year, identifier = "source-a") {
  data.frame(
    lon = lon,
    lat = lat,
    year = year,
    identifier = identifier
  )
}

occupied_repeat_values <- function(result, identifier, period) {
  values <- result$data$vals[
    result$data$identifier == identifier & result$data$Period == period
  ]
  sort(values[!is.na(values)])
}

test_that("occurrence rows are counted within requested-resolution cells", {
  dat <- make_repeat_visit_data(
    lon = c(0.1, 0.8, 1.2),
    lat = c(0.2, 0.7, 0.5),
    year = rep(2000, 3)
  )

  result <- assessRepeatVisits(dat, periods = list(2000), spatRes = 1)

  expect_equal(result$data$vals, c(2, 1))
  expect_equal(result$data$identifier, rep("source-a", 2))
  expect_equal(result$data$Period, rep("p1", 2))
  expect_named(result, c("data", "plot"))
  expect_named(result$data, c("vals", "identifier", "Period"))
})

test_that("period totals are divided by all years specified in each period", {
  dat <- make_repeat_visit_data(
    lon = rep(0.2, 9),
    lat = rep(0.2, 9),
    year = c(2000, 2000, 2001, rep(2002, 6))
  )

  result <- assessRepeatVisits(
    dat,
    periods = list(2000:2001, 2002:2004),
    spatRes = 1
  )

  expect_equal(result$data$vals, c(3 / 2, 6 / 3))
  expect_equal(result$data$Period, c("p1", "p2"))
})

test_that("identifiers are calculated separately on one common grid", {
  dat <- make_repeat_visit_data(
    lon = c(0.1, 0.1, 0.1, 0.1, 1.1, 1.1),
    lat = rep(0.1, 6),
    year = rep(2000, 6),
    identifier = c(rep("source-a", 3), rep("source-b", 3))
  )

  result <- assessRepeatVisits(dat, list(2000), spatRes = 1)

  expect_equal(occupied_repeat_values(result, "source-a", "p1"), 3)
  expect_equal(occupied_repeat_values(result, "source-b", "p1"), c(1, 2))
  expect_equal(nrow(result$data), 4)
})

test_that("grids support geographic and British National Grid coordinates", {
  geographic <- make_repeat_visit_data(
    lon = c(-3.8, -3.2, -2.8),
    lat = rep(55.2, 3),
    year = rep(2000, 3)
  )
  bng <- make_repeat_visit_data(
    lon = c(430100, 430900, 431100),
    lat = rep(350100, 3),
    year = rep(2000, 3)
  )

  geographic_result <- assessRepeatVisits(geographic, list(2000), spatRes = 1)
  bng_result <- assessRepeatVisits(bng, list(2000), spatRes = 1000)

  expect_equal(occupied_repeat_values(geographic_result, "source-a", "p1"), c(1, 2))
  expect_equal(occupied_repeat_values(bng_result, "source-a", "p1"), c(1, 2))
})

test_that("cell boundaries are aligned to zero at the exact resolution", {
  dat <- make_repeat_visit_data(
    lon = c(0.2, 0.2, 1, 0.2, 0.2, 0.2),
    lat = c(0.2, 0.2, 0.2, 1, 1, 1),
    year = rep(2000, 6)
  )

  result <- assessRepeatVisits(dat, list(2000), spatRes = 1)

  # Independent reference for four cells [0,1) x [0,1), [1,2) x [0,1),
  # [0,1) x [1,2), and [1,2) x [1,2).
  expect_equal(result$data$vals, c(3, NA, 2, 1))
  expect_equal(nrow(result$data), 4)

  decimal_boundaries <- make_repeat_visit_data(
    lon = c(0.1, 0.2, 0.3),
    lat = rep(0.05, 3),
    year = rep(2000, 3)
  )
  decimal_result <- assessRepeatVisits(
    decimal_boundaries,
    list(2000),
    spatRes = 0.1
  )

  expect_equal(decimal_result$data$vals, c(1, 1, 1))
})

test_that("missing observations are explicit and out-of-period rows are excluded", {
  missing_year <- make_repeat_visit_data(
    lon = c(0.2, 0.2), lat = c(0.2, 0.2), year = c(2000, NA)
  )
  missing_coordinates <- make_repeat_visit_data(
    lon = c(0.2, NA, 0.2), lat = c(0.2, 0.2, Inf), year = rep(2000, 3)
  )
  outside <- make_repeat_visit_data(
    lon = c(0.2, 1e9), lat = c(0.2, 1e9), year = c(2000, 1999)
  )

  expect_warning(
    year_result <- assessRepeatVisits(missing_year, list(2000), 1),
    "Removing 1 record with a missing or non-finite year",
    fixed = TRUE
  )
  expect_warning(
    coordinate_result <- assessRepeatVisits(missing_coordinates, list(2000), 1),
    "Removing 2 records with missing or non-finite coordinates",
    fixed = TRUE
  )
  outside_result <- assessRepeatVisits(outside, list(2000), 1)

  expect_equal(year_result$data$vals, 1)
  expect_equal(coordinate_result$data$vals, 1)
  expect_equal(outside_result$data$vals, 1)
  expect_equal(nrow(outside_result$data), 1)
})

test_that("empty periods and single-location datasets have defined output", {
  dat <- make_repeat_visit_data(lon = 12, lat = 48, year = 2000)

  result <- assessRepeatVisits(dat, list(2000, 2001:2002), spatRes = 10)

  expect_equal(result$data$vals, c(1, NA))
  expect_equal(result$data$Period, c("p1", "p2"))
  expect_equal(nrow(result$plot$data), 1)
  expect_equal(levels(result$plot$data$Period), c("p1", "p2"))
})

test_that("invalid data, periods, and resolutions fail informatively", {
  dat <- make_repeat_visit_data(lon = 0, lat = 0, year = 2000)
  missing_identifier <- dat
  missing_identifier$identifier <- NA_character_

  expect_error(
    assessRepeatVisits(dat[0, ], list(2000), 1),
    "dat contains no records",
    fixed = TRUE
  )
  expect_error(
    assessRepeatVisits(missing_identifier, list(2000), 1),
    "Missing values are not permitted in identifier",
    fixed = TRUE
  )
  expect_error(
    assessRepeatVisits(dat, list(), 1),
    "periods must be a non-empty list",
    fixed = TRUE
  )
  expect_error(
    assessRepeatVisits(dat, list(numeric()), 1),
    "Each period must be a non-empty finite numeric vector",
    fixed = TRUE
  )
  expect_error(
    assessRepeatVisits(dat, list(c(2000, 2000)), 1),
    "without duplicate years",
    fixed = TRUE
  )
  expect_error(
    assessRepeatVisits(dat, list(2000:2001, 2001:2002), 1),
    "Periods must not overlap",
    fixed = TRUE
  )
  expect_error(
    assessRepeatVisits(dat, list(1990), 1),
    "No records fall within the requested periods",
    fixed = TRUE
  )

  for (bad_resolution in list(0, -1, Inf, NA_real_, c(1, 2))) {
    expect_error(
      assessRepeatVisits(dat, list(2000), bad_resolution),
      "spatRes must be a single positive, finite number",
      fixed = TRUE
    )
  }
})

test_that("impractically large grids are rejected before allocation", {
  dat <- make_repeat_visit_data(
    lon = c(0, 1e6),
    lat = c(0, 1e6),
    year = c(2000, 2000)
  )

  expect_error(
    assessRepeatVisits(dat, list(2000), spatRes = 1),
    "exceed ten million output rows",
    fixed = TRUE
  )
})

test_that("values above twenty remain in both data and plot", {
  dat <- make_repeat_visit_data(
    lon = rep(0.2, 25),
    lat = rep(0.2, 25),
    year = rep(2000, 25)
  )

  result <- assessRepeatVisits(dat, list(2000), spatRes = 1)
  built <- ggplot2::ggplot_build(result$plot)

  expect_equal(result$data$vals, 25)
  expect_s3_class(result$plot, "ggplot")
  expect_equal(result$plot$labels$x, "Occurrence records per grid cell per year in period")
  expect_gte(max(built$layout$panel_params[[1]]$x.range), 25)
})
