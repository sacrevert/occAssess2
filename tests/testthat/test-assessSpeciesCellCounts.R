make_species_cell_data <- function(species, x, y, year,
                                   identifier = "source-a",
                                   uncertainty = 0) {
  data.frame(
    species = species,
    x = x,
    y = y,
    year = year,
    spatialUncertainty = rep(uncertainty, length.out = length(year)),
    identifier = rep(identifier, length.out = length(year)),
    stringsAsFactors = FALSE
  )
}

run_species_cell_counts <- function(dat, periods, res, ...) {
  assessSpeciesCellCounts(
    dat = dat,
    species = "species",
    x = "x",
    y = "y",
    year = "year",
    spatialUncertainty = "spatialUncertainty",
    identifier = "identifier",
    periods = periods,
    res = res,
    ...
  )
}

test_that("independent five-record example gives 5, 4, and 3", {
  dat <- make_species_cell_data(
    species = c("A", "A", "B", "A", NA),
    x = c(0.1, 0.1, 0.1, 1.1, 2.1),
    y = rep(0.1, 5),
    year = rep(2000, 5)
  )

  result <- run_species_cell_counts(dat, list(2000), res = 1)

  expect_equal(result$data$nRecords, 5)
  expect_equal(result$data$nEligibleRecords, 4)
  expect_equal(result$data$nSpeciesCells, 3)
  expect_true(
    result$data$nSpeciesCells <= result$data$nEligibleRecords &&
      result$data$nEligibleRecords <= result$data$nRecords
  )
})

test_that("species and cells contribute independently while repetitions do not", {
  dat <- make_species_cell_data(
    species = c("A", "A", "B", "A", "A", "A"),
    x = c(0.1, 0.1, 0.1, 1.1, 0.1, 0.1),
    y = rep(0.1, 6),
    year = c(2000, 2000, 2000, 2000, 2001, 2002)
  )

  result <- run_species_cell_counts(dat, list(2000:2001, 2002), res = 1)

  expect_equal(result$data$nRecords, c(5, 1))
  expect_equal(result$data$nEligibleRecords, c(5, 1))
  expect_equal(result$data$nSpeciesCells, c(3, 1))
})

test_that("identifier-period combinations are separate and empty periods are zero", {
  dat <- make_species_cell_data(
    species = rep("A", 4),
    x = c(0.1, 0.1, 0.1, 1.1),
    y = rep(0.1, 4),
    year = c(2000, 2000, 2000, 2001),
    identifier = c("first", "first", "second", "second")
  )

  result <- run_species_cell_counts(dat, list(2000, 2001, 2002), res = 1)

  expect_equal(
    result$data,
    data.frame(
      identifier = rep(c("first", "second"), each = 3),
      Period = rep(1:3, 2),
      nRecords = c(2L, 0L, 0L, 1L, 1L, 0L),
      nEligibleRecords = c(2L, 0L, 0L, 1L, 1L, 0L),
      nSpeciesCells = c(1L, 0L, 0L, 1L, 1L, 0L)
    )
  )
})

test_that("missing species and coordinates affect eligibility but not raw counts", {
  dat <- make_species_cell_data(
    species = c("A", " ", NA, "B", "C", "D"),
    x = c(0.1, 0.1, 0.1, NA, Inf, 1.1),
    y = rep(0.1, 6),
    year = rep(2000, 6)
  )

  result <- run_species_cell_counts(dat, list(2000), res = 1)

  expect_equal(result$data$nRecords, 6)
  expect_equal(result$data$nEligibleRecords, 2)
  expect_equal(result$data$nSpeciesCells, 2)
})

test_that("missing years are warned about and outside-period rows are excluded", {
  dat <- make_species_cell_data(
    species = rep("A", 3),
    x = c(0.1, 0.1, 1000000),
    y = c(0.1, 0.1, 1000000),
    year = c(2000, NA, 1999)
  )

  expect_warning(
    result <- run_species_cell_counts(dat, list(2000, 2001), res = 1),
    "Removing 1 record with a missing or non-finite year",
    fixed = TRUE
  )
  expect_equal(result$data$nRecords, c(1, 0))
  expect_equal(result$data$nSpeciesCells, c(1, 0))
})

test_that("spatial uncertainty filtering applies to every metric", {
  dat <- make_species_cell_data(
    species = c("A", "B", "C", "D"),
    x = c(0.1, 0.1, 1.1, 2.1),
    y = rep(0.1, 4),
    year = rep(2000, 4),
    uncertainty = c(0, 10, 11, NA)
  )

  unfiltered <- run_species_cell_counts(dat, list(2000), res = 1)
  filtered <- run_species_cell_counts(
    dat, list(2000), res = 1, maxSpatUncertainty = 10
  )

  expect_equal(unfiltered$data[1, 3:5], data.frame(
    nRecords = 4L, nEligibleRecords = 4L, nSpeciesCells = 4L
  ))
  expect_equal(filtered$data[1, 3:5], data.frame(
    nRecords = 2L, nEligibleRecords = 2L, nSpeciesCells = 2L
  ))
})

test_that("geographic and projected coordinates use one zero-aligned convention", {
  geographic <- make_species_cell_data(
    species = rep("A", 3),
    x = c(-3.8, -3.2, -2.8),
    y = rep(55.2, 3),
    year = rep(2000, 3)
  )
  projected <- make_species_cell_data(
    species = rep("A", 3),
    x = c(430100, 430900, 431100),
    y = rep(350100, 3),
    year = rep(2000, 3)
  )

  expect_equal(
    run_species_cell_counts(geographic, list(2000), 1)$data$nSpeciesCells,
    2
  )
  expect_equal(
    run_species_cell_counts(projected, list(2000), 1000)$data$nSpeciesCells,
    2
  )
})

test_that("fractional resolutions and exact boundaries are deterministic", {
  dat <- make_species_cell_data(
    species = rep("A", 5),
    x = c(-0.5, -0.001, 0, 0.499, 0.5),
    y = rep(0.1, 5),
    year = rep(2000, 5)
  )

  result <- run_species_cell_counts(dat, list(2000), res = 0.5)

  expect_equal(result$data$nEligibleRecords, 5)
  expect_equal(result$data$nSpeciesCells, 3)
})

test_that("period, resolution, identifier, and empty-data errors are explicit", {
  dat <- make_species_cell_data("A", 0, 0, 2000)

  for (badResolution in list(0, -1, Inf, NA_real_, c(1, 2))) {
    expect_error(
      run_species_cell_counts(dat, list(2000), badResolution),
      "res must be a single positive, finite number",
      fixed = TRUE
    )
  }
  expect_error(
    run_species_cell_counts(dat, list(numeric()), 1),
    "Each period must be a non-empty finite numeric vector",
    fixed = TRUE
  )
  expect_error(
    run_species_cell_counts(dat, list(2000:2001, 2001:2002), 1),
    "Periods must not overlap",
    fixed = TRUE
  )
  expect_error(
    run_species_cell_counts(dat[0, ], list(2000), 1),
    "dat contains no records",
    fixed = TRUE
  )
  expect_error(
    run_species_cell_counts(dat, list(1999), 1),
    "No records fall within the requested periods",
    fixed = TRUE
  )

  missingIdentifier <- dat
  missingIdentifier$identifier <- NA_character_
  expect_error(
    run_species_cell_counts(missingIdentifier, list(2000), 1),
    "Missing values are not permitted in identifier",
    fixed = TRUE
  )

  noEligible <- dat
  noEligible$species <- NA_character_
  expect_error(
    run_species_cell_counts(noEligible, list(2000), 1),
    "No records with usable species identities and finite coordinates",
    fixed = TRUE
  )
})

test_that("raw counts agree with assessRecordNumber for equivalent valid input", {
  dat <- make_species_cell_data(
    species = c("A", "A", "B", "A", "B"),
    x = c(0.1, 0.1, 1.1, 0.1, 1.1),
    y = rep(0.1, 5),
    year = c(2000, 2001, 2001, 2002, 2002),
    identifier = c("first", "first", "first", "second", "second")
  )
  periods <- list(2000:2001, 2002, 2003)

  speciesCells <- run_species_cell_counts(dat, periods, res = 1)$data
  records <- assessRecordNumber(
    dat = dat, species = "species", x = "x", y = "y", year = "year",
    spatialUncertainty = "spatialUncertainty", identifier = "identifier",
    periods = periods
  )$data

  expect_equal(speciesCells$nRecords, records$val)
  expect_equal(speciesCells$identifier, records$group)
  expect_equal(speciesCells$Period, records$Period)
})

test_that("return data and plot expose both principal metrics", {
  dat <- make_species_cell_data(
    species = c("A", "A", "B"),
    x = c(0.1, 0.1, 1.1),
    y = rep(0.1, 3),
    year = c(2000, 2000, 2001)
  )
  result <- run_species_cell_counts(dat, list(2000, 2001), res = 1)

  expect_named(result, c("data", "plot"))
  expect_named(
    result$data,
    c("identifier", "Period", "nRecords", "nEligibleRecords", "nSpeciesCells")
  )
  expect_s3_class(result$plot, "ggplot")
  expect_equal(
    levels(result$plot$data$Measure),
    c("Raw records", "Species-cell combinations")
  )
  expect_equal(result$plot$data$Count, c(2L, 1L, 1L, 1L))
  expect_no_error(ggplot2::ggplot_build(result$plot))
})
