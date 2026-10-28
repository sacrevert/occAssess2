make_temporal_occurrences <- function() {
  data.frame(
    species = c("A", "A", "B", "C", "D", "E"),
    x = c(0.1, 0.1, 1.1, 2.1, 3.1, 4.1),
    y = rep(0.1, 6),
    year = c(1995, 2005, 2005, 2005, 2005, 2005),
    start = c(NA, "1999", "2005", "2005", "2005", "2005"),
    end = c(NA, "2003", "2005", "2005", "2005", "2005"),
    spatialUncertainty = rep(0, 6),
    identifier = rep("survey", 6),
    stringsAsFactors = FALSE
  )
}

temporal_args <- list(
  start = "start", end = "end", ambiguous = "maxOverlap"
)

test_that("record, species, and identification summaries use one assignment", {
  dat <- make_temporal_occurrences()
  periods <- list(1990:1999, 2000:2009)

  records <- suppressWarnings(assessRecordNumber(
    dat, "species", "x", "y", "year", "spatialUncertainty", "identifier",
    periods, dateHandling = temporal_args
  ))
  species <- suppressWarnings(assessSpeciesNumber(
    dat, "species", "x", "y", "year", "spatialUncertainty", "identifier",
    periods, dateHandling = temporal_args
  ))
  identified <- suppressWarnings(assessSpeciesID(
    dat, "species", "x", "y", "year", "spatialUncertainty", "identifier",
    periods, type = "count", dateHandling = temporal_args
  ))

  expect_equal(records$data$val, c(1, 5))
  expect_equal(species$data$val, c(1, 5))
  expect_equal(identified$data$prop, c(1, 5))
  expect_equal(attr(records, "temporalDiagnostics")$periodCounts, c(p1 = 1L, p2 = 5L))
})

test_that("uncertain dates are summarised separately for multiple identifiers", {
  dat <- make_temporal_occurrences()[c(1, 2, 2), ]
  dat$identifier <- c("alpha", "alpha", "beta")
  result <- suppressWarnings(assessRecordNumber(
    dat, "species", "x", "y", "year", "spatialUncertainty", "identifier",
    periods = list(1990:1999, 2000:2009), dateHandling = temporal_args
  ))

  expected <- data.frame(
    val = c(1, 1, 0, 1),
    group = c("alpha", "alpha", "beta", "beta"),
    Period = c(1, 2, 1, 2)
  )
  expect_equal(
    result$data[, c("val", "group", "Period")],
    expected,
    ignore_attr = TRUE
  )
  expect_equal(sum(result$data$val), 3)
})

test_that("species-cell counts do not duplicate uncertain observations", {
  dat <- make_temporal_occurrences()
  result <- suppressWarnings(assessSpeciesCellCounts(
    dat, "species", "x", "y", "year", "spatialUncertainty", "identifier",
    periods = list(1990:1999, 2000:2009), res = 1,
    dateHandling = temporal_args
  ))

  expect_equal(result$data$nRecords, c(1L, 5L))
  expect_equal(sum(result$data$nRecords), nrow(dat))
  expect_equal(result$data$nSpeciesCells, c(1L, 5L))
})

test_that("spatial coverage counts an uncertain record in one period only", {
  dat <- make_temporal_occurrences()[1:2, ]
  result <- suppressWarnings(assessSpatialCov(
    dat, "species", "x", "y", "year", "spatialUncertainty", "identifier",
    res = 1, periods = list(1990:1999, 2000:2009), countries = NULL,
    returnRaster = TRUE, dateHandling = temporal_args
  ))

  counts <- vapply(
    seq_len(raster::nlayers(result$rasters)),
    function(layer) sum(raster::getValues(result$rasters[[layer]]), na.rm = TRUE),
    numeric(1)
  )
  expect_equal(counts, c(1, 1))
  expect_equal(sum(counts), 2)
})

test_that("repeat-visit denominator remains the supplied period length", {
  dat <- data.frame(
    lon = c(0.1, 0.1), lat = c(0.1, 0.1),
    start = c("2001", "2004"), end = c("2003", "2004"),
    identifier = c("survey", "survey")
  )
  result <- assessRepeatVisits(
    dat, periods = list(2000:2004), spatRes = 1,
    dateHandling = list(start = "start", end = "end")
  )

  expect_equal(stats::na.omit(result$data$vals), 2 / 5)
})

test_that("spatial-uncertainty proportion uses assigned records only", {
  dat <- data.frame(
    species = "A", x = 1:3, y = 1:3,
    spatialUncertainty = c(10, NA, 30), identifier = "survey",
    start = c("2001", "1999", NA), end = c("2002", "2001", NA)
  )
  result <- suppressWarnings(assessSpatialUncertainty(
    dat, periods = list(2000:2009),
    dateHandling = list(start = "start", end = "end", ambiguous = "drop")
  ))

  expect_equal(unname(result$data[1, 1]), 1)
  expect_equal(attr(result, "temporalDiagnostics")$periodCounts, c(p1 = 1L))
})

test_that("environmental rows remain aligned after uncertain-date filtering", {
  dat <- data.frame(
    species = rep("A", 5), x = 1:5, y = 1:5, year = NA_real_,
    spatialUncertainty = 0, identifier = "survey",
    start = c("2001", "1999", "2003", "2004", "1980"),
    end = c("2001", "2001", "2003", "2004", "1981")
  )
  env <- data.frame(v1 = c(1, 1000, 3, 4, 2000), v2 = c(2, 1000, 6, 8, 2000))
  background <- data.frame(v1 = c(0, 2, 4, 6), v2 = c(0, 4, 8, 12))
  result <- suppressWarnings(assessEnvBias(
    dat, "species", "x", "y", "spatialUncertainty", "year", "identifier",
    periods = list(2000:2009), envDat = env, backgroundEnvDat = background,
    dateHandling = list(start = "start", end = "end", ambiguous = "drop")
  ))

  expected <- stats::predict(stats::prcomp(background), env[c(1, 3, 4), ])
  observed <- result$data[
    result$data$Period != "background", c("scores.PC1", "scores.PC2")
  ]
  expect_equal(unname(as.matrix(observed)), unname(expected))
})

test_that("rarity and spatial-bias results match equivalent exact dates", {
  species <- rep(LETTERS[1:5], each = 3)
  dat <- data.frame(
    species = species,
    x = rep(c(0.1, 1.1, 2.1), 5) + rep(0:4, each = 3) * 0.01,
    y = rep(c(0.1, 1.1, 2.1), 5), year = 2005,
    spatialUncertainty = 0, identifier = "survey",
    start = "2005", end = "2005"
  )
  periods <- list(2000:2009)
  exactRarity <- suppressWarnings(assessRarityBias(
    dat, "species", "x", "y", "year", "spatialUncertainty", "identifier",
    periods, res = 1, prevPerPeriod = TRUE
  ))
  intervalRarity <- suppressWarnings(assessRarityBias(
    dat, "species", "x", "y", NULL, "spatialUncertainty", "identifier",
    periods, res = 1, prevPerPeriod = TRUE,
    dateHandling = list(start = "start", end = "end")
  ))
  expect_equal(intervalRarity$data, exactRarity$data)

  mask <- raster::raster(ncol = 20, nrow = 20, xmn = 0, xmx = 5, ymn = 0, ymx = 5)
  raster::values(mask) <- 1
  set.seed(10)
  exactSpatial <- suppressWarnings(assessSpatialBias(
    dat, "species", "x", "y", "year", "spatialUncertainty", "identifier",
    periods, mask = mask, nSamps = 3, degrade = FALSE
  ))
  set.seed(10)
  intervalSpatial <- suppressWarnings(assessSpatialBias(
    dat, "species", "x", "y", NULL, "spatialUncertainty", "identifier",
    periods, mask = mask, nSamps = 3, degrade = FALSE,
    dateHandling = list(start = "start", end = "end")
  ))
  expect_equal(intervalSpatial$data, exactSpatial$data)
})

test_that("exact-year outputs remain unchanged across simple assessments", {
  dat <- make_temporal_occurrences()
  dat$start <- NULL
  dat$end <- NULL
  periods <- list(1990:1999, 2000:2009)

  records <- assessRecordNumber(
    dat, "species", "x", "y", "year", "spatialUncertainty", "identifier", periods
  )
  species <- assessSpeciesNumber(
    dat, "species", "x", "y", "year", "spatialUncertainty", "identifier", periods
  )
  cells <- assessSpeciesCellCounts(
    dat, "species", "x", "y", "year", "spatialUncertainty", "identifier",
    periods, res = 1
  )

  expect_equal(records$data$val, c(1, 5))
  expect_equal(species$data$val, c(1, 5))
  expect_equal(cells$data$nRecords, c(1L, 5L))
})
