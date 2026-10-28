temporal_periods <- list(1990:1999, 2000:2009, 2010:2019)

assign_intervals <- function(data, ..., periods = temporal_periods,
                             year = "year") {
  occAssess:::assignTemporalPeriods(
    data,
    periods = periods,
    year = year,
    dateHandling = utils::modifyList(
      list(start = "start", end = "end"),
      list(...)
    )
  )
}

test_that("exact years retain supplied period indices", {
  dat <- data.frame(year = c(1990, 2005, 2019, 1989, NA))
  assigned <- occAssess:::assignTemporalPeriods(dat, temporal_periods, "year")

  expect_equal(assigned$period, c(1L, 2L, 3L, NA, NA))
  expect_equal(
    assigned$status,
    c(
      "exact_assigned", "exact_assigned", "exact_assigned",
      "outside_periods", "missing_temporal"
    )
  )
  expect_equal(assigned$diagnostics$periodCounts, c(p1 = 1L, p2 = 1L, p3 = 1L))
})

test_that("same-year and contained multi-year intervals assign automatically", {
  dat <- data.frame(
    year = NA_real_,
    start = c("2005-03-01", "2002", "2000-01-01"),
    end = c("2005-11-30", "2007", "2009-12-31")
  )

  for (rule in c("drop", "start", "end", "maxOverlap")) {
    assigned <- assign_intervals(dat, ambiguous = rule)
    expect_equal(assigned$period, rep(2L, 3))
    expect_equal(assigned$status, rep("interval_unambiguous", 3))
  }
})

test_that("ambiguous bounded intervals obey all allocation rules", {
  dat <- data.frame(year = NA_real_, start = "1999", end = "2003")

  expect_true(is.na(assign_intervals(dat, ambiguous = "drop")$period))
  expect_equal(assign_intervals(dat, ambiguous = "start")$period, 1L)
  expect_equal(assign_intervals(dat, ambiguous = "end")$period, 2L)
  expect_equal(assign_intervals(dat, ambiguous = "maxOverlap")$period, 2L)
})

test_that("maximum overlap uses actual years, unequal periods, and explicit ties", {
  periods <- list(c(1990, 1992, 1994), 1995:2004, 2005:2006)
  dat <- data.frame(
    year = NA_real_,
    start = c("1990", "2002", "2004"),
    end = c("1994", "2006", "2005")
  )

  dropped <- assign_intervals(
    dat, periods = periods, ambiguous = "maxOverlap", ties = "drop"
  )
  earliest <- assign_intervals(
    dat, periods = periods, ambiguous = "maxOverlap", ties = "earliest"
  )
  latest <- assign_intervals(
    dat, periods = periods, ambiguous = "maxOverlap", ties = "latest"
  )

  # The first interval spans gaps, so it is not automatically contained.
  expect_equal(dropped$period, c(1L, 2L, NA))
  expect_equal(dropped$status[[3]], "tie_excluded")
  expect_equal(earliest$period[[3]], 2L)
  expect_equal(latest$period[[3]], 3L)
})

test_that("full-date maximum overlap is calculated in days", {
  dat <- data.frame(
    year = NA_real_,
    start = as.Date("1999-12-31"),
    end = as.Date("2000-01-10")
  )
  assigned <- assign_intervals(dat, ambiguous = "maxOverlap")

  expect_equal(assigned$period, 2L)
})

test_that("partial and wholly outside intervals are distinguished", {
  dat <- data.frame(
    year = NA_real_,
    start = c("1988", "1980", "2025"),
    end = c("1992", "1985", "2030")
  )
  assigned <- assign_intervals(dat, ambiguous = "maxOverlap")

  expect_equal(assigned$period, c(1L, NA, NA))
  expect_equal(
    assigned$status,
    c("ambiguous_assigned", "outside_periods", "outside_periods")
  )
})

test_that("gaps prevent automatic containment", {
  periods <- list(c(2000, 2002), 2003:2005)
  dat <- data.frame(year = NA_real_, start = "2000", end = "2002")

  dropped <- assign_intervals(dat, periods = periods, ambiguous = "drop")
  overlapped <- assign_intervals(
    dat, periods = periods, ambiguous = "maxOverlap"
  )

  expect_true(is.na(dropped$period))
  expect_equal(dropped$status, "ambiguous_excluded")
  expect_equal(overlapped$period, 1L)
})

test_that("chronological rules do not depend on supplied period order", {
  periods <- list(2010:2019, 1990:1999, 2000:2009)
  dat <- data.frame(
    year = NA_real_,
    start = c(NA, "1995", "1999"),
    end = c("2005", NA, "2000")
  )
  assigned <- assign_intervals(
    dat,
    periods = periods,
    ambiguous = "maxOverlap",
    endOnly = "earliest",
    startOnly = "latest",
    ties = "earliest"
  )

  expect_equal(assigned$period, c(2L, 1L, 2L))
})

test_that("open-ended records support drop and chronological allocation", {
  dat <- data.frame(
    year = NA_real_,
    start = c(NA, NA, "1995", "2025"),
    end = c("2005", "1980", NA, NA)
  )

  dropped <- assign_intervals(dat)
  allocated <- assign_intervals(
    dat, endOnly = "earliest", startOnly = "latest"
  )

  expect_equal(
    dropped$status,
    c("end_only_excluded", "outside_periods", "start_only_excluded", "outside_periods")
  )
  expect_equal(allocated$period, c(1L, NA, 3L, NA))
})

test_that("intervals take precedence and exact years provide fallback", {
  dat <- data.frame(
    year = c(2005, 2005, 2015, NA),
    start = c("1995", NA, NA, NA),
    end = c("1996", NA, "2005", NA)
  )
  assigned <- assign_intervals(dat, endOnly = "earliest")

  expect_equal(assigned$period, c(1L, 2L, 1L, NA))
  expect_equal(
    assigned$status,
    c("interval_unambiguous", "exact_assigned", "end_only_assigned", "missing_temporal")
  )
})

test_that("year-only strings, Date columns, and more than nine periods work", {
  periods <- lapply(2001:2011, function(year) year)
  dat <- data.frame(
    year = c("2001", "2010-06-14", "2011"),
    stringsAsFactors = FALSE
  )
  exact <- occAssess:::assignTemporalPeriods(dat, periods, "year")

  expect_equal(exact$period, c(1L, 10L, 11L))

  dates <- data.frame(
    start = as.Date(c("2001-02-01", "2010-01-01")),
    end = as.Date(c("2001-03-01", "2010-12-31"))
  )
  interval <- assign_intervals(dates, periods = periods, year = NULL)
  expect_equal(interval$period, c(1L, 10L))
})

test_that("malformed dates, impossible dates, and reversed intervals fail", {
  expect_error(
    assign_intervals(data.frame(year = NA, start = "05/06/2001", end = "2002")),
    "malformed dates"
  )
  expect_error(
    assign_intervals(data.frame(year = NA, start = "2001-02-30", end = "2002-01-01")),
    "impossible dates"
  )
  expect_error(
    assign_intervals(data.frame(year = NA, start = "2005", end = "2004")),
    "starts after it ends"
  )
})

test_that("period and option validation is authoritative", {
  dat <- data.frame(year = 2000, start = NA, end = NA)

  expect_error(
    occAssess:::assignTemporalPeriods(dat, list(2000:2001, 2001:2002)),
    "Periods must not overlap"
  )
  expect_error(
    occAssess:::assignTemporalPeriods(dat, list(numeric())),
    "non-empty finite numeric vector"
  )
  expect_error(
    assign_intervals(dat, ambiguous = "random"),
    "dateHandling\\$ambiguous"
  )
  expect_error(
    assign_intervals(dat, periods = list(2000.5)),
    "integer calendar years"
  )
})
