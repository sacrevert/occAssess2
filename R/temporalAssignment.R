.validateTemporalPeriods <- function(periods, requireCalendarYears = FALSE) {
  if (!is.list(periods) || length(periods) == 0) {
    stop("periods must be a non-empty list of numeric year vectors.")
  }
  valid <- vapply(
    periods,
    function(period) {
      is.numeric(period) && length(period) > 0 &&
        all(is.finite(period)) && !anyDuplicated(period)
    },
    logical(1)
  )
  if (!all(valid)) {
    stop("Each period must be a non-empty finite numeric vector without duplicate years.")
  }
  if (anyDuplicated(unlist(periods, use.names = FALSE))) {
    stop("Periods must not overlap.")
  }
  if (requireCalendarYears) {
    years <- unlist(periods, use.names = FALSE)
    if (any(years != floor(years)) || any(years < 1) || any(years > 9999)) {
      stop("Interval date handling requires periods containing integer calendar years from 1 to 9999.")
    }
  }
  invisible(periods)
}

.dateHandlingSettings <- function(dateHandling) {
  defaults <- list(
    start = NULL,
    end = NULL,
    ambiguous = "drop",
    endOnly = "drop",
    startOnly = "drop",
    ties = "drop"
  )
  if (is.null(dateHandling)) {
    defaults$enabled <- FALSE
    return(defaults)
  }
  if (!is.list(dateHandling) || is.null(names(dateHandling)) ||
      any(names(dateHandling) == "") || anyDuplicated(names(dateHandling))) {
    stop("dateHandling must be NULL or a named list.")
  }
  unknown <- setdiff(names(dateHandling), names(defaults))
  if (length(unknown) > 0) {
    stop("Unsupported dateHandling field(s): ", paste(unknown, collapse = ", "), ".")
  }
  settings <- utils::modifyList(defaults, dateHandling)
  validateColumn <- function(value, field) {
    if (!is.null(value) &&
        (!is.character(value) || length(value) != 1 || is.na(value) ||
         !nzchar(value))) {
      stop("dateHandling$", field, " must be NULL or one column name.")
    }
  }
  validateColumn(settings$start, "start")
  validateColumn(settings$end, "end")
  if (is.null(settings$start) && is.null(settings$end)) {
    stop("dateHandling must name at least one start or end column.")
  }
  choices <- list(
    ambiguous = c("drop", "start", "end", "maxOverlap"),
    endOnly = c("drop", "earliest"),
    startOnly = c("drop", "latest"),
    ties = c("drop", "earliest", "latest")
  )
  for (field in names(choices)) {
    value <- settings[[field]]
    if (!is.character(value) || length(value) != 1 || is.na(value) ||
        !(value %in% choices[[field]])) {
      stop(
        "dateHandling$", field, " must be one of ",
        paste(sprintf("'%s'", choices[[field]]), collapse = ", "), "."
      )
    }
  }
  settings$enabled <- TRUE
  settings
}

.parseTemporalValues <- function(values, columnName) {
  n <- length(values)
  dates <- as.Date(rep(NA_character_, n))
  years <- rep(NA_integer_, n)
  precision <- rep(NA_character_, n)

  if (is.logical(values) && all(is.na(values))) {
    return(.missingTemporalValues(n))
  }

  if (inherits(values, "Date")) {
    present <- !is.na(values)
    dates[present] <- values[present]
    years[present] <- as.integer(format(values[present], "%Y"))
    precision[present] <- "date"
    return(list(date = dates, year = years, precision = precision,
                present = present))
  }

  if (is.numeric(values)) {
    present <- !is.na(values)
    invalid <- present & (!is.finite(values) | values != floor(values) |
                            values < 1 | values > 9999)
    if (any(invalid)) {
      stop(
        "Column '", columnName,
        "' contains invalid calendar years at row(s): ",
        paste(which(invalid), collapse = ", "), "."
      )
    }
    years[present] <- as.integer(values[present])
    dates[present] <- as.Date(sprintf("%04d-01-01", years[present]))
    precision[present] <- "year"
    return(list(date = dates, year = years, precision = precision,
                present = present))
  }

  if (is.character(values) || is.factor(values)) {
    text <- trimws(as.character(values))
    present <- !is.na(text) & nzchar(text)
    yearOnly <- present & grepl("^[0-9]{4}$", text)
    fullDate <- present & grepl("^[0-9]{4}-[0-9]{2}-[0-9]{2}$", text)
    malformed <- present & !yearOnly & !fullDate
    if (any(malformed)) {
      stop(
        "Column '", columnName,
        "' contains malformed dates at row(s): ",
        paste(which(malformed), collapse = ", "),
        ". Use calendar years or ISO dates in YYYY-MM-DD form."
      )
    }
    if (any(yearOnly)) {
      parsedYears <- as.integer(text[yearOnly])
      invalidYears <- parsedYears < 1 | parsedYears > 9999
      if (any(invalidYears)) {
        stop("Column '", columnName, "' contains invalid calendar years.")
      }
      years[yearOnly] <- parsedYears
      dates[yearOnly] <- as.Date(sprintf("%04d-01-01", parsedYears))
      precision[yearOnly] <- "year"
    }
    if (any(fullDate)) {
      parsed <- suppressWarnings(as.Date(text[fullDate], format = "%Y-%m-%d"))
      canonical <- ifelse(
        is.na(parsed), NA_character_, format(parsed, "%Y-%m-%d")
      )
      invalidDates <- is.na(parsed) | canonical != text[fullDate]
      if (any(invalidDates)) {
        rows <- which(fullDate)[invalidDates]
        stop(
          "Column '", columnName,
          "' contains impossible dates at row(s): ",
          paste(rows, collapse = ", "), "."
        )
      }
      dates[fullDate] <- parsed
      years[fullDate] <- as.integer(format(parsed, "%Y"))
      precision[fullDate] <- "date"
    }
    return(list(date = dates, year = years, precision = precision,
                present = present))
  }

  stop(
    "Column '", columnName,
    "' must contain numeric years, Date values, or ISO date strings."
  )
}

.missingTemporalValues <- function(n) {
  list(
    date = as.Date(rep(NA_character_, n)),
    year = rep(NA_integer_, n),
    precision = rep(NA_character_, n),
    present = rep(FALSE, n)
  )
}

.periodContainingYear <- function(year, periods) {
  matches <- which(vapply(periods, function(period) year %in% period, logical(1)))
  if (length(matches) == 0) NA_integer_ else matches[[1]]
}

.intervalOverlap <- function(startDate, endDate, startYear, endYear,
                             yearPrecision, periods) {
  if (yearPrecision) {
    return(vapply(
      periods,
      function(period) sum(period >= startYear & period <= endYear),
      numeric(1)
    ))
  }
  vapply(periods, function(period) {
    sum(vapply(period, function(periodYear) {
      periodStart <- as.Date(sprintf("%04d-01-01", as.integer(periodYear)))
      periodEnd <- as.Date(sprintf("%04d-12-31", as.integer(periodYear)))
      overlapStart <- max(startDate, periodStart)
      overlapEnd <- min(endDate, periodEnd)
      if (overlapStart > overlapEnd) 0 else as.numeric(overlapEnd - overlapStart) + 1
    }, numeric(1)))
  }, numeric(1))
}

.containedPeriod <- function(startYear, endYear, periods) {
  span <- endYear - startYear + 1
  matches <- which(vapply(periods, function(period) {
    sum(period >= startYear & period <= endYear) == span
  }, logical(1)))
  if (length(matches) == 1) matches[[1]] else NA_integer_
}

.temporalDiagnostics <- function(period, status, periods) {
  count <- function(value) sum(status %in% value)
  periodCounts <- tabulate(period, nbins = length(periods))
  names(periodCounts) <- paste0("p", seq_along(periods))
  list(
    records = length(period),
    assignedUnambiguously = count(c("exact_assigned", "interval_unambiguous")),
    exactAssigned = count("exact_assigned"),
    intervalUnambiguous = count("interval_unambiguous"),
    ambiguousAssigned = count("ambiguous_assigned"),
    ambiguousExcluded = count("ambiguous_excluded"),
    outsidePeriods = count("outside_periods"),
    missingTemporal = count("missing_temporal"),
    endOnlyAssigned = count("end_only_assigned"),
    endOnlyExcluded = count("end_only_excluded"),
    startOnlyAssigned = count("start_only_assigned"),
    startOnlyExcluded = count("start_only_excluded"),
    tiesExcluded = count("tie_excluded"),
    periodCounts = periodCounts
  )
}

# Internal authoritative temporal-assignment implementation.
assignTemporalPeriods <- function(data, periods, year = "year",
                                  dateHandling = NULL) {
  if (!is.data.frame(data)) {
    stop("data must be a data frame.")
  }
  settings <- .dateHandlingSettings(dateHandling)
  .validateTemporalPeriods(periods, requireCalendarYears = settings$enabled)
  n <- nrow(data)

  if (is.null(year)) {
    exact <- .missingTemporalValues(n)
  } else {
    if (!is.character(year) || length(year) != 1 || is.na(year) ||
        !(year %in% names(data))) {
      stop("year must be NULL or name a column in data.")
    }
    exact <- .parseTemporalValues(data[[year]], year)
  }
  if (!settings$enabled && is.null(year)) {
    stop("year cannot be NULL when dateHandling is NULL.")
  }

  parseBound <- function(column) {
    if (is.null(column)) return(.missingTemporalValues(n))
    if (!(column %in% names(data))) {
      stop("dateHandling column '", column, "' does not exist in data.")
    }
    .parseTemporalValues(data[[column]], column)
  }
  start <- if (settings$enabled) parseBound(settings$start) else
    .missingTemporalValues(n)
  end <- if (settings$enabled) parseBound(settings$end) else
    .missingTemporalValues(n)

  startDate <- start$date
  endDate <- end$date
  endYearPrecision <- end$present & end$precision == "year"
  endDate[endYearPrecision] <- as.Date(sprintf(
    "%04d-12-31", end$year[endYearPrecision]
  ))
  reversed <- start$present & end$present & startDate > endDate
  if (any(reversed)) {
    stop(
      "Temporal interval starts after it ends at row(s): ",
      paste(which(reversed), collapse = ", "), "."
    )
  }

  assigned <- rep(NA_integer_, n)
  status <- rep("missing_temporal", n)
  chronologicalStart <- vapply(periods, min, numeric(1))
  chronologicalEnd <- vapply(periods, max, numeric(1))

  for (row in seq_len(n)) {
    hasStart <- start$present[[row]]
    hasEnd <- end$present[[row]]

    if (hasStart && hasEnd) {
      contained <- .containedPeriod(start$year[[row]], end$year[[row]], periods)
      if (!is.na(contained)) {
        assigned[[row]] <- contained
        status[[row]] <- "interval_unambiguous"
        next
      }

      overlap <- .intervalOverlap(
        startDate[[row]], endDate[[row]], start$year[[row]], end$year[[row]],
        start$precision[[row]] == "year" && end$precision[[row]] == "year",
        periods
      )
      if (all(overlap == 0)) {
        status[[row]] <- "outside_periods"
        next
      }
      if (settings$ambiguous == "drop") {
        status[[row]] <- "ambiguous_excluded"
        next
      }
      if (settings$ambiguous == "start") {
        chosen <- .periodContainingYear(start$year[[row]], periods)
      } else if (settings$ambiguous == "end") {
        chosen <- .periodContainingYear(end$year[[row]], periods)
      } else {
        candidates <- which(overlap == max(overlap))
        if (length(candidates) > 1) {
          if (settings$ties == "drop") {
            status[[row]] <- "tie_excluded"
            next
          }
          chosen <- if (settings$ties == "earliest") {
            candidates[[which.min(chronologicalStart[candidates])]]
          } else {
            candidates[[which.max(chronologicalEnd[candidates])]]
          }
        } else {
          chosen <- candidates[[1]]
        }
      }
      if (is.na(chosen)) {
        status[[row]] <- "ambiguous_excluded"
      } else {
        assigned[[row]] <- chosen
        status[[row]] <- "ambiguous_assigned"
      }
      next
    }

    if (hasEnd) {
      compatible <- which(vapply(
        periods,
        function(period) min(period) <= end$year[[row]],
        logical(1)
      ))
      if (length(compatible) == 0) {
        status[[row]] <- "outside_periods"
      } else if (settings$endOnly == "earliest") {
        chosen <- compatible[[which.min(chronologicalStart[compatible])]]
        assigned[[row]] <- chosen
        status[[row]] <- "end_only_assigned"
      } else {
        status[[row]] <- "end_only_excluded"
      }
      next
    }

    if (hasStart) {
      compatible <- which(vapply(
        periods,
        function(period) max(period) >= start$year[[row]],
        logical(1)
      ))
      if (length(compatible) == 0) {
        status[[row]] <- "outside_periods"
      } else if (settings$startOnly == "latest") {
        chosen <- compatible[[which.max(chronologicalEnd[compatible])]]
        assigned[[row]] <- chosen
        status[[row]] <- "start_only_assigned"
      } else {
        status[[row]] <- "start_only_excluded"
      }
      next
    }

    if (exact$present[[row]]) {
      chosen <- .periodContainingYear(exact$year[[row]], periods)
      if (is.na(chosen)) {
        status[[row]] <- "outside_periods"
      } else {
        assigned[[row]] <- chosen
        status[[row]] <- "exact_assigned"
      }
    }
  }

  diagnostics <- .temporalDiagnostics(assigned, status, periods)
  structure(
    list(
      period = assigned,
      status = status,
      diagnostics = diagnostics,
      settings = settings
    ),
    class = "occAssessTemporalAssignment"
  )
}

.finalizeTemporalAssignment <- function(period, status, periods,
                                        dateHandling = NULL,
                                        exactMissingWarning = NULL) {
  diagnostics <- .temporalDiagnostics(period, status, periods)
  if (is.null(dateHandling) && diagnostics$missingTemporal > 0 &&
      !is.null(exactMissingWarning)) {
    message <- if (is.function(exactMissingWarning)) {
      exactMissingWarning(diagnostics$missingTemporal)
    } else {
      exactMissingWarning
    }
    warning(message, call. = FALSE)
  }
  excluded <- diagnostics$ambiguousExcluded + diagnostics$missingTemporal +
    diagnostics$endOnlyExcluded + diagnostics$startOnlyExcluded +
    diagnostics$tiesExcluded
  if (is.null(dateHandling)) {
    excluded <- excluded - diagnostics$missingTemporal
  }
  assumed <- diagnostics$ambiguousAssigned + diagnostics$endOnlyAssigned +
    diagnostics$startOnlyAssigned
  if (excluded > 0 || assumed > 0) {
    details <- character()
    if (assumed > 0) details <- c(details, paste(assumed, "assigned by allocation rules"))
    if (excluded > 0) details <- c(details, paste(excluded, "excluded as temporally unassignable"))
    warning(
      "Temporal assignment: ", paste(details, collapse = "; "), ".",
      call. = FALSE
    )
  }
  list(keep = !is.na(period), diagnostics = diagnostics)
}

.addTemporalDiagnostics <- function(result, diagnostics) {
  attr(result, "temporalDiagnostics") <- diagnostics
  result
}
