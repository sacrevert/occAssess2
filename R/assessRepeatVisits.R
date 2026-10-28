#' Summarise occurrence-record counts within grid cells
#'
#' Counts occurrence records in grid cells and divides each count by the number
#' of years specified in its period. Despite the historical function name, the
#' result is not a count of independent visits: no visit identifier is supplied
#' and coincident records are not deduplicated.
#'
#' @param dat A data frame containing numeric \code{lon}, \code{lat}, and
#'   an optional exact \code{year} column and an \code{identifier} column. Coordinates may be
#'   longitude/latitude or projected coordinates, but both coordinate columns
#'   and \code{spatRes} must use the same coordinate system and units.
#' @param periods A non-empty list of non-empty, non-overlapping numeric vectors
#'   defining the years in each period. Duplicate years, within or between
#'   periods, are not permitted. The number of entries in each vector is the
#'   denominator for that period.
#' @param spatRes A single positive, finite number giving the grid-cell width and
#'   height in the units of \code{lon} and \code{lat}.
#' @param dateHandling Optional uncertain-date configuration. When non-`NULL`,
#'   start/end columns may replace `year`; see Details.
#'
#' @details For a grid cell, identifier, and period, the returned value is
#'
#'   \code{number of occurrence rows / length(periods[[period]])}.
#'
#'   It is therefore occurrence records per grid cell per specified period-year.
#'   It can be used as a proxy for repeated spatial sampling, but it does not
#'   identify independent recording visits. In particular, the denominator is
#'   not the number of years in which a cell was actually sampled.
#'
#'   A single rectangular grid spanning the retained coordinates is constructed
#'   for all periods and identifiers. Cell boundaries are multiples of
#'   \code{spatRes} relative to coordinate zero; a coordinate exactly on a
#'   boundary is assigned to the cell beginning at that boundary. This makes
#'   every cell exactly \code{spatRes} units wide and high and makes grid
#'   alignment independent of row order, period, and identifier. Only records
#'   in the requested periods determine the grid.
#'
#'   In exact-year mode, records with missing or non-finite years or coordinates
#'   are removed with a warning. Interval-mode exclusions are reported compactly.
#'   Records outside all requested periods are excluded. Missing
#'   identifiers are not permitted. A defined period with no observations is
#'   retained in \code{data} with missing values for every grid cell. Unoccupied
#'   cells are likewise represented by \code{NA_real_}; the frequency polygon
#'   displays the distribution across occupied cells only. Values are not
#'   truncated in the returned data or plot.
#'
#'   To prevent accidental excessive allocation, an error is raised when the
#'   grid or the expanded period-by-identifier output would exceed ten million
#'   rows; increase \code{spatRes} in that case.
#'
#'   Temporal assignment follows the package-wide rules in the README and does
#'   not change the denominator: counts are still divided by
#'   `length(periods[[period]])`. Diagnostics are stored in the result's
#'   `temporalDiagnostics` attribute.
#'
#' @return A list with \code{data} and \code{plot}. \code{data} is a data frame
#'   with columns \code{vals}, \code{identifier}, and \code{Period}. It contains
#'   one row for every grid-cell, identifier, and period combination;
#'   \code{vals} contains occurrence records per grid cell per specified
#'   period-year and is missing for unoccupied cells. \code{plot} is a ggplot
#'   frequency polygon of non-missing \code{vals}, coloured by period and
#'   faceted by identifier.
#' @seealso [temporal-assignment]
#' @export

assessRepeatVisits <- function(dat, periods, spatRes, dateHandling = NULL) {
  yearColumn <- if ("year" %in% names(dat)) "year" else NULL
  required <- c("lon", "lat", "identifier", if (is.null(dateHandling)) "year")

  if (!is.data.frame(dat)) {
    stop("dat must be a data frame.")
  }
  if (any(!(required %in% names(dat)))) {
    stop("dat must contain lon, lat, identifier, and an exact year unless interval mode is used.")
  }
  if (nrow(dat) == 0) {
    stop("dat contains no records.")
  }
  if (!is.numeric(dat$lon) || !is.numeric(dat$lat)) {
    stop("lon and lat must be numeric.")
  }
  if (any(is.na(dat$identifier))) {
    stop("Missing values are not permitted in identifier.")
  }

  temporal <- assignTemporalPeriods(dat, periods, yearColumn, dateHandling)
  dat$.Period <- temporal$period
  dat$.temporalStatus <- temporal$status

  if (!is.numeric(spatRes) || length(spatRes) != 1 ||
      is.na(spatRes) || !is.finite(spatRes) || spatRes <= 0) {
    stop("spatRes must be a single positive, finite number.")
  }

  invalidCoordinates <- !is.finite(dat$lon) | !is.finite(dat$lat)
  if (any(invalidCoordinates)) {
    warning(sprintf(
      "Removing %d record%s with missing or non-finite coordinates.",
      sum(invalidCoordinates),
      if (sum(invalidCoordinates) == 1) "" else "s"
    ))
    dat <- dat[!invalidCoordinates, , drop = FALSE]
  }

  selection <- .finalizeTemporalAssignment(
    dat$.Period, dat$.temporalStatus, periods, dateHandling,
    exactMissingWarning = function(n) sprintf(
      "Removing %d record%s with a missing or non-finite year.",
      n, if (n == 1) "" else "s"
    )
  )
  diagnostics <- selection$diagnostics
  dat <- dat[selection$keep, , drop = FALSE]

  if (nrow(dat) == 0) {
    stop("No records fall within the requested periods.")
  }

  cellBin <- function(coordinate) {
    scaled <- coordinate / spatRes
    if (any(!is.finite(scaled))) {
      stop("Coordinates are too large relative to spatRes to construct a grid.")
    }
    nearestInteger <- round(scaled)
    tolerance <- pmin(
      sqrt(.Machine$double.eps),
      8 * .Machine$double.eps * pmax(1, abs(scaled))
    )
    onBoundary <- abs(scaled - nearestInteger) <= tolerance
    scaled[onBoundary] <- nearestInteger[onBoundary]
    floor(scaled)
  }

  xBin <- cellBin(dat$lon)
  yBin <- cellBin(dat$lat)
  if (any(!is.finite(xBin)) || any(!is.finite(yBin))) {
    stop("Coordinates are too large relative to spatRes to construct a grid.")
  }

  xRange <- range(xBin)
  yRange <- range(yBin)
  nCols <- xRange[2] - xRange[1] + 1
  nRows <- yRange[2] - yRange[1] + 1
  nCells <- nCols * nRows
  identifiers <- unique(dat$identifier)
  nOutputRows <- nCells * length(periods) * length(identifiers)
  maxRows <- 1e7

  if (!is.finite(nCells) || nCells > maxRows ||
      !is.finite(nOutputRows) || nOutputRows > maxRows) {
    stop(
      "The requested grid would exceed ten million output rows; increase spatRes."
    )
  }

  nCols <- as.integer(nCols)
  nCells <- as.integer(nCells)
  xCell <- as.integer(xBin - xRange[1] + 1)
  yCell <- as.integer(yRange[2] - yBin + 1)
  dat$.cell <- (yCell - 1L) * nCols + xCell
  dat$.period <- dat$.Period

  output <- lapply(identifiers, function(id) {
    periodOutput <- lapply(seq_along(periods), function(p) {
      selected <- dat$identifier == id & dat$.period == p
      counts <- tabulate(dat$.cell[selected], nbins = nCells)
      values <- rep(NA_real_, nCells)
      occupied <- counts > 0
      values[occupied] <- counts[occupied] / length(periods[[p]])

      data.frame(
        vals = values,
        identifier = rep(id, nCells),
        Period = rep(paste0("p", p), nCells)
      )
    })

    do.call("rbind", periodOutput)
  })

  out <- do.call("rbind", output)
  rownames(out) <- NULL

  plotData <- out[!is.na(out$vals), , drop = FALSE]
  plotData$Period <- factor(
    plotData$Period,
    levels = paste0("p", seq_along(periods))
  )

  p <- ggplot2::ggplot(
    data = plotData,
    ggplot2::aes(x = vals, colour = Period, group = Period)
  ) +
    ggplot2::geom_freqpoly(bins = 30) +
    ggplot2::xlab("Occurrence records per grid cell per year in period") +
    ggplot2::ylab("Number of occupied grid cells") +
    ggplot2::theme_linedraw() +
    ggplot2::scale_colour_discrete(drop = FALSE) +
    ggplot2::facet_wrap(~identifier)

  return(.addTemporalDiagnostics(
    list(data = out, plot = p), diagnostics
  ))
}
