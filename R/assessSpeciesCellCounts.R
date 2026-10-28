#' Count records and unique species-grid-cell combinations
#'
#' Summarises raw occurrence-record volume and the number of distinct
#' species-grid-cell combinations for each identifier and user-supplied period.
#'
#' @param dat A data frame containing species, coordinate,
#'   spatial-uncertainty, identifier, and exact-year and/or interval-bound
#'   columns.
#' @param species Character string naming the species-identity column in
#'   `dat`.
#' @param x Character string naming the x-coordinate column in `dat`.
#' @param y Character string naming the y-coordinate column in `dat`.
#' @param year Character string naming the exact-year column, or `NULL` in
#'   interval mode without an exact-year fallback.
#' @param spatialUncertainty Character string naming the spatial-uncertainty
#'   column in `dat`.
#' @param identifier Character string naming the grouping column in `dat`.
#'   Counts are calculated separately for each identifier.
#' @param periods A non-empty list of non-empty, non-overlapping numeric year
#'   vectors. Periods are reported in the supplied order, including periods
#'   containing no records.
#' @param res A single positive, finite number giving grid-cell width and height
#'   in the units of the supplied coordinates. Use degrees for longitude and
#'   latitude, or the corresponding linear units for projected coordinates.
#' @param maxSpatUncertainty Optional single finite maximum spatial uncertainty.
#'   When supplied, records with missing uncertainty or uncertainty greater
#'   than this value are excluded from all three counts.
#' @param dateHandling Optional uncertain-date configuration. `NULL` preserves
#'   exact-year assignment; see Details.
#'
#' @details For each identifier-period combination, `nRecords` is the number of
#'   rows remaining after period selection and optional spatial-uncertainty
#'   filtering. It includes rows whose species or coordinates cannot be used
#'   for spatial aggregation. `nEligibleRecords` is the subset having a
#'   non-missing, non-blank species identity and finite x and y coordinates.
#'   `nSpeciesCells` is the number of distinct species identity and grid-cell
#'   pairs among those eligible rows. Repeated records of a species in one cell
#'   count once within a period, including repetitions in different years;
#'   occurrences in another cell or another period count separately.
#'
#'   Cells are aligned to multiples of `res` relative to coordinate zero, using
#'   the same convention as [assessSpatialCov()]. A coordinate exactly on a
#'   boundary belongs to the cell beginning at that boundary. The convention is
#'   shared across every identifier and period and supports negative,
#'   geographic, projected, and non-integer coordinates.
#'
#'   Temporal assignment follows the package-wide rules documented in the
#'   README, and every occurrence contributes to at most one period. Rows
#'   without an assigned period are excluded before counting. Missing identifiers are not
#'   permitted. Missing or non-finite coordinates and missing or blank species
#'   identities remain in `nRecords` but not `nEligibleRecords` or
#'   `nSpeciesCells`. Missing spatial uncertainty is retained when
#'   `maxSpatUncertainty` is `NULL` and excluded when a threshold is supplied.
#'   If no eligible row remains in any requested period, the function stops
#'   because no species-cell statistic can be calculated.
#'
#'   This statistic differs from species richness because one species can
#'   contribute in several cells, and from occupied-cell count because several
#'   species can contribute in one cell. It is less sensitive than raw record
#'   count to repeated recording of the same species and cell, but it is not a
#'   count of independent visits or a direct measure of effort or data quality.
#'   Finer grid resolutions will generally produce at least as many
#'   species-cell combinations as coarser resolutions, subject to grid
#'   alignment.
#'
#' @return A list with `data` and `plot`. `data` contains one row for every
#'   represented identifier and requested period, with columns `identifier`,
#'   `Period`, `nRecords`, `nEligibleRecords`, and `nSpeciesCells`. `Period` is
#'   the position of the period in `periods`. `plot` is a ggplot object comparing
#'   raw record and species-cell counts; eligible-record counts remain available
#'   in `data` but are not plotted.
#' @seealso [assessRecordNumber()], [assessSpeciesNumber()],
#'   [assessSpatialCov()], and [temporal-assignment]
#' @examples
#' records <- data.frame(
#'   taxon = c("A", "A", "B", "A", NA),
#'   longitude = c(0.1, 0.1, 0.1, 1.1, 2.1),
#'   latitude = rep(0.1, 5),
#'   year = rep(2000, 5),
#'   uncertainty = rep(0, 5),
#'   source = rep("survey", 5)
#' )
#' result <- assessSpeciesCellCounts(
#'   records, species = "taxon", x = "longitude", y = "latitude",
#'   year = "year", spatialUncertainty = "uncertainty",
#'   identifier = "source", periods = list(2000), res = 1
#' )
#' result$data
#' result$plot
#' @export
assessSpeciesCellCounts <- function(dat,
                                    species,
                                    x,
                                    y,
                                    year,
                                    spatialUncertainty,
                                    identifier,
                                    periods,
                                    res,
                                    maxSpatUncertainty = NULL,
                                    dateHandling = NULL) {
  required <- c(species, x, y, spatialUncertainty, identifier, year)
  if (any(!(required %in% colnames(dat)))) {
    stop("You have specified columns that do not exist in dat.")
  }
  temporal <- assignTemporalPeriods(dat, periods, year, dateHandling)
  if (!is.numeric(res) || length(res) != 1 || is.na(res) ||
      !is.finite(res) || res <= 0) {
    stop("res must be a single positive, finite number.")
  }
  if (!is.null(maxSpatUncertainty) &&
      (!is.numeric(maxSpatUncertainty) ||
       length(maxSpatUncertainty) != 1 ||
       is.na(maxSpatUncertainty) || !is.finite(maxSpatUncertainty))) {
    stop("maxSpatUncertainty must be a single finite number.")
  }

  dat <- createData(
    data = dat, species, x, y, year, spatialUncertainty, identifier
  )
  dat$Period <- temporal$period
  dat$.temporalStatus <- temporal$status
  if (nrow(dat) == 0) {
    stop("dat contains no records.")
  }
  if (!is.numeric(dat$x) || !is.numeric(dat$y)) {
    stop("Coordinate columns must be numeric.")
  }
  if (any(is.na(dat$identifier))) {
    stop("Missing values are not permitted in identifier.")
  }

  if (!is.null(maxSpatUncertainty)) {
    if (!is.numeric(dat$spatialUncertainty)) {
      stop("spatialUncertainty must be numeric when maxSpatUncertainty is supplied.")
    }
    keep <- !is.na(dat$spatialUncertainty) &
      is.finite(dat$spatialUncertainty) &
      dat$spatialUncertainty <= maxSpatUncertainty
    dat <- dat[keep, , drop = FALSE]
  }
  if (nrow(dat) == 0) {
    stop("No records remain after spatial-uncertainty filtering.")
  }

  selection <- .finalizeTemporalAssignment(
    dat$Period, dat$.temporalStatus, periods, dateHandling,
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

  identifiers <- unique(dat$identifier)
  identifierIndex <- match(dat$identifier, identifiers)
  nPeriods <- length(periods)
  groupIndex <- (identifierIndex - 1L) * nPeriods + dat$Period
  nGroups <- length(identifiers) * nPeriods

  usableSpecies <- !is.na(dat$species)
  if (is.character(dat$species) || is.factor(dat$species)) {
    usableSpecies <- usableSpecies & nzchar(trimws(as.character(dat$species)))
  }
  finiteCoordinates <- is.finite(dat$x) & is.finite(dat$y)
  eligible <- usableSpecies & finiteCoordinates
  if (!any(eligible)) {
    stop("No records with usable species identities and finite coordinates remain in the requested periods.")
  }

  cellBins <- .occAssessCellBins(dat$x[eligible], dat$y[eligible], res)
  combinations <- data.frame(
    group = groupIndex[eligible],
    species = dat$species[eligible],
    xBin = cellBins$x,
    yBin = cellBins$y,
    check.names = FALSE
  )
  distinctCombinations <- combinations[!duplicated(combinations), , drop = FALSE]

  resultData <- data.frame(
    identifier = rep(identifiers, each = nPeriods),
    Period = rep(seq_len(nPeriods), times = length(identifiers)),
    nRecords = tabulate(groupIndex, nbins = nGroups),
    nEligibleRecords = tabulate(groupIndex[eligible], nbins = nGroups),
    nSpeciesCells = tabulate(distinctCombinations$group, nbins = nGroups)
  )

  measures <- c("Raw records", "Species-cell combinations")
  plotData <- rbind(
    data.frame(
      identifier = resultData$identifier,
      Period = resultData$Period,
      Measure = measures[[1]],
      Count = resultData$nRecords
    ),
    data.frame(
      identifier = resultData$identifier,
      Period = resultData$Period,
      Measure = measures[[2]],
      Count = resultData$nSpeciesCells
    )
  )
  plotData$Measure <- factor(plotData$Measure, levels = measures)

  p <- ggplot2::ggplot(
    plotData,
    ggplot2::aes(
      x = Period, y = Count, colour = Measure, linetype = Measure,
      shape = Measure, group = interaction(identifier, Measure)
    )
  ) +
    ggplot2::geom_line() +
    ggplot2::geom_point() +
    ggplot2::facet_wrap(~identifier) +
    ggplot2::scale_x_continuous(breaks = seq_along(periods)) +
    ggplot2::theme_linedraw() +
    ggplot2::labs(
      x = "Period", y = "Count", colour = "Measure",
      linetype = "Measure", shape = "Measure"
    )

  .addTemporalDiagnostics(
    list(data = resultData, plot = p), diagnostics
  )
}
