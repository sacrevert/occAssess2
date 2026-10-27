#' Map gridded spatial coverage
#'
#' Counts occurrence records on a common spatial grid and maps either those
#' counts, the number of sampled periods, or cells sampled in at least a chosen
#' number of periods.
#'
#' @param dat A data frame containing species, coordinate, year,
#'   spatial-uncertainty, and identifier columns.
#' @param species Character string naming the species column in dat.
#' @param x Character string naming the x-coordinate column in dat.
#' @param y Character string naming the y-coordinate column in dat.
#' @param year Character string naming the numeric year column in dat.
#' @param spatialUncertainty Character string naming the spatial-uncertainty
#'   column in dat.
#' @param identifier Character string naming the grouping column in dat. Results
#'   are calculated separately for each identifier.
#' @param res A single positive, finite number giving grid-cell width and height
#'   in the units of the supplied coordinates. Geographic coordinates therefore
#'   use degrees, while projected coordinates commonly use metres.
#' @param logCount Logical. For density output, whether plotted and returned cell
#'   counts are transformed with log10. It has no effect on nPeriods or overlap;
#'   period presence is always calculated from untransformed counts.
#' @param countries Character vector of country names accepted by
#'   ggplot2::map_data("world"), or a named list with one entry per identifier.
#'   Use only with longitude/latitude coordinates. Optional.
#' @param shp A spatial polygon object, or a named list with one object per
#'   identifier, used only to add boundaries to plots. Do not supply both
#'   countries and shp.
#' @param periods A non-empty list of non-empty, non-overlapping numeric vectors
#'   defining periods in the requested order. Every supplied period is retained,
#'   including periods containing no observations.
#' @param maxSpatUncertainty Optional maximum permitted spatial uncertainty.
#'   Records with missing uncertainty or values above this threshold are removed.
#' @param output Character string selecting "density", "nPeriods", or
#'   "overlap". See Details.
#' @param minPeriods For overlap output, a single integer from one to
#'   length(periods) giving the minimum sampled-period count. The default NULL
#'   requires sampling in every user-supplied period.
#' @param returnRaster Logical. Whether to return the rasters underpinning the
#'   maps in addition to the named list of plots.
#'
#' @details A single rectangular grid spanning all retained coordinates is used
#'   for every identifier and requested period. Cell boundaries are multiples of
#'   res relative to coordinate zero. A coordinate on a boundary is placed in
#'   the cell beginning at that boundary. This zero-aligned convention makes cell
#'   width, height, and alignment independent of record order, identifier, and
#'   period, but can change results from earlier package versions whose grid
#'   began at the minimum observed coordinate and did not use res exactly.
#'
#'   Records with missing or non-finite years or coordinates are removed with a
#'   warning. Records outside periods are excluded and do not determine the grid.
#'   Missing identifiers are not permitted. A grid or expanded
#'   identifier-by-period stack exceeding ten million cell values is rejected to
#'   prevent accidental excessive allocation.
#'
#'   The output modes are:
#'
#'   - "density": number of occurrence rows in each grid cell, separately for
#'     every identifier and period. Repeated coordinates count repeatedly. This
#'     is records per grid cell, not an area-standardised or spatially continuous
#'     density. Unoccupied cells are NA. If logCount is TRUE, log10(count) is
#'     used in both the map and returned density rasters; thus a raw count of one
#'     is returned as zero. Period facets share one fill scale for comparison,
#'     but their values are neither pooled nor summed across periods.
#'   - "nPeriods": number of supplied periods containing at least one record in
#'     each cell, calculated from raw count presence. Cells never sampled in any
#'     period are NA.
#'   - "overlap": logical coverage. Qualifying cells are TRUE and other cells
#'     are NA; qualification means that the raw-count sampled-period total is at
#'     least minPeriods.
#'
#' @return By default, a named list of ggplot objects, one per identifier.
#'   Density plots contain facets p1, p2, and so on in the supplied period order.
#'   If returnRaster is TRUE, a list with plots and rasters is returned. rasters
#'   is a RasterStack ordered first by identifier and then, for density output,
#'   by supplied period. Density layer names have the form id.ID.period.pN,
#'   where ID is the identifier and N is the period index; other modes use
#'   id.ID. Names are made syntactically valid and unique.
#' @seealso assessSpatialBias(), which measures deviation from a random spatial
#'   distribution.
#' @importFrom rasterVis gplot
#' @importFrom ggplot2 map_data
#' @importFrom maps map
#' @export

assessSpatialCov <- function(dat,
                             species,
                             x,
                             y,
                             year,
                             spatialUncertainty,
                             identifier,
                             res,
                             logCount = FALSE,
                             countries = NULL,
                             shp = NULL,
                             periods,
                             maxSpatUncertainty = NULL,
                             output = "density",
                             minPeriods = NULL,
                             returnRaster = FALSE) {
  required <- c(species, x, y, year, spatialUncertainty, identifier)
  if (any(!(required %in% colnames(dat)))) {
    stop("You have specified columns that do not exist in dat.")
  }

  if (!is.numeric(res) || length(res) != 1 || is.na(res) ||
      !is.finite(res) || res <= 0) {
    stop("res must be a single positive, finite number.")
  }
  if (!is.logical(logCount) || length(logCount) != 1 || is.na(logCount)) {
    stop("logCount must be TRUE or FALSE.")
  }
  if (!is.logical(returnRaster) || length(returnRaster) != 1 ||
      is.na(returnRaster)) {
    stop("returnRaster must be TRUE or FALSE.")
  }
  if (!is.character(output) || length(output) != 1 ||
      !(output %in% c("density", "nPeriods", "overlap"))) {
    stop("output must be one of 'density', 'nPeriods', or 'overlap'.")
  }
  if (!is.list(periods) || length(periods) == 0) {
    stop("periods must be a non-empty list of numeric year vectors.")
  }
  validPeriods <- vapply(
    periods,
    function(period) {
      is.numeric(period) && length(period) > 0 &&
        all(is.finite(period)) && !anyDuplicated(period)
    },
    logical(1)
  )
  if (!all(validPeriods)) {
    stop("Each period must be a non-empty finite numeric vector without duplicate years.")
  }
  if (anyDuplicated(unlist(periods, use.names = FALSE))) {
    stop("Periods must not overlap.")
  }

  if (!is.null(minPeriods)) {
    if (!is.numeric(minPeriods) || length(minPeriods) != 1 ||
        is.na(minPeriods) || !is.finite(minPeriods) ||
        minPeriods != floor(minPeriods) || minPeriods < 1 ||
        minPeriods > length(periods)) {
      stop("minPeriods must be a single integer between 1 and length(periods).")
    }
  }
  overlapThreshold <- if (is.null(minPeriods)) length(periods) else minPeriods

  if (!is.null(shp) && !is.null(countries)) {
    stop("Only one of countries and shp may be specified.")
  }
  if (is.list(countries) || is.list(shp)) {
    message(
      "A list was supplied for shp or countries; list names should match identifiers."
    )
  }

  dat <- createData(data = dat,
                    species,
                    x,
                    y,
                    year,
                    spatialUncertainty,
                    identifier)

  if (nrow(dat) == 0) {
    stop("dat contains no records.")
  }
  if (!is.numeric(dat$x) || !is.numeric(dat$y) || !is.numeric(dat$year)) {
    stop("Coordinate and year columns must be numeric.")
  }
  if (any(is.na(dat$identifier))) {
    stop("Missing values are not permitted in identifier.")
  }

  if (!is.null(maxSpatUncertainty)) {
    if (!is.numeric(maxSpatUncertainty) || length(maxSpatUncertainty) != 1 ||
        is.na(maxSpatUncertainty) || !is.finite(maxSpatUncertainty)) {
      stop("maxSpatUncertainty must be a single finite number.")
    }
    keep <- !is.na(dat$spatialUncertainty) &
      dat$spatialUncertainty <= maxSpatUncertainty
    dat <- dat[keep, , drop = FALSE]
  }
  if (nrow(dat) == 0) {
    stop("No records remain after spatial-uncertainty filtering.")
  }

  invalidYear <- !is.finite(dat$year)
  if (any(invalidYear)) {
    warning(sprintf(
      "Removing %d record%s with a missing or non-finite year.",
      sum(invalidYear),
      if (sum(invalidYear) == 1) "" else "s"
    ))
    dat <- dat[!invalidYear, , drop = FALSE]
  }

  invalidCoordinates <- !is.finite(dat$x) | !is.finite(dat$y)
  if (any(invalidCoordinates)) {
    warning(sprintf(
      "Removing %d record%s with missing or non-finite coordinates.",
      sum(invalidCoordinates),
      if (sum(invalidCoordinates) == 1) "" else "s"
    ))
    dat <- dat[!invalidCoordinates, , drop = FALSE]
  }

  periodIndex <- rep(NA_integer_, nrow(dat))
  for (period in seq_along(periods)) {
    periodIndex[dat$year %in% periods[[period]]] <- period
  }
  retained <- !is.na(periodIndex)
  dat <- dat[retained, , drop = FALSE]
  periodIndex <- periodIndex[retained]
  if (nrow(dat) == 0) {
    stop("No records fall within the requested periods.")
  }

  cellBin <- function(coordinate) {
    scaled <- coordinate / res
    if (any(!is.finite(scaled))) {
      stop("Coordinates are too large relative to res to construct a grid.")
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

  xBin <- cellBin(dat$x)
  yBin <- cellBin(dat$y)
  xRange <- range(xBin)
  yRange <- range(yBin)
  nCols <- xRange[2] - xRange[1] + 1
  nRows <- yRange[2] - yRange[1] + 1
  nCells <- nCols * nRows
  identifiers <- unique(dat$identifier)
  identifierLabels <- as.character(identifiers)
  nExpandedCells <- nCells * length(periods) * length(identifiers)
  maxCells <- 1e7

  if (!is.finite(nCells) || nCells > maxCells ||
      !is.finite(nExpandedCells) || nExpandedCells > maxCells) {
    stop("The requested grid would exceed ten million cell values; increase res.")
  }

  nCols <- as.integer(nCols)
  nRows <- as.integer(nRows)
  nCells <- as.integer(nCells)
  xCell <- as.integer(xBin - xRange[1] + 1)
  yCell <- as.integer(yRange[2] - yBin + 1)
  dat$.cell <- (yCell - 1L) * nCols + xCell
  dat$.period <- periodIndex
  periodLabels <- paste0("p", seq_along(periods))

  template <- raster::raster(
    ncol = nCols,
    nrow = nRows,
    xmn = xRange[1] * res,
    xmx = (xRange[2] + 1) * res,
    ymn = yRange[1] * res,
    ymx = (yRange[2] + 1) * res
  )

  rawByIdentifier <- lapply(seq_along(identifiers), function(idIndex) {
    id <- identifiers[[idIndex]]
    layers <- lapply(seq_along(periods), function(period) {
      selected <- dat$identifier == id & dat$.period == period
      counts <- tabulate(dat$.cell[selected], nbins = nCells)
      values <- rep(NA_real_, nCells)
      occupied <- counts > 0
      values[occupied] <- counts[occupied]
      layer <- raster::setValues(template, values)
      names(layer) <- periodLabels[[period]]
      layer
    })
    result <- raster::stack(layers)
    names(result) <- periodLabels
    result
  })
  names(rawByIdentifier) <- identifierLabels

  rastersByIdentifier <- lapply(seq_along(rawByIdentifier), function(idIndex) {
    raw <- rawByIdentifier[[idIndex]]

    if (output == "density") {
      result <- if (logCount) log10(raw) else raw
      names(result) <- periodLabels
      return(result)
    }

    presence <- vapply(
      seq_len(raster::nlayers(raw)),
      function(layer) !is.na(raster::getValues(raw[[layer]])),
      logical(nCells)
    )
    if (length(periods) == 1) {
      presence <- matrix(presence, ncol = 1)
    }
    sampledPeriods <- rowSums(presence)
    sampledPeriods[sampledPeriods == 0] <- NA_real_

    if (output == "nPeriods") {
      result <- raster::setValues(template, sampledPeriods)
    } else {
      qualifies <- !is.na(sampledPeriods) & sampledPeriods >= overlapThreshold
      values <- rep(NA, nCells)
      values[qualifies] <- TRUE
      result <- raster::setValues(template, values)
    }
    names(result) <- identifierLabels[[idIndex]]
    result
  })
  names(rastersByIdentifier) <- identifierLabels

  boundaryForIdentifier <- function(identifierLabel) {
    if (is.null(countries) && is.null(shp)) {
      return(NULL)
    }

    if (!is.null(countries)) {
      selectedCountries <- if (is.list(countries)) {
        countries[[identifierLabel]]
      } else {
        countries
      }
      validRegions <- unique(ggplot2::map_data("world")$region)
      if (is.null(selectedCountries) ||
          any(!(selectedCountries %in% validRegions))) {
        warning("Some or all supplied country names are not recognised.")
        return(NULL)
      }
      return(ggplot2::map_data("world", regions = selectedCountries))
    }

    selectedShape <- if (is.list(shp)) shp[[identifierLabel]] else shp
    if (is.null(selectedShape)) {
      warning("No shp entry was found for identifier ", identifierLabel, ".")
      return(NULL)
    }
    ggplot2::fortify(selectedShape)
  }

  missingColour <- grDevices::rgb(
    255, 255, 255, max = 255, alpha = 125, names = "white50"
  )

  plots <- lapply(seq_along(identifierLabels), function(idIndex) {
    identifierLabel <- identifierLabels[[idIndex]]
    rasterResult <- rastersByIdentifier[[idIndex]]
    p <- rasterVis::gplot(rasterResult) + ggplot2::theme_linedraw()

    if (output == "density") {
      p$data$variable <- factor(
        as.character(p$data$variable),
        levels = periodLabels
      )
      legend <- if (logCount) "log10(record count)" else "Record count"
      p <- p +
        ggplot2::geom_tile(ggplot2::aes(fill = value)) +
        ggplot2::facet_wrap(~variable) +
        ggplot2::scale_fill_gradient(
          low = "red",
          high = "blue",
          na.value = missingColour,
          name = legend
        )
    } else {
      observedValues <- sort(unique(stats::na.omit(p$data$value)))
      colours <- grDevices::terrain.colors(
        max(1, length(observedValues)), rev = TRUE
      )
      if (length(observedValues) > 0) {
        names(colours) <- as.character(observedValues)
      }
      legend <- if (output == "nPeriods") {
        "Number of periods sampled"
      } else {
        paste0("Sampled in at least ", overlapThreshold, " periods")
      }
      p <- p +
        ggplot2::geom_tile(ggplot2::aes(fill = factor(value))) +
        ggplot2::scale_fill_manual(
          values = colours,
          na.value = missingColour,
          name = legend,
          na.translate = FALSE
        ) +
        ggplot2::ggtitle(identifierLabel)
      if (output == "overlap") {
        p <- p + ggplot2::theme(legend.position = "none")
      }
    }

    map <- boundaryForIdentifier(identifierLabel)
    if (!is.null(map)) {
      p <- p + ggplot2::geom_polygon(
        data = map,
        ggplot2::aes(x = long, y = lat, group = group),
        colour = "black",
        fill = missingColour,
        inherit.aes = FALSE
      )
    }
    p
  })
  names(plots) <- identifierLabels

  if (!returnRaster) {
    return(plots)
  }

  rasterLayers <- list()
  rasterNames <- character()
  for (idIndex in seq_along(rastersByIdentifier)) {
    identifierName <- make.names(identifierLabels[[idIndex]])
    idRasters <- rastersByIdentifier[[idIndex]]
    for (layer in seq_len(raster::nlayers(idRasters))) {
      rasterLayers[[length(rasterLayers) + 1]] <- idRasters[[layer]]
      if (output == "density") {
        rasterNames <- c(
          rasterNames,
          paste0("id.", identifierName, ".period.", periodLabels[[layer]])
        )
      } else {
        rasterNames <- c(rasterNames, paste0("id.", identifierName))
      }
    }
  }
  rasterNames <- make.unique(make.names(rasterNames))
  returnedRasters <- raster::stack(rasterLayers)
  names(returnedRasters) <- rasterNames

  list(plots = plots, rasters = returnedRasters)
}
