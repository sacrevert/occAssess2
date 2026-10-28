#' Assess association between recorded range size and record count
#'
#' For each period and identifier, regresses species record count on the number
#' of occupied grid cells, used as a proxy for recorded range size. The summary
#' is either the coefficient of determination or Pearson correlation. No index
#' is calculated for combinations containing fewer than five species.
#'
#' @param dat A data frame containing species, coordinate,
#'   spatial-uncertainty, identifier, and exact-year and/or interval-bound
#'   columns.
#' @param species Character string naming the species column in `dat`.
#' @param x Character string naming the x-coordinate column in `dat`.
#' @param y Character string naming the y-coordinate column in `dat`.
#' @param year Character string naming the exact-year column, or `NULL` in
#'   interval mode without an exact-year fallback.
#' @param spatialUncertainty Character string naming the spatial-uncertainty
#'   column in `dat`.
#' @param identifier Character string naming the grouping column in `dat`.
#' @param periods A list of non-overlapping numeric year vectors. Results use
#'   the supplied period order.
#' @param res A single positive, finite number giving the width and height of
#'   grid cells in the units of `x` and `y`.
#' @param prevPerPeriod Logical. If `TRUE`, occupied cells are calculated
#'   separately for each period. If `FALSE`, occupied cells for a species are
#'   calculated across all selected periods, while record counts remain
#'   period-specific.
#' @param maxSpatUncertainty Optional maximum permitted spatial uncertainty.
#'   Records with missing uncertainty or values above this threshold are
#'   dropped. Units must match the supplied data.
#' @param metric Character string. `"r2"` returns the coefficient of
#'   determination; `"cor"` returns Pearson correlation. Both diagnostics are
#'   calculated from the same untransformed OLS relationship.
#' @param dateHandling Optional uncertain-date configuration. `NULL` preserves
#'   exact-year assignment; see [temporal-assignment].
#' @param returnDiagnostics Logical. If `TRUE`, append `speciesDiagnostics` and
#'   `regressionDiagnostics` to the established result.
#' @param plotRegressions Logical. If `TRUE`, append `regressionPlots`, a named
#'   list with one period-faceted ggplot per identifier. This also makes both
#'   diagnostic data frames available.
#' @param labelTopResiduals Non-negative integer. In regression plots, label up
#'   to this many species with the largest finite absolute residual in each
#'   identifier-period panel. Zero, the default, draws no species labels.
#'
#' @details For every identifier-period combination, the fitted model is
#'
#' `recordCount = intercept + slope * occupiedCells + residual`.
#'
#' Counts and models use their original, untransformed scales. The raw residual
#' is `recordCount - fittedCount`; a positive value means that a species has
#' more records than predicted conditional on its estimated occupied-cell
#' count. It is not proof of over-recording, a significance classification, or
#' a direct measurement of sampling bias.
#'
#' Grid cells are aligned to multiples of `res` relative to coordinate zero,
#' using the same convention as [assessSpatialCov()]. Boundaries are therefore
#' deterministic and cells have exactly the requested width and height. Records
#' with missing or non-finite coordinates remain in `recordCount` but cannot
#' contribute an occupied cell. Missing species identities are excluded.
#'
#' At least five represented species and variation in `occupiedCells` are
#' required for a valid two-parameter regression. Constant record counts permit
#' fitted counts but make R-squared and correlation undefined. Non-finite
#' standardised residuals, leverage values, or Cook's distances are returned as
#' `NA`, rather than as misleading infinite values.
#'
#' `residualRMSE` is `sqrt(mean(residual^2))` across the species fitted in that
#' identifier-period combination. Occupied-cell count is only a proxy for
#' commonness and is estimated from the same records as the response. A strong
#' linear relationship therefore does not establish unbiased sampling, and
#' highly recorded or influential species can strongly affect the fitted line.
#'
#' Temporal assignment is shared with the other period-based assessments. Each
#' uncertain record is allocated to at most one period before record and cell
#' counts are calculated. Temporal diagnostics remain in the result's
#' `temporalDiagnostics` attribute.
#'
#' @return By default, the established list with `plot`, the index time series,
#'   and `data`, a data frame with `period`, `id`, and `index`. With diagnostics
#'   requested, `speciesDiagnostics` contains `identifier`, `period`, `species`,
#'   `occupiedCells`, `recordCount`, `fittedCount`, `residual`,
#'   `standardizedResidual`, `leverage`, and `cooksDistance`;
#'   `regressionDiagnostics` contains one row per identifier-period combination
#'   with `nSpecies`, coefficients, both indices, residual RMSE, and `status`.
#'   `regressionPlots` is added only when `plotRegressions = TRUE`.
#' @references Boyd et al. (2021), *Ecology and Evolution*.
#'   \doi{10.1002/ece3.8299}
#' @seealso [assessSpatialCov()] and [temporal-assignment]
#' @export

assessRarityBias <- function(dat,
                             species,
                             x,
                             y,
                             year,
                             spatialUncertainty,
                             identifier,
                             periods,
                             res,
                             prevPerPeriod,
                             maxSpatUncertainty = NULL,
                             metric = "r2",
                             dateHandling = NULL,
                             returnDiagnostics = FALSE,
                             plotRegressions = FALSE,
                             labelTopResiduals = 0) {
  required <- c(species, x, y, spatialUncertainty, identifier, year)
  if (any(!(required %in% colnames(dat)))) {
    stop("You have specified columns that don't exist in dat.")
  }
  if (!(is.character(metric) && length(metric) == 1 &&
        !is.na(metric) && metric %in% c("r2", "cor"))) {
    stop("metric must be one of r2 or cor")
  }
  if (!is.numeric(res) || length(res) != 1 || is.na(res) ||
      !is.finite(res) || res <= 0) {
    stop("res must be a single positive, finite number.")
  }
  .rarityLogicalArgument(prevPerPeriod, "prevPerPeriod")
  .rarityLogicalArgument(returnDiagnostics, "returnDiagnostics")
  .rarityLogicalArgument(plotRegressions, "plotRegressions")
  if (!is.numeric(labelTopResiduals) || length(labelTopResiduals) != 1 ||
      is.na(labelTopResiduals) || !is.finite(labelTopResiduals) ||
      labelTopResiduals < 0 || labelTopResiduals != floor(labelTopResiduals)) {
    stop("labelTopResiduals must be a single non-negative integer.")
  }

  temporal <- assignTemporalPeriods(dat, periods, year, dateHandling)
  dat <- createData(
    data = dat, species, x, y, year, spatialUncertainty, identifier
  )
  if (!is.numeric(dat$x) || !is.numeric(dat$y)) {
    stop("The x and y coordinate columns must be numeric.")
  }
  dat$Period <- temporal$period
  dat$.temporalStatus <- temporal$status

  if (any(is.na(dat$identifier))) {
    stop("One or more NAs in the identifier field. NAs are not permitted.")
  }
  if (!is.null(maxSpatUncertainty)) {
    dat <- dat[
      !is.na(dat$spatialUncertainty) &
        dat$spatialUncertainty <= maxSpatUncertainty,
      , drop = FALSE
    ]
  }
  if (nrow(dat) == 0) {
    stop("No records remain after spatial-uncertainty filtering.")
  }

  selection <- .finalizeTemporalAssignment(
    dat$Period, dat$.temporalStatus, periods, dateHandling
  )
  temporalDiagnostics <- selection$diagnostics
  dat <- dat[selection$keep, , drop = FALSE]
  if (nrow(dat) == 0) stop("No records fall within the requested periods.")
  identifiers <- unique(dat$identifier)

  missingSpecies <- is.na(dat$species)
  if (any(missingSpecies)) {
    warning(
      "Removing ", sum(missingSpecies),
      " record(s) without a species identity.",
      call. = FALSE
    )
    dat <- dat[!missingSpecies, , drop = FALSE]
  }

  invalidCoordinates <- !is.finite(dat$x) | !is.finite(dat$y)
  if (any(invalidCoordinates)) {
    warning(
      sum(invalidCoordinates),
      " record(s) with missing or non-finite coordinates cannot contribute ",
      "to occupied-cell counts but remain in record counts.",
      call. = FALSE
    )
  }
  dat$.cellX <- rep(NA_real_, nrow(dat))
  dat$.cellY <- rep(NA_real_, nrow(dat))
  if (any(!invalidCoordinates)) {
    bins <- .occAssessCellBins(
      dat$x[!invalidCoordinates], dat$y[!invalidCoordinates], res
    )
    dat$.cellX[!invalidCoordinates] <- bins$x
    dat$.cellY[!invalidCoordinates] <- bins$y
  }

  combinations <- lapply(identifiers, function(id) {
    lapply(seq_along(periods), function(period) {
      result <- .rarityBiasCombination(dat, id, period, prevPerPeriod)
      if (result$regression$status == "insufficient_species") {
        warning(
          "Fewer than five species in period ", period,
          " for identifier ", id,
          ". No index will be calculated for this period/identifier combination.",
          call. = FALSE
        )
      }
      result
    })
  })
  combinations <- unlist(combinations, recursive = FALSE)

  regressionDiagnostics <- do.call(
    rbind, lapply(combinations, `[[`, "regression")
  )
  speciesRows <- lapply(combinations, `[[`, "species")
  speciesRows <- speciesRows[vapply(speciesRows, nrow, integer(1)) > 0]
  speciesDiagnostics <- if (length(speciesRows) == 0) {
    .emptyRaritySpeciesDiagnostics()
  } else {
    do.call(rbind, speciesRows)
  }
  rownames(regressionDiagnostics) <- NULL
  rownames(speciesDiagnostics) <- NULL

  out <- data.frame(
    period = regressionDiagnostics$period,
    id = rep(identifiers, each = length(periods)),
    index = if (metric == "r2") {
      regressionDiagnostics$rSquared
    } else {
      regressionDiagnostics$correlation
    },
    stringsAsFactors = FALSE
  )

  p <- ggplot2::ggplot(
    data = out,
    ggplot2::aes(x = period, y = index, colour = id, group = id)
  ) +
    ggplot2::theme_linedraw() +
    ggplot2::geom_point() +
    ggplot2::geom_line() +
    ggplot2::xlab("Period") +
    ggplot2::ylab("Taxonomic bias index") +
    ggplot2::labs(colour = "")

  result <- list(plot = p, data = out)
  if (returnDiagnostics || plotRegressions) {
    result$speciesDiagnostics <- speciesDiagnostics
    result$regressionDiagnostics <- regressionDiagnostics
  }
  if (plotRegressions) {
    result$regressionPlots <- .rarityRegressionPlots(
      speciesDiagnostics,
      regressionDiagnostics,
      identifiers,
      length(periods),
      as.integer(labelTopResiduals)
    )
  }

  .addTemporalDiagnostics(result, temporalDiagnostics)
}

.rarityLogicalArgument <- function(value, name) {
  if (!is.logical(value) || length(value) != 1 || is.na(value)) {
    stop(name, " must be TRUE or FALSE.")
  }
  invisible(value)
}

.emptyRaritySpeciesDiagnostics <- function() {
  data.frame(
    identifier = character(), period = integer(), species = character(),
    occupiedCells = integer(), recordCount = integer(), fittedCount = numeric(),
    residual = numeric(), standardizedResidual = numeric(), leverage = numeric(),
    cooksDistance = numeric(), stringsAsFactors = FALSE
  )
}

.rarityFiniteOrNA <- function(values) {
  values[!is.finite(values)] <- NA_real_
  as.numeric(values)
}

.rarityBiasCombination <- function(dat, id, period, prevPerPeriod) {
  inPeriod <- dat$identifier == id & dat$Period == period
  species <- unique(dat$species[inPeriod])

  if (length(species) == 0) {
    regression <- data.frame(
      identifier = as.character(id), period = period, nSpecies = 0L,
      intercept = NA_real_, slope = NA_real_, rSquared = NA_real_,
      correlation = NA_real_, residualRMSE = NA_real_, status = "no_species",
      stringsAsFactors = FALSE
    )
    return(list(species = .emptyRaritySpeciesDiagnostics(),
                regression = regression))
  }

  stats <- do.call(rbind, lapply(species, function(currentSpecies) {
    speciesRows <- dat$identifier == id & dat$species == currentSpecies
    rangeRows <- speciesRows
    if (prevPerPeriod) rangeRows <- rangeRows & dat$Period == period
    validCells <- rangeRows & !is.na(dat$.cellX) & !is.na(dat$.cellY)
    occupied <- if (any(validCells)) {
      nrow(unique(dat[validCells, c(".cellX", ".cellY"), drop = FALSE]))
    } else {
      0L
    }
    data.frame(
      identifier = as.character(id),
      period = period,
      species = as.character(currentSpecies),
      occupiedCells = as.integer(occupied),
      recordCount = as.integer(sum(inPeriod & dat$species == currentSpecies)),
      fittedCount = NA_real_, residual = NA_real_,
      standardizedResidual = NA_real_, leverage = NA_real_,
      cooksDistance = NA_real_, stringsAsFactors = FALSE
    )
  }))

  nSpecies <- nrow(stats)
  status <- "ok"
  model <- NULL
  if (nSpecies < 5) {
    status <- "insufficient_species"
  } else if (length(unique(stats$occupiedCells)) < 2) {
    status <- "constant_occupied_cells"
  } else {
    model <- stats::lm(recordCount ~ occupiedCells, data = stats)
    if (length(unique(stats$recordCount)) < 2) status <- "constant_record_count"
  }

  intercept <- slope <- rSquared <- correlation <- residualRMSE <- NA_real_
  if (!is.null(model)) {
    coefficients <- stats::coef(model)
    intercept <- unname(coefficients[["(Intercept)"]])
    slope <- unname(coefficients[["occupiedCells"]])
    stats$fittedCount <- .rarityFiniteOrNA(stats::fitted(model))
    stats$residual <- .rarityFiniteOrNA(stats::residuals(model))
    stats$standardizedResidual <- .rarityFiniteOrNA(
      suppressWarnings(stats::rstandard(model))
    )
    stats$leverage <- .rarityFiniteOrNA(stats::hatvalues(model))
    stats$cooksDistance <- .rarityFiniteOrNA(
      suppressWarnings(stats::cooks.distance(model))
    )
    residualRMSE <- sqrt(mean(stats$residual^2))
    if (status == "constant_record_count") {
      stats$standardizedResidual <- NA_real_
      stats$cooksDistance <- NA_real_
    } else {
      rSquared <- .rarityFiniteOrNA(
        suppressWarnings(summary(model)$r.squared)
      )
      correlation <- .rarityFiniteOrNA(
        stats::cor(stats$recordCount, stats$occupiedCells)
      )
    }
  }

  regression <- data.frame(
    identifier = as.character(id), period = period,
    nSpecies = as.integer(nSpecies), intercept = intercept, slope = slope,
    rSquared = rSquared, correlation = correlation,
    residualRMSE = residualRMSE, status = status,
    stringsAsFactors = FALSE
  )
  list(species = stats, regression = regression)
}

.rarityRegressionPlots <- function(speciesDiagnostics,
                                   regressionDiagnostics,
                                   identifiers,
                                   nPeriods,
                                   labelTopResiduals) {
  periodLevels <- seq_len(nPeriods)
  plots <- lapply(identifiers, function(id) {
    points <- speciesDiagnostics[speciesDiagnostics$identifier == id, , drop = FALSE]
    summaries <- regressionDiagnostics[
      regressionDiagnostics$identifier == id, , drop = FALSE
    ]
    points$period <- factor(points$period, levels = periodLevels)
    summaries$period <- factor(summaries$period, levels = periodLevels)

    lines <- summaries[
      is.finite(summaries$intercept) & is.finite(summaries$slope),
      , drop = FALSE
    ]
    summaries$annotation <- paste0(
      "n = ", summaries$nSpecies, "\nR-squared = ",
      ifelse(is.finite(summaries$rSquared),
             format(round(summaries$rSquared, 3), nsmall = 3), "NA")
    )

    plot <- ggplot2::ggplot(
      points, ggplot2::aes(x = occupiedCells, y = recordCount)
    ) +
      ggplot2::geom_point() +
      ggplot2::geom_abline(
        data = lines,
        mapping = ggplot2::aes(intercept = intercept, slope = slope),
        inherit.aes = FALSE,
        colour = "#3366AA"
      ) +
      ggplot2::geom_text(
        data = summaries,
        mapping = ggplot2::aes(x = -Inf, y = Inf, label = annotation),
        inherit.aes = FALSE,
        hjust = -0.05,
        vjust = 1.1
      ) +
      ggplot2::facet_wrap(~period, drop = FALSE) +
      ggplot2::theme_linedraw() +
      ggplot2::labs(
        x = "Occupied grid cells",
        y = "Occurrence records",
        title = as.character(id)
      )

    if (labelTopResiduals > 0 && nrow(points) > 0) {
      candidates <- points[is.finite(points$residual), , drop = FALSE]
      if (nrow(candidates) > 0) {
        groups <- split(candidates, candidates$period, drop = TRUE)
        labels <- do.call(rbind, lapply(groups, function(group) {
          group[
            utils::head(order(abs(group$residual), decreasing = TRUE),
                        labelTopResiduals),
            , drop = FALSE
          ]
        }))
        plot <- plot + ggplot2::geom_text(
          data = labels,
          mapping = ggplot2::aes(label = species),
          vjust = -0.6,
          check_overlap = TRUE
        )
      }
    }
    plot
  })
  names(plots) <- as.character(identifiers)
  plots
}
