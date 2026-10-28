#' Assess reference environmental coverage with KDE highest-density regions
#'
#' Fits one standardised principal component analysis (PCA) to reference
#' environmental observations, projects occurrence environments into that
#' common space, and measures how much of the reference sample lies inside
#' sample-derived kernel-density highest-density regions (HDRs).
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
#' @param periods A list of non-overlapping numeric year vectors. Results retain
#'   the supplied period order.
#' @param envDat A numeric data frame with exactly one row per row of `dat`, in
#'   the same order, containing occurrence environmental measurements.
#' @param backgroundEnvDat A numeric data frame representing the reference
#'   environmental domain. Its names and column order must exactly match
#'   `envDat`.
#' @param xPC,yPC Distinct positive integers selecting the two PCA dimensions
#'   used for KDE and plotting. Defaults are PC1 and PC2.
#' @param maxSpatUncertainty Optional maximum permitted spatial uncertainty.
#'   Records with missing uncertainty or values above this threshold are
#'   removed before temporal and environmental assessment.
#' @param dateHandling Optional uncertain-date configuration. `NULL` preserves
#'   exact-year assignment; see [temporal-assignment].
#' @param hdrProbs Numeric vector of unique probability masses strictly between
#'   zero and one. Defaults to 50, 80, and 95 percent HDRs.
#' @param bandwidthAdjust A positive finite scalar multiplying bandwidth
#'   standard deviations. The KDE bandwidth matrix is consequently multiplied
#'   by `bandwidthAdjust^2`.
#' @param gridSize Integer from 25 to 501 giving the number of KDE evaluation
#'   points along each PCA axis. HDR thresholds use numerical integration over
#'   this grid.
#' @param bandwidth Optional fixed 2 by 2, symmetric, positive-definite KDE
#'   bandwidth matrix in the selected PCA coordinates. If `NULL`, one plug-in
#'   matrix is estimated with [ks::Hpi()] from all retained occurrence scores
#'   pooled across identifiers and periods. The same adjusted matrix is then
#'   used for every comparison.
#'
#' @details The PCA is fitted once to complete, finite rows of
#' `backgroundEnvDat` with `center = TRUE` and `scale. = TRUE`. This differs
#' intentionally from [assessEnvBias()], whose existing PCA defaults are not
#' changed. A common scaled reference PCA makes period and identifier scores
#' directly comparable when environmental variables use different units.
#'
#' For a sample KDE \eqn{\hat f}, the probability-mass HDR with requested mass
#' \eqn{p} is the level set \eqn{\{z: \hat f(z) \geq t_p\}}, where \eqn{t_p}
#' is chosen by accumulating gridded density from highest to lowest until mass
#' \eqn{p} is reached. This is not a contour of density divided by its maximum.
#' `backgroundCoverage` is the number of reference scores whose sample-KDE
#' density is at least \eqn{t_p}, divided by the number of reference scores.
#' Ties on the numerical grid can make the represented HDR mass slightly exceed
#' the requested value. Increasing `gridSize` reduces discretisation error at
#' additional computational cost.
#'
#' Grid limits are common to every comparison and include all reference and
#' retained sample scores plus five marginal bandwidth standard deviations.
#' The Gaussian KDE is evaluated with a ten-bandwidth-standard-deviation
#' computational support, making omitted kernel mass negligible, and without
#' stochastic jitter. Combinations
#' with fewer than five observations, rank-deficient sample scores, unavailable
#' bandwidths, or failed KDE estimation return `NA` coverage with an explicit
#' `status`. No alternative estimator is substituted.
#'
#' The background must represent the environmental domain of interest.
#' Coverage depends on its sampling design: for example, uniformly sampled
#' longitude-latitude raster cells need not represent equal land area.
#' Occurrence and background variables must be measured and processed
#' comparably. Repeated occurrences influence the KDE; records are not
#' automatically deduplicated because the appropriate sampling unit may be a
#' record, visit, site, grid cell, or site-year. Prepare any scientifically
#' appropriate aggregation before calling the function.
#'
#' KDE describes sampled observations in the selected two-dimensional
#' environmental projection. Coverage is not a probability of occurrence, a
#' causal estimate of selection bias, or an automatic correction. It depends on
#' sample size, bandwidth, PCA axes, background definition, and sampling unit.
#' Similar coverage does not rule out spatial, temporal, taxonomic, or
#' higher-dimensional environmental bias.
#'
#' Temporal and spatial-uncertainty filtering preserve exact row alignment
#' between `dat` and `envDat`. Every uncertain occurrence contributes to at most
#' one period. Temporal diagnostics are stored in the result's
#' `temporalDiagnostics` attribute.
#'
#' @return A list with `data`, `plot`, and `metadata`. `data` has one row per
#'   identifier, period, and HDR probability, with sample and background sizes,
#'   background count inside the HDR, `backgroundCoverage`, `hdrThreshold`,
#'   integrated `gridMass`, and `status`. `plot` shows the reference scores,
#'   retained sample scores, and the same HDR thresholds used numerically.
#'   `metadata` contains the reference PCA, selected components and explained
#'   variance, common bandwidth and source, grid bounds, retained row-aligned
#'   sample scores, and effective input sizes. Large density arrays are not
#'   returned.
#' @references Exploratory methodological provenance:
#'   <https://github.com/sacrevert/npmsEnvBiasChecks>.
#' @seealso [assessEnvBias()] and [temporal-assignment]
#' @export

assessEnvCoverage <- function(dat,
                              species,
                              x,
                              y,
                              year,
                              spatialUncertainty,
                              identifier,
                              periods,
                              envDat,
                              backgroundEnvDat,
                              xPC = 1,
                              yPC = 2,
                              maxSpatUncertainty = NULL,
                              dateHandling = NULL,
                              hdrProbs = c(0.50, 0.80, 0.95),
                              bandwidthAdjust = 1,
                              gridSize = 151,
                              bandwidth = NULL) {
  required <- c(species, x, y, spatialUncertainty, identifier, year)
  if (!is.data.frame(dat) || any(!(required %in% names(dat)))) {
    stop("You have specified columns that do not exist in dat.")
  }
  if (!is.data.frame(envDat) || !is.data.frame(backgroundEnvDat)) {
    stop("envDat and backgroundEnvDat must be data frames.")
  }
  if (nrow(envDat) != nrow(dat)) {
    stop("envDat must contain exactly one row per row of dat, in the same order.")
  }
  if (!identical(names(envDat), names(backgroundEnvDat))) {
    stop("envDat and backgroundEnvDat must have identical column names and order.")
  }
  if (ncol(envDat) < 2 || any(!vapply(envDat, is.numeric, logical(1))) ||
      any(!vapply(backgroundEnvDat, is.numeric, logical(1)))) {
    stop("Environmental data must contain at least two numeric variables.")
  }
  selectedPCs <- .envCoveragePCSelection(xPC, yPC)
  hdrProbs <- .envCoverageProbabilities(hdrProbs)
  if (!is.numeric(bandwidthAdjust) || length(bandwidthAdjust) != 1 ||
      is.na(bandwidthAdjust) || !is.finite(bandwidthAdjust) ||
      bandwidthAdjust <= 0) {
    stop("bandwidthAdjust must be a single positive, finite number.")
  }
  if (!is.numeric(gridSize) || length(gridSize) != 1 || is.na(gridSize) ||
      !is.finite(gridSize) || gridSize != floor(gridSize) ||
      gridSize < 25 || gridSize > 501) {
    stop("gridSize must be an integer from 25 to 501.")
  }
  if (!is.null(bandwidth)) .envCoverageValidateBandwidth(bandwidth)

  backgroundOriginal <- nrow(backgroundEnvDat)
  backgroundComplete <- .envCoverageFiniteRows(backgroundEnvDat)
  if (any(!backgroundComplete)) {
    warning(
      "Removing ", sum(!backgroundComplete),
      " background row(s) with missing or non-finite environmental values.",
      call. = FALSE
    )
    backgroundEnvDat <- backgroundEnvDat[backgroundComplete, , drop = FALSE]
  }
  if (nrow(backgroundEnvDat) < 3) {
    stop("At least three complete background environmental rows are required.")
  }
  backgroundSD <- vapply(backgroundEnvDat, stats::sd, numeric(1))
  if (any(!is.finite(backgroundSD)) || any(backgroundSD <= 0)) {
    stop("Every background environmental variable must have positive variance.")
  }

  pca <- stats::prcomp(
    backgroundEnvDat, center = TRUE, scale. = TRUE, retx = FALSE
  )
  estimable <- sum(
    pca$sdev > max(pca$sdev) * sqrt(.Machine$double.eps)
  )
  if (max(selectedPCs) > estimable) {
    stop("The selected PCA dimensions are not estimable from backgroundEnvDat.")
  }
  varianceExplained <- pca$sdev^2 / sum(pca$sdev^2)
  backgroundScoresAll <- stats::predict(pca, newdata = backgroundEnvDat)
  backgroundScores <- backgroundScoresAll[, selectedPCs, drop = FALSE]
  colnames(backgroundScores) <- c("xScore", "yScore")

  temporal <- assignTemporalPeriods(dat, periods, year, dateHandling)
  occurrence <- createData(
    data = dat, species, x, y, year, spatialUncertainty, identifier
  )
  occurrence$.row <- seq_len(nrow(dat))
  occurrence$period <- temporal$period
  occurrence$.temporalStatus <- temporal$status
  if (any(is.na(occurrence$identifier))) {
    stop("Missing values are not permitted in identifier.")
  }

  if (!is.null(maxSpatUncertainty)) {
    keep <- !is.na(occurrence$spatialUncertainty) &
      occurrence$spatialUncertainty <= maxSpatUncertainty
    occurrence <- occurrence[keep, , drop = FALSE]
    envDat <- envDat[keep, , drop = FALSE]
  }
  if (nrow(occurrence) == 0) {
    stop("No records remain after spatial-uncertainty filtering.")
  }

  occurrenceComplete <- .envCoverageFiniteRows(envDat)
  if (any(!occurrenceComplete)) {
    warning(
      "Removing ", sum(!occurrenceComplete),
      " occurrence row(s) with missing or non-finite environmental values.",
      call. = FALSE
    )
    occurrence <- occurrence[occurrenceComplete, , drop = FALSE]
    envDat <- envDat[occurrenceComplete, , drop = FALSE]
  }
  if (nrow(occurrence) == 0) {
    stop("No occurrence rows have complete finite environmental data.")
  }

  selection <- .finalizeTemporalAssignment(
    occurrence$period,
    occurrence$.temporalStatus,
    periods,
    dateHandling,
    exactMissingWarning = "Removing data without a specified year"
  )
  temporalDiagnostics <- selection$diagnostics
  occurrence <- occurrence[selection$keep, , drop = FALSE]
  envDat <- envDat[selection$keep, , drop = FALSE]
  if (nrow(occurrence) == 0) {
    stop("No records fall within the requested periods.")
  }

  sampleScoresAll <- stats::predict(pca, newdata = envDat)
  sampleScores <- sampleScoresAll[, selectedPCs, drop = FALSE]
  colnames(sampleScores) <- c("xScore", "yScore")
  sampleScoreData <- data.frame(
    row = occurrence$.row,
    identifier = occurrence$identifier,
    period = occurrence$period,
    xScore = sampleScores[, 1],
    yScore = sampleScores[, 2],
    stringsAsFactors = FALSE
  )

  bandwidthResult <- .envCoverageBandwidth(
    sampleScores, bandwidth, bandwidthAdjust
  )
  commonBandwidth <- bandwidthResult$bandwidth
  gridBounds <- .envCoverageGridBounds(
    backgroundScores, sampleScores, commonBandwidth
  )

  estimatedWork <- nrow(sampleScores) *
    (as.double(gridSize)^2 + nrow(backgroundScores))
  if (estimatedWork > 2e8) {
    stop(
      "The requested exact KDE evaluations exceed the computational safeguard; ",
      "reduce gridSize or aggregate the occurrence/background samples."
    )
  }

  identifiers <- unique(occurrence$identifier)
  combinations <- lapply(identifiers, function(currentIdentifier) {
    lapply(seq_along(periods), function(currentPeriod) {
      rows <- occurrence$identifier == currentIdentifier &
        occurrence$period == currentPeriod
      .envCoverageGroup(
        points = sampleScores[rows, , drop = FALSE],
        backgroundScores = backgroundScores,
        identifier = currentIdentifier,
        period = currentPeriod,
        hdrProbs = hdrProbs,
        bandwidth = commonBandwidth,
        bandwidthStatus = bandwidthResult$status,
        gridSize = as.integer(gridSize),
        gridBounds = gridBounds
      )
    })
  })
  combinations <- unlist(combinations, recursive = FALSE)
  coverage <- do.call(rbind, lapply(combinations, `[[`, "data"))
  rownames(coverage) <- NULL
  contourRows <- lapply(combinations, `[[`, "contours")
  contourRows <- contourRows[vapply(contourRows, nrow, integer(1)) > 0]
  contours <- if (length(contourRows) == 0) {
    .emptyEnvCoverageContours()
  } else {
    do.call(rbind, contourRows)
  }
  rownames(contours) <- NULL

  plot <- .envCoveragePlot(
    backgroundScores,
    sampleScoreData,
    contours,
    identifiers,
    length(periods),
    hdrProbs,
    selectedPCs,
    varianceExplained
  )

  metadata <- list(
    pca = pca,
    selectedPCs = selectedPCs,
    varianceExplained = varianceExplained,
    bandwidth = commonBandwidth,
    bandwidthSource = bandwidthResult$source,
    bandwidthStatus = bandwidthResult$status,
    bandwidthAdjust = bandwidthAdjust,
    gridSize = as.integer(gridSize),
    gridBounds = gridBounds,
    nBackgroundOriginal = backgroundOriginal,
    nBackground = nrow(backgroundScores),
    nSample = nrow(sampleScores),
    sampleScores = sampleScoreData
  )

  .addTemporalDiagnostics(
    list(data = coverage, plot = plot, metadata = metadata),
    temporalDiagnostics
  )
}

.envCoveragePCSelection <- function(xPC, yPC) {
  pcs <- c(xPC, yPC)
  if (!is.numeric(pcs) || length(pcs) != 2 || any(is.na(pcs)) ||
      any(!is.finite(pcs)) || any(pcs != floor(pcs)) || any(pcs < 1) ||
      pcs[[1]] == pcs[[2]]) {
    stop("xPC and yPC must be distinct positive integers.")
  }
  as.integer(pcs)
}

.envCoverageProbabilities <- function(probabilities) {
  if (!is.numeric(probabilities) || length(probabilities) == 0 ||
      any(is.na(probabilities)) || any(!is.finite(probabilities)) ||
      any(probabilities <= 0 | probabilities >= 1) ||
      anyDuplicated(probabilities)) {
    stop("hdrProbs must contain unique finite probabilities between zero and one.")
  }
  as.numeric(probabilities)
}

.envCoverageFiniteRows <- function(data) {
  if (nrow(data) == 0) return(logical())
  stats::complete.cases(data) &
    apply(data, 1, function(row) all(is.finite(row)))
}

.envCoverageValidateBandwidth <- function(bandwidth) {
  if (!is.matrix(bandwidth) || !identical(dim(bandwidth), c(2L, 2L)) ||
      !is.numeric(bandwidth) || any(!is.finite(bandwidth)) ||
      !isTRUE(all.equal(bandwidth, t(bandwidth), tolerance = 1e-10)) ||
      any(eigen(bandwidth, symmetric = TRUE, only.values = TRUE)$values <= 0)) {
    stop("bandwidth must be a finite symmetric positive-definite 2 by 2 matrix.")
  }
  invisible(bandwidth)
}

.envCoverageMatrixRank <- function(points) {
  centered <- scale(points, center = TRUE, scale = FALSE)
  singularValues <- svd(centered, nu = 0, nv = 0)$d
  if (length(singularValues) == 0 || max(singularValues) == 0) return(0L)
  as.integer(sum(
    singularValues > max(singularValues) * sqrt(.Machine$double.eps)
  ))
}

.envCoverageBandwidth <- function(points, bandwidth, adjust) {
  if (!is.null(bandwidth)) {
    adjusted <- bandwidth * adjust^2
    .envCoverageValidateBandwidth(adjusted)
    return(list(
      bandwidth = adjusted,
      source = "user",
      status = "ok"
    ))
  }
  if (nrow(points) < 5 || .envCoverageMatrixRank(points) < 2) {
    return(list(
      bandwidth = NULL,
      source = "pooled_Hpi",
      status = "pooled_sample_degenerate"
    ))
  }
  selected <- tryCatch(
    suppressWarnings(ks::Hpi(x = points)),
    error = function(error) NULL
  )
  if (is.null(selected)) {
    return(list(
      bandwidth = NULL,
      source = "pooled_Hpi",
      status = "bandwidth_selection_failed"
    ))
  }
  adjusted <- selected * adjust^2
  valid <- tryCatch({
    .envCoverageValidateBandwidth(adjusted)
    TRUE
  }, error = function(error) FALSE)
  if (!valid) {
    return(list(
      bandwidth = NULL,
      source = "pooled_Hpi",
      status = "bandwidth_selection_failed"
    ))
  }
  list(bandwidth = adjusted, source = "pooled_Hpi", status = "ok")
}

.envCoverageGridBounds <- function(background, samples, bandwidth) {
  combined <- rbind(background, samples)
  ranges <- apply(combined, 2, range)
  if (is.null(bandwidth)) {
    span <- pmax(ranges[2, ] - ranges[1, ], 1)
    margin <- 0.05 * span
  } else {
    margin <- 5 * sqrt(diag(bandwidth))
  }
  list(min = ranges[1, ] - margin, max = ranges[2, ] + margin)
}

.envCoverageHdrThresholds <- function(density, evalPoints, probabilities) {
  spacing <- vapply(evalPoints, function(axis) axis[[2]] - axis[[1]], numeric(1))
  cellArea <- prod(spacing)
  values <- as.vector(density)
  gridMass <- sum(values) * cellArea
  if (!is.finite(gridMass) || gridMass <= 0) {
    return(list(thresholds = rep(NA_real_, length(probabilities)),
                gridMass = gridMass))
  }
  orderHigh <- order(values, decreasing = TRUE)
  cumulative <- cumsum(values[orderHigh] * cellArea) / gridMass
  thresholds <- vapply(probabilities, function(probability) {
    values[orderHigh[[which(cumulative >= probability)[[1]]]]]
  }, numeric(1))
  list(thresholds = thresholds, gridMass = gridMass)
}

.emptyEnvCoverageContours <- function() {
  data.frame(
    xScore = numeric(), yScore = numeric(), hdrProbability = numeric(),
    hdrLabel = character(), contourGroup = character(),
    identifier = character(), period = integer(), stringsAsFactors = FALSE
  )
}

.envCoverageContourData <- function(fit, thresholds, probabilities,
                                    identifier, period) {
  pieces <- lapply(seq_along(probabilities), function(index) {
    lines <- suppressWarnings(grDevices::contourLines(
      fit$eval.points[[1]], fit$eval.points[[2]], fit$estimate,
      levels = thresholds[[index]]
    ))
    if (length(lines) == 0) return(NULL)
    rows <- lapply(seq_along(lines), function(piece) {
      data.frame(
        xScore = lines[[piece]]$x,
        yScore = lines[[piece]]$y,
        hdrProbability = probabilities[[index]],
        hdrLabel = paste0(format(100 * probabilities[[index]], trim = TRUE), "%"),
        contourGroup = paste(identifier, period, index, piece, sep = "."),
        identifier = as.character(identifier),
        period = period,
        stringsAsFactors = FALSE
      )
    })
    do.call(rbind, rows)
  })
  pieces <- pieces[!vapply(pieces, is.null, logical(1))]
  if (length(pieces) == 0) .emptyEnvCoverageContours() else do.call(rbind, pieces)
}

.envCoverageStatusRows <- function(identifier, period, probabilities,
                                   nSample, nBackground, status) {
  data.frame(
    identifier = rep(as.character(identifier), length(probabilities)),
    period = rep(period, length(probabilities)),
    hdrProbability = probabilities,
    nSample = rep(as.integer(nSample), length(probabilities)),
    nBackground = rep(as.integer(nBackground), length(probabilities)),
    nBackgroundInside = rep(NA_integer_, length(probabilities)),
    backgroundCoverage = rep(NA_real_, length(probabilities)),
    hdrThreshold = rep(NA_real_, length(probabilities)),
    gridMass = rep(NA_real_, length(probabilities)),
    status = rep(status, length(probabilities)),
    stringsAsFactors = FALSE
  )
}

.envCoverageGroup <- function(points,
                              backgroundScores,
                              identifier,
                              period,
                              hdrProbs,
                              bandwidth,
                              bandwidthStatus,
                              gridSize,
                              gridBounds) {
  nSample <- nrow(points)
  nBackground <- nrow(backgroundScores)
  invalid <- NULL
  if (nSample == 0) {
    invalid <- "no_observations"
  } else if (nSample < 5) {
    invalid <- "insufficient_observations"
  } else if (.envCoverageMatrixRank(points) < 2) {
    invalid <- "degenerate_sample"
  } else if (is.null(bandwidth)) {
    invalid <- bandwidthStatus
  }
  if (!is.null(invalid)) {
    return(list(
      data = .envCoverageStatusRows(
        identifier, period, hdrProbs, nSample, nBackground, invalid
      ),
      contours = .emptyEnvCoverageContours()
    ))
  }

  fit <- tryCatch(
    ks::kde(
      x = points,
      H = bandwidth,
      gridsize = rep(gridSize, 2),
      xmin = gridBounds$min,
      xmax = gridBounds$max,
      supp = 10,
      binned = FALSE,
      compute.cont = FALSE
    ),
    error = function(error) NULL
  )
  backgroundDensity <- tryCatch(
    ks::kde(
      x = points,
      H = bandwidth,
      eval.points = backgroundScores,
      supp = 10,
      binned = FALSE,
      compute.cont = FALSE
    )$estimate,
    error = function(error) NULL
  )
  if (is.null(fit) || is.null(backgroundDensity) ||
      any(!is.finite(fit$estimate)) || any(!is.finite(backgroundDensity))) {
    return(list(
      data = .envCoverageStatusRows(
        identifier, period, hdrProbs, nSample, nBackground, "kde_failed"
      ),
      contours = .emptyEnvCoverageContours()
    ))
  }

  hdr <- .envCoverageHdrThresholds(
    fit$estimate, fit$eval.points, hdrProbs
  )
  if (!is.finite(hdr$gridMass) || hdr$gridMass < 0.99 ||
      any(!is.finite(hdr$thresholds))) {
    rows <- .envCoverageStatusRows(
      identifier, period, hdrProbs, nSample, nBackground,
      "incomplete_grid_mass"
    )
    rows$gridMass <- hdr$gridMass
    return(list(data = rows, contours = .emptyEnvCoverageContours()))
  }

  inside <- vapply(hdr$thresholds, function(threshold) {
    sum(backgroundDensity >= threshold)
  }, integer(1))
  rows <- data.frame(
    identifier = rep(as.character(identifier), length(hdrProbs)),
    period = rep(period, length(hdrProbs)),
    hdrProbability = hdrProbs,
    nSample = rep(as.integer(nSample), length(hdrProbs)),
    nBackground = rep(as.integer(nBackground), length(hdrProbs)),
    nBackgroundInside = as.integer(inside),
    backgroundCoverage = inside / nBackground,
    hdrThreshold = hdr$thresholds,
    gridMass = rep(hdr$gridMass, length(hdrProbs)),
    status = rep("ok", length(hdrProbs)),
    stringsAsFactors = FALSE
  )
  contours <- .envCoverageContourData(
    fit, hdr$thresholds, hdrProbs, identifier, period
  )
  list(data = rows, contours = contours)
}

.envCoveragePlot <- function(backgroundScores,
                             sampleScores,
                             contours,
                             identifiers,
                             nPeriods,
                             hdrProbs,
                             selectedPCs,
                             varianceExplained) {
  background <- as.data.frame(backgroundScores)
  periods <- paste0("p", seq_len(nPeriods))
  sampleScores$period <- factor(
    paste0("p", sampleScores$period), levels = periods
  )
  if (nrow(contours) == 0) {
    contours$period <- factor(character(), levels = periods)
    contours$hdrLabel <- factor(
      character(),
      levels = paste0(format(100 * hdrProbs, trim = TRUE), "%")
    )
  } else {
    contours$period <- factor(paste0("p", contours$period), levels = periods)
    contours$hdrLabel <- factor(
      contours$hdrLabel,
      levels = paste0(format(100 * hdrProbs, trim = TRUE), "%")
    )
  }
  anchors <- expand.grid(
    identifier = as.character(identifiers),
    period = periods,
    stringsAsFactors = FALSE
  )
  anchors$period <- factor(anchors$period, levels = periods)
  anchors$xScore <- mean(background$xScore)
  anchors$yScore <- mean(background$yScore)

  ggplot2::ggplot() +
    ggplot2::geom_point(
      data = background,
      mapping = ggplot2::aes(x = xScore, y = yScore),
      colour = "grey55", alpha = 0.18, size = 0.45
    ) +
    ggplot2::geom_point(
      data = sampleScores,
      mapping = ggplot2::aes(x = xScore, y = yScore),
      colour = "#C44E52", alpha = 0.65, size = 0.75
    ) +
    ggplot2::geom_path(
      data = contours,
      mapping = ggplot2::aes(
        x = xScore, y = yScore, group = contourGroup, linetype = hdrLabel
      ),
      colour = "black", linewidth = 0.55
    ) +
    ggplot2::geom_blank(
      data = anchors,
      mapping = ggplot2::aes(x = xScore, y = yScore)
    ) +
    ggplot2::facet_grid(identifier ~ period, drop = FALSE) +
    ggplot2::coord_equal() +
    ggplot2::theme_linedraw() +
    ggplot2::labs(
      x = sprintf(
        "PC%d (%.1f%%)", selectedPCs[[1]],
        100 * varianceExplained[[selectedPCs[[1]]]]
      ),
      y = sprintf(
        "PC%d (%.1f%%)", selectedPCs[[2]],
        100 * varianceExplained[[selectedPCs[[2]]]]
      ),
      linetype = "Sample HDR probability mass"
    )
}
