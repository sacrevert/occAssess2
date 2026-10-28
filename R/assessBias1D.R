#' Compare one-dimensional sample and population distributions
#'
#' Compares binned relative-frequency distributions for a reference population,
#' one or more samples, and optionally weighted versions of those samples. When
#' requested, total variation distance and mean absolute error quantify each
#' sample-to-population comparison.
#'
#' @param pop A data frame with one row for every spatial unit in the reference
#'   population. It must contain the focal variable and one or more
#'   sample-inclusion columns.
#' @param R Character vector naming columns in `pop` that indicate sample
#'   inclusion. Values must be `0`, `1`, or `NA`; only `1` denotes inclusion.
#' @param x Character string naming the numeric variable whose distribution is
#'   assessed.
#' @param weights Optional weights. The documented and preferred format is one
#'   numeric vector of length `nrow(pop)`, or a list of such vectors, aligned to
#'   population rows. For compatibility, a list of weighting schemes may
#'   instead contain one list per scheme, with one sample-aligned numeric vector
#'   per element of `R`; each inner vector must have the length of that
#'   indicator's included rows before missing `x` values are removed. Mixing or
#'   ambiguously nesting these formats is an error.
#' @param breaks Either one positive integer giving the number of bins or a
#'   strictly increasing finite numeric vector of explicit bin boundaries.
#'   Explicit boundaries must cover every finite population value of `x`.
#' @param RNames Unique, non-missing character labels for the sample-inclusion
#'   columns, with one label per element of `R`.
#' @param WNames Unique, non-missing character labels for weighting schemes,
#'   with one label per scheme. Required when `weights` is supplied.
#' @param returnMetrics Logical. If `FALSE` (the default), preserve the
#'   established `list(plot = ..., data = ...)` result. If `TRUE`, append a
#'   `metrics` data frame.
#'
#' @details Bin boundaries are defined once from the finite population values
#' of `x`, using [cut()] with right-closed intervals and `include.lowest = TRUE`,
#' and then shared by every sample and weighting scheme. All bins are retained,
#' including bins with zero observations. Missing and non-finite `x` values are
#' excluded consistently from population, sample, and weighted denominators.
#' Missing sample-inclusion values are not sampled; other values besides zero
#' and one are rejected.
#'
#' For population proportions \eqn{q_k} and sample proportions \eqn{p_k} over
#' \eqn{K} common bins, total variation distance is
#' \deqn{TV = 0.5 \sum_k |p_k - q_k|}
#' and mean absolute error is
#' \deqn{MAE = K^{-1} \sum_k |p_k - q_k|.}
#' Consequently, \eqn{MAE = 2 TV / K}; these are two presentations of the same
#' binned discrepancies, not independent measures. TV ranges from zero for
#' identical distributions to one for disjoint binned support. Both measures
#' use proportions rather than raw spatial-unit counts.
#'
#' Weights are applied only after selecting the relevant sample. They must be
#' finite and non-negative; zero individual weights are allowed, but eligible
#' sampled observations must have positive total weight. Weighted bin totals
#' are divided by that total and never alter the reference-population
#' distribution.
#'
#' A smaller weighted discrepancy means improved agreement for this reference,
#' variable, and binning scheme only. Weighting can also increase discrepancy.
#' Neither outcome establishes representativeness for other variables,
#' geographical domains, or biological responses. TV measures a difference
#' between discrete binned distributions, not continuous densities. Changing
#' bin boundaries or bin count can change TV, while MAE additionally changes
#' mechanically through its \eqn{1/K} factor. The metrics are descriptive, not
#' formal tests or universal estimates of sampling bias.
#'
#' @return A list containing `plot`, the established ggplot relative-frequency
#'   display, and `data`, its long-format data frame with every bin represented.
#'   When `returnMetrics = TRUE`, the list also contains `metrics`, with one
#'   unweighted row per sampling indicator and one weighted row per indicator
#'   and weighting scheme. Its columns are `id`, `comparison`, `weightType`,
#'   `totalVariation`, `meanAbsoluteError`, `nBins`, `nPopulation`, `nSample`,
#'   and `status`. An empty or entirely unusable sample has `NA` discrepancies
#'   and status `no_eligible_sample`; the function never reports it as a zero
#'   discrepancy.
#'
#' @examples
#' population <- data.frame(
#'   environment = c(1, 1, 2, 2, 3, 3, 4, 4),
#'   sampled = c(1, 0, 1, 0, 0, 0, 1, 0)
#' )
#' result <- assessBias1D(
#'   population, R = "sampled", x = "environment", breaks = 4,
#'   RNames = "survey", returnMetrics = TRUE
#' )
#' result$metrics
#'
#' @export
assessBias1D <- function(pop,
                         R,
                         x,
                         weights = NULL,
                         breaks,
                         RNames,
                         WNames = NULL,
                         returnMetrics = FALSE) {
  if (!is.data.frame(pop)) stop("pop must be a data.frame")
  if (nrow(pop) == 0) stop("pop must contain at least one row")
  if (!is.character(R) || length(R) == 0 || anyNA(R) ||
      any(!nzchar(R)) || anyDuplicated(R)) {
    stop("R must contain unique, non-missing column names")
  }
  if (!is.character(x) || length(x) != 1 || is.na(x) || !nzchar(x)) {
    stop("x must be one column name")
  }
  if (any(!R %in% names(pop))) stop("R must be in colnames(pop)")
  if (!x %in% names(pop)) stop("x must be in colnames(pop)")
  if (!is.numeric(pop[[x]])) stop("x must name a numeric column")
  if (!is.character(RNames) || length(RNames) != length(R) ||
      anyNA(RNames) || any(!nzchar(RNames)) || anyDuplicated(RNames)) {
    stop("RNames must contain one unique, non-missing label per element of R")
  }
  if (!is.logical(returnMetrics) || length(returnMetrics) != 1 ||
      is.na(returnMetrics)) {
    stop("returnMetrics must be TRUE or FALSE")
  }

  sampleRows <- lapply(R, function(column) {
    indicator <- pop[[column]]
    invalid <- !is.na(indicator) & !(indicator %in% c(0, 1))
    if (any(invalid)) {
      stop("Sampling indicators must contain only 0, 1, or NA")
    }
    !is.na(indicator) & indicator == 1
  })

  binning <- .bias1DBins(pop[[x]], breaks)
  weightInfo <- .bias1DWeights(weights, WNames, nrow(pop), sampleRows, RNames)
  nPopulation <- sum(binning$valid)
  populationDistribution <- tabulate(
    binning$bin[binning$valid], nbins = binning$nBins
  ) / nPopulation

  dataRows <- list()
  metricRows <- list()
  dataIndex <- 0L
  metricIndex <- 0L

  for (indicatorIndex in seq_along(R)) {
    eligible <- sampleRows[[indicatorIndex]] & binning$valid
    nSample <- sum(eligible)
    sampleDistribution <- if (nSample == 0) {
      rep(NA_real_, binning$nBins)
    } else {
      tabulate(binning$bin[eligible], nbins = binning$nBins) / nSample
    }
    sampleStatus <- if (nSample == 0) "no_eligible_sample" else "ok"

    metricIndex <- metricIndex + 1L
    metricRows[[metricIndex]] <- .bias1DMetricRow(
      id = RNames[[indicatorIndex]], comparison = "sample",
      weightType = NA_character_, sampleDistribution = sampleDistribution,
      populationDistribution = populationDistribution, nBins = binning$nBins,
      nPopulation = nPopulation, nSample = nSample, status = sampleStatus
    )

    if (weightInfo$nSchemes == 0) {
      dataIndex <- dataIndex + 1L
      dataRows[[dataIndex]] <- .bias1DDataRows(
        sampleDistribution, populationDistribution,
        weightedDistribution = NULL, id = RNames[[indicatorIndex]],
        weightType = NA_character_, variableIndex = 1L
      )
    } else {
      for (schemeIndex in seq_len(weightInfo$nSchemes)) {
        sampleWeights <- .bias1DSampleWeights(
          weightInfo, schemeIndex, indicatorIndex,
          sampleRows[[indicatorIndex]], eligible
        )
        weightedDistribution <- if (nSample == 0) {
          rep(NA_real_, binning$nBins)
        } else {
          totalWeight <- sum(sampleWeights)
          if (!is.finite(totalWeight) || totalWeight <= 0) {
            stop(
              "Each weighting scheme must have positive total weight among ",
              "eligible sampled observations"
            )
          }
          weightedCounts <- rowsum(
            sampleWeights, binning$bin[eligible], reorder = FALSE
          )
          totals <- numeric(binning$nBins)
          totals[as.integer(rownames(weightedCounts))] <- weightedCounts[, 1]
          totals / totalWeight
        }

        dataIndex <- dataIndex + 1L
        dataRows[[dataIndex]] <- .bias1DDataRows(
          sampleDistribution, populationDistribution, weightedDistribution,
          id = RNames[[indicatorIndex]],
          weightType = weightInfo$names[[schemeIndex]],
          variableIndex = schemeIndex
        )
        metricIndex <- metricIndex + 1L
        metricRows[[metricIndex]] <- .bias1DMetricRow(
          id = RNames[[indicatorIndex]], comparison = "weighted_sample",
          weightType = weightInfo$names[[schemeIndex]],
          sampleDistribution = weightedDistribution,
          populationDistribution = populationDistribution,
          nBins = binning$nBins, nPopulation = nPopulation,
          nSample = nSample, status = sampleStatus
        )
      }
    }
  }

  dat <- do.call(rbind, dataRows)
  rownames(dat) <- NULL
  metrics <- do.call(rbind, metricRows)
  rownames(metrics) <- NULL

  p <- ggplot2::ggplot(
    data = dat,
    ggplot2::aes(y = value, x = bin, colour = variable)
  ) +
    ggplot2::geom_line() +
    ggplot2::theme_linedraw() +
    ggplot2::labs(colour = "", x = "", y = "Relative frequency")

  if (weightInfo$nSchemes > 0 || length(R) > 1) {
    p <- p + ggplot2::facet_grid(id ~ weightType, scales = "free_y")
  }

  result <- list(plot = p, data = dat)
  if (returnMetrics) result$metrics <- metrics
  result
}

.bias1DBins <- function(values, breaks) {
  valid <- !is.na(values) & is.finite(values)
  if (!any(valid)) stop("pop contains no finite values of x")
  usable <- values[valid]

  if (!is.numeric(breaks) || length(breaks) == 0 || anyNA(breaks) ||
      any(!is.finite(breaks))) {
    stop("breaks must be a positive integer or finite numeric boundaries")
  }
  if (length(breaks) == 1) {
    if (breaks != floor(breaks) || breaks <= 0) {
      stop("A scalar breaks value must be a positive integer")
    }
    nBins <- as.integer(breaks)
    bins <- if (nBins == 1L) {
      rep(1L, length(usable))
    } else {
      as.integer(cut(usable, breaks = nBins, include.lowest = TRUE))
    }
  } else {
    if (any(diff(breaks) <= 0)) {
      stop("Explicit breaks must be strictly increasing")
    }
    if (min(usable) < breaks[[1]] || max(usable) > breaks[[length(breaks)]]) {
      stop("Explicit breaks must cover every finite population value of x")
    }
    nBins <- length(breaks) - 1L
    bins <- as.integer(cut(
      usable, breaks = breaks, include.lowest = TRUE, right = TRUE
    ))
  }
  if (anyNA(bins)) stop("Unable to assign every finite population value to a bin")
  index <- rep(NA_integer_, length(values))
  index[valid] <- bins
  list(bin = index, valid = valid, nBins = as.integer(nBins))
}

.bias1DWeights <- function(weights, WNames, nPopulation, sampleRows, RNames) {
  if (is.null(weights)) {
    return(list(
      nSchemes = 0L, names = character(), format = "none", values = list()
    ))
  }
  if (!is.character(WNames) || length(WNames) == 0 || anyNA(WNames) ||
      any(!nzchar(WNames)) || anyDuplicated(WNames)) {
    stop("WNames must uniquely name every weighting scheme")
  }

  if (is.numeric(weights)) {
    schemes <- list(weights)
    format <- "population"
  } else if (is.list(weights) && length(weights) > 0 &&
             all(vapply(weights, is.numeric, logical(1)))) {
    schemes <- weights
    format <- "population"
  } else if (is.list(weights) && length(weights) > 0 &&
             all(vapply(weights, is.list, logical(1)))) {
    schemes <- weights
    format <- "nested_sample"
  } else {
    stop(
      "weights must be a population-row vector, a list of such vectors, or ",
      "a list of scheme-specific sample-vector lists"
    )
  }
  if (length(schemes) != length(WNames)) {
    stop("length of weighting schemes must equal length(WNames)")
  }

  if (format == "population") {
    validLengths <- vapply(schemes, length, integer(1)) == nPopulation
    if (!all(validLengths)) {
      stop("Population-row-aligned weight vectors must have nrow(pop) values")
    }
  } else {
    for (scheme in schemes) {
      if (length(scheme) != length(sampleRows) ||
          !all(vapply(scheme, is.numeric, logical(1)))) {
        stop("Each nested weighting scheme must contain one numeric vector per R")
      }
      lengths <- vapply(scheme, length, integer(1))
      expected <- vapply(sampleRows, sum, integer(1))
      if (!identical(lengths, expected)) {
        stop(
          "Nested sample-aligned weights must match the included-row count for ",
          paste(RNames, collapse = ", ")
        )
      }
    }
  }

  allWeights <- if (format == "population") {
    unlist(schemes, use.names = FALSE)
  } else {
    unlist(unlist(schemes, recursive = FALSE), use.names = FALSE)
  }
  if (anyNA(allWeights) || any(!is.finite(allWeights)) || any(allWeights < 0)) {
    stop("weights must be finite and non-negative")
  }
  list(
    nSchemes = length(schemes), names = WNames, format = format,
    values = schemes
  )
}

.bias1DSampleWeights <- function(weightInfo,
                                 schemeIndex,
                                 indicatorIndex,
                                 sampleRows,
                                 eligible) {
  if (weightInfo$format == "population") {
    return(weightInfo$values[[schemeIndex]][eligible])
  }
  aligned <- rep(NA_real_, length(sampleRows))
  aligned[sampleRows] <- weightInfo$values[[schemeIndex]][[indicatorIndex]]
  aligned[eligible]
}

.bias1DDataRows <- function(sampleDistribution,
                            populationDistribution,
                            weightedDistribution,
                            id,
                            weightType,
                            variableIndex) {
  wide <- data.frame(
    sample = sampleDistribution,
    population = populationDistribution,
    bin = seq_along(populationDistribution),
    var = variableIndex,
    id = id,
    weightType = weightType,
    stringsAsFactors = FALSE
  )
  if (!is.null(weightedDistribution)) {
    wide$weighted_sample <- weightedDistribution
    wide <- wide[, c(
      "weighted_sample", "sample", "population", "bin", "var", "id",
      "weightType"
    )]
  }
  reshape2::melt(
    wide,
    id = c("bin", "var", "id", "weightType"),
    variable.name = "variable",
    value.name = "value"
  )
}

.bias1DMetricRow <- function(id,
                             comparison,
                             weightType,
                             sampleDistribution,
                             populationDistribution,
                             nBins,
                             nPopulation,
                             nSample,
                             status) {
  if (status == "ok") {
    absoluteDifference <- abs(sampleDistribution - populationDistribution)
    totalVariation <- 0.5 * sum(absoluteDifference)
    meanAbsoluteError <- mean(absoluteDifference)
  } else {
    totalVariation <- NA_real_
    meanAbsoluteError <- NA_real_
  }
  data.frame(
    id = id,
    comparison = comparison,
    weightType = weightType,
    totalVariation = totalVariation,
    meanAbsoluteError = meanAbsoluteError,
    nBins = as.integer(nBins),
    nPopulation = as.integer(nPopulation),
    nSample = as.integer(nSample),
    status = status,
    stringsAsFactors = FALSE
  )
}
