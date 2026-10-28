#' Assess reporting of spatial uncertainty
#'
#' Calculates, for each identifier, the proportion of retained records for which
#' spatial uncertainty is specified and plots the distributions of the specified
#' uncertainty values by period.
#'
#' @param dat A data frame containing columns named `species`, `x`, `y`, optionally `year`,
#'   `spatialUncertainty`, and `identifier`. `identifier` must not contain missing
#'   values. Missing `spatialUncertainty` values indicate that uncertainty was not
#'   specified. The units of non-missing uncertainty values must be consistent,
#'   but may be chosen by the user.
#' @param periods A list of vectors defining the years in each time period. For
#'   example, `list(1950:1990, 1991:2019)` defines two periods.
#' @param dateHandling Optional uncertain-date configuration. When non-`NULL`,
#'   start/end columns may replace `year`; see Details.
#'
#' @details Records without an assigned period are excluded using the shared
#' temporal rules documented in the README. For each identifier, the
#' returned value is pooled across all selected periods and is calculated as the
#' number of retained records with non-missing `spatialUncertainty` divided by
#' the total number of retained records. Thus, records with missing uncertainty
#' are included in the denominator but not the numerator.
#'
#' The histogram includes every non-missing spatial uncertainty value from the
#' retained records and overlays the distributions for the selected periods. No
#' fixed plotting range is applied. Its annotation is the pooled proportion, not
#' a period-specific proportion.
#'
#' @return A list with two elements: `data`, a one-column numeric matrix named
#'   `proportion_specified` with identifiers as row names and the supplied
#'   `periods` stored as an attribute; and `plot`, a ggplot object containing the
#'   period-grouped histograms and pooled annotations. Temporal diagnostics are
#'   stored in the result's `temporalDiagnostics` attribute.
#' @seealso [temporal-assignment]
#' @export


assessSpatialUncertainty <- function(dat, periods, dateHandling = NULL) {
  yearColumn <- if ("year" %in% names(dat)) "year" else NULL
  required <- c(
    "species", "x", "y", "spatialUncertainty", "identifier",
    if (is.null(dateHandling)) "year"
  )

  if (any(!(required %in% colnames(dat)))) {
    stop("Data must include species, x, y, spatialUncertainty, identifier, and an exact year unless interval mode is used")
  }

  temporal <- assignTemporalPeriods(dat, periods, yearColumn, dateHandling)
  dat$Period <- temporal$period
  dat$.temporalStatus <- temporal$status
  
  if (any(is.na(dat$identifier))) stop("One or more NAs in the identifier field. NAs are not permitted.")
  
  selection <- .finalizeTemporalAssignment(
    dat$Period, dat$.temporalStatus, periods, dateHandling,
    exactMissingWarning = function(n) sprintf(
      "Removing %d %s because %s not have a year associated.",
      n,
      if (n == 1) "record" else "records",
      if (n == 1) "it does" else "they do"
    )
  )
  diagnostics <- selection$diagnostics
  dat <- dat[selection$keep, , drop = FALSE]
  if (nrow(dat) == 0) stop("No records fall within the requested periods.")
  if (!is.null(yearColumn)) dat <- dat[order(dat[[yearColumn]]), , drop = FALSE]
  dat$Period <- paste0("p", dat$Period)

  identifiers <- unique(dat$identifier)
  proportions <- vapply(identifiers, function(id) {
    uncertainty <- dat$spatialUncertainty[dat$identifier == id]
    sum(!is.na(uncertainty)) / length(uncertainty)
  }, numeric(1))

  props <- matrix(
    unname(proportions),
    ncol = 1,
    dimnames = list(identifiers, "proportion_specified")
  )
  attr(props, "periods") <- periods

  text <- data.frame(
    identifier = identifiers,
    prop = unname(proportions)
  )

  p <- ggplot2::ggplot(
    data = dat,
    ggplot2::aes(
      x = spatialUncertainty,
      group = Period,
      colour = Period,
      fill = Period
    )
  ) +
    ggplot2::geom_histogram(position = "identity", alpha = 0.35) +
    ggplot2::theme_linedraw() +
    ggplot2::facet_wrap(~identifier) +
    ggplot2::xlab("Spatial uncertainty") +
    ggplot2::geom_text(
      data = text,
      mapping = ggplot2::aes(
        x = -Inf,
        y = Inf,
        label = sprintf("Proportion specified (pooled) = %.2f", prop)
      ),
      inherit.aes = FALSE,
      hjust = -0.1,
      vjust = 1.1,
      colour = "black"
    )

  return(.addTemporalDiagnostics(
    list(data = props, plot = p), diagnostics
  ))
}


