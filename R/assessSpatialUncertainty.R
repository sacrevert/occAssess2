#' Assess reporting of spatial uncertainty
#'
#' Calculates, for each identifier, the proportion of retained records for which
#' spatial uncertainty is specified and plots the distributions of the specified
#' uncertainty values by period.
#'
#' @param dat A data frame containing columns named `species`, `x`, `y`, `year`,
#'   `spatialUncertainty`, and `identifier`. `identifier` must not contain missing
#'   values. Missing `spatialUncertainty` values indicate that uncertainty was not
#'   specified. The units of non-missing uncertainty values must be consistent,
#'   but may be chosen by the user.
#' @param periods A list of vectors defining the years in each time period. For
#'   example, `list(1950:1990, 1991:2019)` defines two periods.
#'
#' @details Records with missing years are removed with a warning. Records whose
#' years do not occur in `periods` are also excluded. For each identifier, the
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
#'   period-grouped histograms and pooled annotations.
#' @export


assessSpatialUncertainty <- function(dat, periods) {
  required <- c("species", "x", "y", "year", "spatialUncertainty", "identifier")

  if (any(!(required %in% colnames(dat)))) {
    stop("Data must include columns for species, x, y, year, spatialUncertainty and identifier")
  }
  
  if (any(is.na(dat$identifier))) stop("One or more NAs in the identifier field. NAs are not permitted.")
  
  if (any(is.na(dat$year))) {
    nMissingYear <- sum(is.na(dat$year))
    warning(sprintf(
      "Removing %d %s because %s not have a year associated.",
      nMissingYear,
      if (nMissingYear == 1) "record" else "records",
      if (nMissingYear == 1) "it does" else "they do"
    ))

    dat <- dat[!is.na(dat$year), , drop = FALSE]
  }

  dat <- dat[order(dat$year), ]

  dat <- dat[dat$year %in% unlist(periods), , drop = FALSE]

  if (nrow(dat) == 0) stop("No records fall within the requested periods.")

  dat$Period <- NA_character_

  for (i in seq_along(periods)) {
    dat$Period <- ifelse(dat$year %in% periods[[i]], paste0("p", i), dat$Period)
  }

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

  return(list(data = props,
              plot = p))
}


