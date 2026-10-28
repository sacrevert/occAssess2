#' Summarise records identified to species level
#'
#' Calculates the proportion or count of records with a non-missing species
#' value in each period and identifier group.
#' @param dat A data frame containing species, coordinate,
#'   spatial-uncertainty, identifier, and exact-year and/or interval-bound
#'   columns.
#' @param species Character string. column name in dat giving species names.
#' @param x string. Column name in dat giving x coordinates. Any coordinate and spatial reference systems are permitted.
#' @param y string. Column name in dat giving y coordinates. Any coordinate and spatial reference systems are permitted.
#' @param year Character string naming the exact-year column, or `NULL` in
#'   interval mode without an exact-year fallback.
#' @param spatialUncertainty String. Column name in dat giving uncertainty associated with x and y. Any units are permitted. 
#' @param identifier String. Column name in dat giving record "identifiers". Identifiers are used to group the data; heuristics will be calculated separately for each group.
#' @param periods String. A list of time periods. For example, for two periods, the first spanning 1950 to 1990, and the second 1991 to 2019: periods = list(1950:1990, 1991:2019).
#' @param type Character string: \code{"proportion"} or \code{"count"}.
#' @param maxSpatUncertainty Numeric. Maximum permitted spatial uncertainty. All records more uncertain than this value will be dropped. Units must match the units in your data.
#' @param dateHandling Optional uncertain-date configuration. `NULL` preserves
#'   exact-year assignment; see Details.
#' @details In interval mode, `dateHandling` names start/end columns and the
#'   allocation rules documented in the README. Diagnostics are stored in the
#'   result's `temporalDiagnostics` attribute.
#' @return A list with \code{data}, a data frame containing \code{prop},
#'   \code{group}, and \code{Period}; and \code{plot}, a ggplot object. The
#'   \code{prop} column contains proportions or counts according to \code{type}.
#' @seealso [temporal-assignment]
#' @export

assessSpeciesID <- function(dat, 
                            species,
                            x,
                            y,
                            year,
                            spatialUncertainty,
                            identifier,
                            periods,
                            type,
                            maxSpatUncertainty = NULL,
                            dateHandling = NULL) {

  required <- c(species, x, y, spatialUncertainty, identifier, year)
  if (any(!(required %in% colnames(dat)))) stop("You have specified columns that don't exist in dat.")

  temporal <- assignTemporalPeriods(dat, periods, year, dateHandling)
  
  dat <- createData(data = dat,
                    species,
                    x,
                    y,
                    year,
                    spatialUncertainty,
                    identifier)
  dat$Period <- temporal$period
  dat$.temporalStatus <- temporal$status
  
  if (any(is.na(dat$identifier))) stop("One or more NAs in the identifier field. NAs are not permitted.")
  
  if (!is.null(maxSpatUncertainty)) dat <- dat[!is.na(dat$spatialUncertainty) & dat$spatialUncertainty <= maxSpatUncertainty, ]
  
  if (nrow(dat) == 0) stop("No records with with spatialUncertainty < maxSpatUncertainty")
  
  selection <- .finalizeTemporalAssignment(
    dat$Period, dat$.temporalStatus, periods, dateHandling
  )
  diagnostics <- selection$diagnostics
  dat <- dat[selection$keep, , drop = FALSE]
  if (nrow(dat) == 0) stop("No records fall within the requested periods.")
  dat <- dat[order(dat$year), , drop = FALSE]

  groups <- unique(dat$identifier)

  data <- lapply(groups,
              function(i) {

                dat <- dat[dat$identifier == i, ]

                assign(paste0("props_", i), lapply(1:length(periods),
                                                   function(x) {

                                                     if (type == "proportion") {
                                                       
                                                       prop <- ifelse(length(dat$species[dat$Period == x]) > 0, 1 - (length(dat$species[is.na(dat$species) & dat$Period == x]) /
                                                                                                                       length(dat$species[dat$Period == x])), NA)
                                                       
                                                       data.frame(prop = prop,
                                                                  group = i,
                                                                  Period = x)
                                                       
                                                     } else {
                                                       
                                                       prop <- ifelse(length(dat$species[dat$Period == x]) > 0, length(dat$species[!is.na(dat$species) & dat$Period == x]), 0)
                                                       
                                                       data.frame(prop = prop,
                                                                  group = i,
                                                       
                                                                  Period = x)
                                                       
                                                     }
                                                     })
                )

                assign(paste0("props_", i), do.call("rbind", get(paste0("props_", i))))

              }
  )

  data <- do.call("rbind", data)

  #data <- data[order(data$year), ]
  
  ylab <- ifelse(type == "proportion", "Proportion identified to species level", "Number of records identified to species level")
  
  p <- ggplot2::ggplot(data = data, ggplot2::aes(y = prop, x = Period, colour = group, group = group)) +
    ggplot2::geom_point() +
    ggplot2::geom_line() +
    ggplot2::theme_linedraw() +
    ggplot2::ylab(ylab) +
    ggplot2::xlab("Period")
    ggplot2::labs(colour = "")

  return(.addTemporalDiagnostics(
    list(data = data, plot = p), diagnostics
  ))


}

