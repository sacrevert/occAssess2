#' Count recorded species by period
#'
#' Calculates the number of unique, non-missing species recorded in each
#' user-defined period and identifier group.
#' @param dat A data frame containing species, coordinate,
#'   spatial-uncertainty, identifier, and exact-year and/or interval-bound
#'   columns.
#' @param species Character string. column name in dat giving species names.
#' @param x string. Column name in dat giving x coordinates. Any coordinate and spatial reference systems are permitted.
#' @param y string. Column name in dat giving y coordinates. Any coordinate and spatial reference systems are permitted.
#' @param year Character string naming the exact-year column, or `NULL` in
#'   interval mode when no exact-year fallback is available.
#' @param spatialUncertainty String. Column name in dat giving uncertainty associated with x and y. Any units are permitted. 
#' @param identifier String. Column name in dat giving record "identifiers". Identifiers are used to group the data; heuristics will be calculated separately for each group.
#' @param periods String. A list of time periods. For example, for two periods, the first spanning 1950 to 1990, and the second 1991 to 2019: periods = list(1950:1990, 1991:2019).
#' @param maxSpatUncertainty Numeric. Maximum permitted spatial uncertainty. All records more uncertain than this value will be dropped. Units must match the units in your data.
#' @param normalize Logical. Whether or not to normalize the counts for each level of identifier. This is helpful where the number of records varies widely among levels of identifier.
#' @param dateHandling Optional uncertain-date configuration; see Details.
#' @details `dateHandling = NULL` preserves exact-year assignment. In interval
#'   mode, name start/end columns and allocation rules as documented in the
#'   README. Diagnostics are stored in the result's `temporalDiagnostics`
#'   attribute.
#' @seealso [assessSpeciesID()] and [temporal-assignment]
#' @return A list with \code{data}, a data frame containing \code{val},
#'   \code{group}, and \code{Period}; and \code{plot}, a ggplot object.
#' @export

assessSpeciesNumber <- function(dat, 
                                species,
                                x,
                                y,
                                year,
                                spatialUncertainty,
                                identifier,
                                periods,
                                maxSpatUncertainty = NULL,
                                normalize = FALSE,
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
                                                        
                                                        val <- ifelse(length(dat$species[dat$Period == x & !is.na(dat$species)]) > 0, length(unique(dat$species[dat$Period == x & !is.na(dat$species)])), 0)
                                                        
                                                        data.frame(val = val,
                                                                   group = i,
                                                                   Period = x)
                                                        
                                                      })
                   )
                   
                   assign(paste0("props_", i), do.call("rbind", get(paste0("props_", i))))
                   
                 }
  )
  
  data <- do.call("rbind", data)
  
  if (normalize) {
    
    for (i in groups) {
      
      data$val[data$group == i] <- data$val[data$group == i] / max(data$val[data$group == i])
      
    }
    
  }

  #data <- data[order(data$Period), ]

  
  p <- ggplot2::ggplot(data = data, ggplot2::aes(y = val, x = Period, colour = group, group = group)) +
    ggplot2::geom_point() +
    ggplot2::geom_line() +
    ggplot2::theme_linedraw() +
    ggplot2::ylab("Number of species recorded") +
    ggplot2::xlab("Period") +
    ggplot2::labs(colour = "")

  return(.addTemporalDiagnostics(
    list(data = data, plot = p), diagnostics
  ))
  
  
}
