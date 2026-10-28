#' Assess association between recorded range size and record count
#'
#' For each period and identifier, regresses species record count on the number
#' of occupied grid cells, used as a proxy for recorded range size. The summary
#' is either the coefficient of determination or Pearson correlation. No index
#' is calculated for combinations containing fewer than five species.
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
#' @param res Numeric. Spatial resolution at which to grid the data for estimation of range size.
#' @return A list with \code{plot}, a ggplot object, and \code{data}, a data
#'   frame containing \code{period}, \code{id}, and \code{index}.
#' @param prevPerPeriod Logical. If \code{TRUE}, occupied cells are calculated
#'   separately for each period. If \code{FALSE}, occupied cells for a species
#'   are calculated across all selected periods.
#' @param maxSpatUncertainty Numeric. Maximum permitted spatial uncertainty. All records more uncertain than this value will be dropped. Units must match the units in your data.
#' @param metric Character string. \code{"r2"} returns the coefficient of
#'   determination from the linear regression; \code{"cor"} returns Pearson's
#'   correlation coefficient.
#' @param dateHandling Optional uncertain-date configuration. `NULL` preserves
#'   exact-year assignment; see Details.
#' @details Temporal assignment is shared with the other period-based
#'   assessments. An uncertain record is allocated to at most one period before
#'   both record counts and occupied-cell calculations, including when
#'   `prevPerPeriod = FALSE`. Diagnostics are stored in the result's
#'   `temporalDiagnostics` attribute.
#' @seealso [temporal-assignment]
#' @export
#' 

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
                             dateHandling = NULL) {
  
  required <- c(species, x, y, spatialUncertainty, identifier, year)
  if (any(!(required %in% colnames(dat)))) stop("You have specified columns that don't exist in dat.")
  
  if (metric != "r2" & metric != "cor") stop("metric must be one of r2 or cor")

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

  if (any(is.na(dat$species))) {
    
    warning(paste("Removing", nrow(dat[is.na(dat$species), ]), "records because they are do not identified to species level."))
    
    dat <- dat[-which(is.na(dat$species)), ]
    
  }
  
  
  xmin <- min(dat$x, na.rm = T)
  
  xmax <- max(dat$x, na.rm = T)
  
  ymin <- min(dat$y, na.rm = T)
  
  ymax <- max(dat$y, na.rm = T)
  
  rast <- raster::raster(ncol=length(seq(xmin, xmax, res)),
                         nrow=length(seq(ymin, ymax, res)),
                         xmn=xmin,
                         xmx=xmax,
                         ymn=ymin,
                         ymx=ymax)

  for (i in unique(dat$identifier)) {

    x <- lapply(1:length(periods),
                function(y) {

                  spp <- unique(dat$species[dat$Period == y & dat$identifier == i])

                  if (length(spp) >= 5) {
                    
                    stats <- lapply(spp, 
                                    function(x) {
                                      if (prevPerPeriod == TRUE) {
                                        r <- raster::rasterize(cbind(dat$x[dat$species == x & dat$identifier == i & dat$Period == y], dat$y[dat$species == x& dat$identifier == i & dat$Period == y]), rast)
                                      } else {
                                        r <- raster::rasterize(cbind(dat$x[dat$species == x & dat$identifier == i], dat$y[dat$species == x& dat$identifier == i]), rast)
                                      }
                                      
                                      cells <- raster::getValues(r)
                                      cells <- length(cells[!is.na(cells)])
                                      recs <- nrow(dat[dat$species == x & dat$Period == y & dat$identifier == i, ])
                                      
                                      
                                      
                                      data.frame(species = x,
                                                 cells = cells,
                                                 recs = recs)
                                      
                                    })
                    
                    stats <- do.call("rbind", stats)

                    if(metric == "r2") {
                      
                      mod <- summary(stats::lm(stats$recs ~ stats$cells))$r.squared
                      
                    } else {
                      
                      mod <- stats::cor(stats$recs, stats$cells)
                      
                    }
                    
                    
                  } else {
                    
                    warning(paste("Fewer than five species in period", y, "for identifier", i,". No index will be calculated for this period/ identifier combination."))
                    
                    mod <- NA
                    
                  }
                  
                  
                  data.frame(period = y, 
                             id = i,
                             index = mod)
                })

    assign(paste0("out", i), do.call("rbind", x))
 
  }
  
  out <- lapply(unique(dat$identifier), 
                function(x) {get(paste0("out", x))})
  
  out <- do.call("rbind", out)
  
  p <- ggplot2::ggplot(data = out, ggplot2::aes(x = period, y = index, colour = id, group = id)) +
    ggplot2::theme_linedraw() +
    ggplot2::geom_point() +
    ggplot2::geom_line() + 
    ggplot2::xlab("Period") +
    ggplot2::ylab("Taxonomic bias index") +
    ggplot2::labs(colour = "")
  
  return(.addTemporalDiagnostics(
    list(plot = p, data = out), diagnostics
  ))
  
}

