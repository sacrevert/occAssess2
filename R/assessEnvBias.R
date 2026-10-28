#' Plot occurrences in principal-component environmental space
#'
#' Uses principal component analysis (PCA) to represent occurrence records in
#' environmental space and plots selected component pairs by period and
#' identifier.
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
#' @param envDat A data frame of numeric environmental variables with exactly
#'   one row per row of \code{dat}, in the same order. Filtering is applied to
#'   both objects so retained occurrence and environmental rows remain paired.
#' @param backgroundEnvDat Optional data frame containing a background sample of
#'   environmental space. Its column names and order must match \code{envDat};
#'   its number of rows need not match \code{dat}.
#' @param xPC Numeric. Which principal component to use as the x axis.
#' @param yPC Numeric. As xPC but for the y axis. 
#' @param maxSpatUncertainty Numeric. Maximum permitted spatial uncertainty. All records more uncertain than this value will be dropped. Units must match the units in your data.
#' @param dateHandling Optional uncertain-date configuration. `NULL` preserves
#'   exact-year assignment; see Details.
#' @details Records without an assigned period, records outside \code{periods},
#'   and rows with incomplete occurrence or environmental data are omitted. If
#'   \code{backgroundEnvDat} is supplied, one PCA is fitted to the background
#'   data and occurrence rows are projected onto it. Otherwise, a separate PCA
#'   is fitted to the retained occurrence environmental data for each identifier.
#'   Background scores are repeated in the returned data for each identifier so
#'   each plot facet has its own background reference.
#'
#'   Temporal assignment is performed row-for-row before occurrence and
#'   environmental rows are filtered together. Each uncertain occurrence is
#'   allocated to at most one period, and temporal diagnostics are stored in the
#'   result's `temporalDiagnostics` attribute.
#' @return A list with \code{data}, a data frame containing period, identifier,
#'   PCA scores, and percentages of variance explained by the selected axes; and
#'   \code{plot}, a ggplot object.
#' @seealso [temporal-assignment]
#' @export

assessEnvBias <- function(dat,
                          species,
                          x,
                          y,
                          spatialUncertainty, 
                          year, 
                          identifier,
                          periods,
                          envDat,
                          backgroundEnvDat = NULL,
                          xPC = 1,
                          yPC = 2, 
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
  
  if (nrow(envDat) != nrow(dat)) stop("nrow of environmental data does not equal nrow of species occurrence data.")
  
  if (xPC > ncol(envDat) | yPC > ncol(envDat)) stop("You have chosen a principal component that doesn't exist for one of the x or y axes")
  
  if (!is.null(maxSpatUncertainty)) {
    keep <- !is.na(dat$spatialUncertainty) & dat$spatialUncertainty <= maxSpatUncertainty
    dat <- dat[keep, , drop = FALSE]
    envDat <- envDat[keep, , drop = FALSE]
  }
  
  if (nrow(dat) == 0) stop("No records with with spatialUncertainty < maxSpatUncertainty")
  
  if (!is.null(backgroundEnvDat) & any(colnames(envDat) != colnames(backgroundEnvDat))) stop("Column names of envDat must match column names of backgrounEnvDat")
  
  dat <- cbind(dat, envDat)
  
  envCols <- ((ncol(dat) - ncol(envDat)) + 1):ncol(dat)
  
  completeEnv <- stats::complete.cases(dat[, envCols, drop = FALSE])
  if (any(!completeEnv)) dat <- dat[completeEnv, , drop = FALSE]

  selection <- .finalizeTemporalAssignment(
    dat$Period, dat$.temporalStatus, periods, dateHandling,
    exactMissingWarning = "Removing data without a specified year"
  )
  diagnostics <- selection$diagnostics
  dat <- dat[selection$keep, , drop = FALSE]
  if (nrow(dat) == 0) stop("No records fall within the requested periods.")
  dat <- dat[order(dat$year), , drop = FALSE]
  dat$Period <- paste0("p", dat$Period)
  
  if (is.null(backgroundEnvDat)) {
    
    plotDat <- lapply(unique(dat$identifier),
                      function(x) { 
                        pca <- stats::prcomp(dat[dat$identifier == x, envCols])
                        scores <- pca$x
                        data.frame(Period = dat$Period[dat$identifier == x],
                                   identifier = x,
                                   scores = scores,
                                   xVar = summary(pca)$importance[2,xPC] * 100,
                                   yVar = summary(pca)$importance[2,yPC] * 100)})
    
    plotDat <- do.call("rbind", plotDat)
    
  } else {
    
    plotDat <- lapply(unique(dat$identifier),
                      function(x) { 
                        pca <- stats::prcomp(backgroundEnvDat)
                        pca2 <- stats::predict(pca, dat[dat$identifier == x, envCols])
                        scores <- rbind(pca2, pca$x)
                        data.frame(Period = c(dat$Period[dat$identifier == x], rep("background", nrow(backgroundEnvDat))),
                                   identifier = x,
                                   scores = scores,
                                   xVar = summary(pca)$importance[2,xPC] * 100,
                                   yVar = summary(pca)$importance[2,yPC] * 100)})
    
    plotDat <- do.call("rbind", plotDat)
    
  }
  
  if (!is.null(backgroundEnvDat)) {
    
    p <- ggplot2::ggplot(data = plotDat, ggplot2::aes(x = plotDat[, (2 + xPC)], y = plotDat[, (2+yPC)], 
                                                      colour = Period, group = Period,
                                                      fill = factor(ifelse(Period == "background", "background", "samples")))) + 
      ggplot2::stat_ellipse(type = "norm", geom = "polygon", alpha = 0.25) +
      ggplot2::scale_fill_manual(aesthetics = "fill", values = c("red", "white"))
    
  } else {
    
    p <- ggplot2::ggplot(data = plotDat, ggplot2::aes(x = plotDat[, (2 + xPC)], y = plotDat[, (2+yPC)], 
                                                      colour = Period, group = Period)) + 
      ggplot2::stat_ellipse(type = "norm") 
    
  } 
  
  p <- p + ggplot2::facet_wrap(~identifier) +
           ggplot2::labs(x = paste0("PC", xPC, " (", round(plotDat$xVar[1], 2), "%)"),
                  y = paste0("PC", yPC, " (", round(plotDat$yVar[1], 2), "%)"),
                  fill = "") +
           ggplot2::theme_linedraw()
  
  
  return(.addTemporalDiagnostics(
    list(data = plotDat, plot = p), diagnostics
  ))
  
}

