#### internal function to convert users' data to occAssess-friendly format
createData <- function(data, 
                       species,
                       x,
                       y,
                       year,
                       spatialUncertainty,
                       identifier) {

  yearValues <- if (is.null(year)) {
    rep(NA_real_, nrow(data))
  } else {
    dplyr::pull(data, year)
  }

  dat <- data.frame(species = data[, species],
                    x = dplyr::pull(data, x),
                    y = dplyr::pull(data, y),
                    year = yearValues,
                    spatialUncertainty = dplyr::pull(data, spatialUncertainty),
                    identifier = dplyr::pull(data, identifier))

    return(dat)
}

