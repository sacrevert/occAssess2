#### internal function to convert users' data to occAssess-friendly format
createData <- function(data, 
                       species,
                       x,
                       y,
                       year,
                       spatialUncertainty,
                       identifier) {

  dat <- data.frame(species = data[, species],
                    x = dplyr::pull(data, x),
                    y = dplyr::pull(data, y),
                    year = dplyr::pull(data, year),
                    spatialUncertainty = dplyr::pull(data, spatialUncertainty),
                    identifier = dplyr::pull(data, identifier))

    return(dat)
}

