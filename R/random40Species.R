#' Simulated occurrence data.
#'
#' Simulated occurrence data for 40 species in the UK in the format required by occAssess. The data are simulated to have random distributions in space and time. 
#' Record prevalence varies among species by a factor of five with a mean of
#' approximately 1,000 records. All records are resolved to species level.
#'
#' @docType data
#'
#' @usage data(random40Species)
#'
#' @format A data frame with 40,000 rows and six columns: \code{species},
#'   \code{x}, \code{y}, \code{year}, \code{spatialUncertainty}, and
#'   \code{identifier}. Coordinates use a projected UK grid in metres.
#'
#' @keywords datasets
#'
#' @examples
#' data(random40Species)

"random40Species"
