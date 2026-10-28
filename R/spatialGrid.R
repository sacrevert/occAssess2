.occAssessCellBins <- function(x, y, res) {
  cellBin <- function(coordinate) {
    scaled <- coordinate / res
    if (any(!is.finite(scaled))) {
      stop("Coordinates are too large relative to res to construct a grid.")
    }
    nearestInteger <- round(scaled)
    tolerance <- pmin(
      sqrt(.Machine$double.eps),
      8 * .Machine$double.eps * pmax(1, abs(scaled))
    )
    onBoundary <- abs(scaled - nearestInteger) <= tolerance
    scaled[onBoundary] <- nearestInteger[onBoundary]
    floor(scaled)
  }

  list(x = cellBin(x), y = cellBin(y))
}
