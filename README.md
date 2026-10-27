![R-CMD-check](https://github.com/sacrevert/occAssess2/actions/workflows/check-standard.yml/badge.svg)

# occAssess

<code>occAssess</code> screens species-occurrence data for common forms of
temporal, taxonomic, spatial, and environmental bias. The GitHub repository is
named **occAssess2**, but the installed R package and the name used with
<code>library()</code> are both **occAssess**.

## Installation

The package is not yet on CRAN. Install the current GitHub version with:

~~~r
install.packages("remotes") # needed only if remotes is not already installed
remotes::install_github(
  "sacrevert/occAssess2",
  dependencies = TRUE,
  build_vignettes = TRUE
)
library(occAssess)
~~~

The `build_vignettes = TRUE` option is required for
`browseVignettes("occAssess")` after installation from GitHub because
`remotes` does not build package vignettes by default. Building the vignettes
also requires Pandoc; it is included with RStudio and can otherwise be
[installed separately](https://pandoc.org/installing.html).

## Quick, reproducible example

This example uses <code>simDat</code>, a simulated occurrence dataset bundled
with the package. It requires no local files, downloads, or network access after
installation.

~~~r
library(occAssess)
data("simDat", package = "occAssess")

periods <- list(
  2001:2002,
  2003:2004,
  2005:2006,
  2007:2008,
  2009:2010
)

record_count <- assessRecordNumber(
  dat = simDat,
  species = "species",
  x = "x",
  y = "y",
  year = "year",
  spatialUncertainty = "spatialUncertainty",
  identifier = "identifier",
  periods = periods
)

record_count$data
#>     val   group Period
#> 1  1054 survey1      1
#> 2  1062 survey1      2
#> 3  1129 survey1      3
#> 4  1088 survey1      4
#> 5  1019 survey1      5
#> 6  1064 survey2      1
#> 7  1072 survey2      2
#> 8  1126 survey2      3
#> 9  1050 survey2      4
#> 10 1105 survey2      5

record_count$plot
~~~

<code>assessSpeciesNumber()</code> uses the same period and identifier
structure:

~~~r
species_count <- assessSpeciesNumber(
  dat = simDat,
  species = "species",
  x = "x",
  y = "y",
  year = "year",
  spatialUncertainty = "spatialUncertainty",
  identifier = "identifier",
  periods = periods
)

species_count$data
#>    val   group Period
#> 1   11 survey1      1
#> 2   11 survey1      2
#> 3   11 survey1      3
#> 4   11 survey1      4
#> 5   11 survey1      5
#> 6   11 survey2      1
#> 7   11 survey2      2
#> 8   11 survey2      3
#> 9   11 survey2      4
#> 10  11 survey2      5

species_count$plot
~~~

Both functions return a list containing the numerical result in
<code>$data</code> and a ggplot object in <code>$plot</code>. Here,
<code>identifier</code> separates two simulated surveys and
<code>Period</code> refers to the position of the corresponding vector in
<code>periods</code>.

## What can the package assess?

- <code>assessRecordNumber()</code>: number of occurrence records per period.
- <code>assessSpeciesNumber()</code>: number of recorded species per period.
- <code>assessSpeciesID()</code>: count or proportion of records identified to species.
- <code>assessRarityBias()</code>: association between recorded range size and
  record count.
- <code>assessSpatialCov()</code>: gridded spatial coverage by period.
- <code>assessSpatialBias()</code>: nearest-neighbour index relative to random
  samples drawn within a supplied raster mask.
- <code>assessEnvBias()</code>: occurrence locations in principal-component
  environmental space, optionally relative to background environmental data.
- <code>assessSpatialUncertainty()</code>: distribution of reported spatial
  uncertainty and the pooled proportion of records for which it is specified.
- <code>assessBias1D()</code>: comparison of sample and population distributions
  for one variable, optionally after weighting.
- <code>assessRepeatVisits()</code>: occurrence-record counts per grid cell per
  specified period-year, as a proxy for repeated spatial sampling.

Most occurrence-data functions require columns representing species,
coordinates, year, spatial uncertainty, and an identifier. Coordinate and
uncertainty units must be internally consistent. Spatial resolutions and masks
must use the same coordinate reference system as the occurrence coordinates.
See each function's help page for its exact missing-data and return-value
behavior.

## Vignettes and case studies

- [Fully reproducible simulated-data vignette](vignettes/occAssessReproducible.Rmd):
  an executable walkthrough using only bundled data.
- [Historical South American case study](vignettes/occAssess.Rmd): a static
  example based on GBIF occurrence data available from
  [GBIF download 10.15468/dl.wqz6z3](https://doi.org/10.15468/dl.wqz6z3).
  Its spatial and environmental sections additionally require user-supplied
  boundary, mask, and [WorldClim](https://www.worldclim.org/) data.
- [Historical random-UK case study](vignettes/occAssessUKRandom40SpeciesVaryingRecNumbers.Rmd):
  a static example whose occurrence and environmental tables are bundled, but
  whose original mapping sections used an unavailable private boundary package
  and an external raster mask.

After installation, list the built vignettes with:

~~~r
browseVignettes("occAssess")
~~~

Current CRAN-preparation results and known submission notes are recorded in
[cran-comments.md](cran-comments.md).

## License

MIT © UK Centre for Ecology and Hydrology.
