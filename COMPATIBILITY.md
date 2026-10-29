# occAssess 1.4.0 compatibility report

This audit compares the release candidate with the published
`robboyd/occAssess` tag `v1.3.0` (`aacc515046a5e889693c20bd1ef38bc46f08e672`).
It covers exported names, formal arguments and defaults, default return
contracts, and known numerical changes.

## Public API

No v1.3.0 export has been removed or renamed. The following exports are new:

- `assessBias1D()`;
- `assessEnvCoverage()`;
- `assessSpeciesCellCounts()`.

All arguments added to an existing exported function are optional and follow
the v1.3.0 arguments:

| Function(s) | Optional arguments added | Default effect |
| --- | --- | --- |
| All ten period-based occurrence assessments | `dateHandling = NULL` | Retains exact-year assignment |
| `assessSpatialCov()` | `returnRaster = FALSE` | Retains plotted return mode |
| `assessRarityBias()` | `returnDiagnostics = FALSE`, `plotRegressions = FALSE`, `labelTopResiduals = 0` | Retains the period-index plot and data |

The new `assessBias1D()` has `returnMetrics = FALSE`, preserving the two-component
result used by its post-v1.3.0 development implementation. No established
v1.3.0 argument default has changed.

## Return compatibility

- Existing default calls continue to return their established named list
  components. Optional rarity diagnostics, regression plots, rasters and
  one-dimensional discrepancy metrics appear only when requested.
- Period-based results gain a `temporalDiagnostics` attribute. This does not add
  or reorder list components.
- `assessSpatialUncertainty()` retains its one-column numeric matrix and basic
  two-component list; descriptive dimnames and a periods attribute were added.
- The three newly exported functions establish new return contracts and do not
  replace an existing assessment.

Code that requires objects to have no extra attributes, rather than using their
documented components, should be reviewed. No other material default-structure
break was identified.

## Numerical compatibility

The following corrections can change results for identical inputs and should be
treated as corrected numerical behavior rather than interface breakage:

- `assessSpatialCov()` now counts occurrence rows, uses exact zero-aligned
  `res` cells, includes requested empty periods, separates raw presence from
  log-transformed display counts, and bases default overlap on all requested
  periods.
- `assessRarityBias()` uses exact zero-aligned `res` cells for its recorded-range
  proxy. The OLS response, predictor, transformation and five-species minimum
  are otherwise unchanged.
- `assessRepeatVisits()` uses an input-derived common grid rather than a fixed
  South American extent and approximate raster dimensions.
- `assessEnvBias()` applies occurrence filters to environmental rows in exact
  correspondence and removes incomplete environmental rows across all
  variables.
- `assessBias1D()` uses a complete common set of bins, valid denominators and
  row-aligned weights. This function was not exported by the published v1.3.0
  tag, but results from later unreleased code may differ.
- `assessSpatialUncertainty()` retains the historical pooled proportion
  specified; its corrected plot label, annotation placement and plot range do
  not redefine that statistic.

Input validation is intentionally stricter for invalid periods, resolutions,
coordinates, indicators, weights and empty selections. Inputs that previously
failed obscurely, recycled values, or produced undefined outputs may now stop
with an explicit message or return a documented unavailable status.

## Version assessment

Version 1.4.0 is appropriate because the release predominantly adds optional
functionality while retaining exported names, argument order, defaults and
default list components. The spatial and weighting fixes materially affect
some numerical results, but correct implementations that did not match their
documented scientific definitions. No unintended public-interface break
requiring a major version was found.
