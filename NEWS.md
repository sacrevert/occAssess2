# occAssess 1.4.0

This is the first release prepared from the maintained
[`sacrevert/occAssess2`](https://github.com/sacrevert/occAssess2) repository. It
collects the changes made since the original `robboyd/occAssess` v1.3.0 release.
The installed R package remains named `occAssess`.

## New functionality

- Added `assessSpeciesCellCounts()` to compare raw recording volume with the
  number of unique species–grid-cell combinations by identifier and period.
- Added shared, user-configurable handling of uncertain recording dates to all
  period-based assessments. Bounded and open-ended intervals can be allocated
  by documented rules, while each occurrence contributes to at most one period.
- Extended `assessRarityBias()` with optional species-level residuals,
  influence diagnostics, regression summaries and period-specific regression
  plots. Its established untransformed OLS method and default output are
  retained.
- Added `assessEnvCoverage()` for reference-environment coverage using a common
  background PCA, bivariate kernel density estimates and probability-mass
  highest-density regions.
- Exported `assessBias1D()` and added optional total variation distance and mean
  absolute error for unweighted and weighted common-bin distributions.
- Added optional raster returns to `assessSpatialCov()`.

## Corrections

- Corrected `assessSpatialCov()` so density layers contain occurrence-record
  counts, spatial cells use the requested resolution, requested empty periods
  are retained, temporal presence is calculated from raw counts, and overlap
  defaults to all requested periods.
- Replaced the fixed South American extent and approximate cell dimensions in
  `assessRepeatVisits()` with a common, exact-resolution, zero-aligned grid
  derived from retained coordinates. Records above the former plotting limit
  are no longer hidden.
- Preserved exact occurrence–environment row alignment in `assessEnvBias()`
  during spatial-uncertainty, temporal and complete-case filtering.
- Corrected the grid used for the occupied-cell proxy in
  `assessRarityBias()` so cells have the requested resolution and consistent
  alignment.
- Corrected `assessBias1D()` to retain empty bins, use valid common
  denominators, align weights with the intended rows, and reject undefined
  weighted distributions.
- Corrected `assessSpatialUncertainty()` labels and pooled annotation, retained
  observations above 15,000 in plots, and documented the actual pooled
  proportion specified.
- Added explicit validation and informative handling for missing, empty,
  degenerate and impractically large inputs across the corrected assessments.

## Important changes to numerical results

These are scientific corrections, not cosmetic changes. Analyses that must
reproduce v1.3.0 exactly should retain that version and record it in their
workflow.

- **`assessSpatialCov()`**: v1.3.0 rasterisation did not reliably count records,
  and `length(seq(..., res))` did not create cells of width `res`. Version 1.4.0
  returns actual record counts on exact-resolution, zero-aligned cells. Repeated
  coordinates, boundary observations, requested periods without records,
  `logCount`, `nPeriods`, and default overlap can therefore all produce
  different—and corrected—values.
- **`assessRarityBias()`**: occupied-cell counts now use the same exact,
  zero-aligned grid convention. Species recorded near historical grid
  boundaries can have different occupied-cell counts, changing fitted slopes,
  correlations and R-squared values while leaving the established OLS estimand
  unchanged. For the bundled regression fixture documented in the tests, one
  R-squared value changes from approximately `0.9992043` to
  `0.9985977582` because of corrected grid geometry.
- **`assessBias1D()`**: this function was not exported in the published v1.3.0
  tag, but existed in later unreleased development code. Version 1.4.0 includes
  all requested bins, excludes missing/non-finite focal values consistently,
  and correctly aligns and normalises both documented population-row weights
  and the legacy nested sample-row format. Results from that unreleased code can
  differ where bins were empty, values were missing, or weights were used.
- **`assessRepeatVisits()`**: results can change because records are no longer
  restricted to a hard-coded geographic extent and cell width now equals
  `spatRes`. The denominator remains the number of years specified in each
  period.
- **`assessEnvBias()`**: results with `maxSpatUncertainty` can change where the
  previous code paired filtered occurrences with the wrong environmental rows.
  The PCA method itself is unchanged.

## Testing, documentation and maintenance

- Replaced the obsolete CI workflow with the current five-job r-lib/actions
  matrix: R release on Linux, Windows and macOS, plus R-devel and oldrel-1 on
  Linux. Superseded runs on the same pull request are cancelled automatically.
- Added deterministic numerical regression suites for the corrected and new
  assessments, including temporal integration and row-alignment tests.
- Made the simulated-data vignette fully executable, retained historical case
  studies with explicit data requirements, and corrected GitHub installation
  instructions so built vignettes are available to `browseVignettes()`.
- Expanded help pages and the README to document estimands, missing-data rules,
  return structures, uncertainty assumptions and scientific limitations.
- Added R and GitHub citation metadata for the original package paper.
- Improved DESCRIPTION, dependency declarations, namespace generation, license
  metadata and CRAN preparation records without renaming the package.
