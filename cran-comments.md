## Release-candidate check results (2026-10-29)

Local environment:

- Ubuntu 22.04.5 LTS
- R 4.6.1 (R-devel, 2026-06-24)

Results:

- `R CMD build`: passed, including all three vignettes.
- `R CMD check --as-cran`: 0 ERRORs, 0 WARNINGs, 1 NOTE with the CRAN incoming
  service enabled. All package, example, test, vignette, PDF-manual and
  HTML-manual checks passed.
- Complete testthat suite: 461 expectations passed.
- Installation of the built source package into a clean library: passed.
- Package loading, the README example, all three installed vignettes,
  `citation("occAssess")`, and `toBibtex(citation("occAssess"))`: passed.

The single NOTE is expected for an initial CRAN submission:

1. CRAN incoming feasibility identifies `occAssess` as a new submission.

GitHub Actions results for the release-preparation commit will be recorded in
the pull request after the workflow completes. The workflow continues to cover
R release on Linux, Windows, and macOS, plus R-devel and oldrel-1 on Linux.

## Initial check record (unmodified master, 2026-10-26)

The first local `R CMD check --as-cran` completed with 0 ERRORs, 0 WARNINGs,
and 6 NOTEs:

1. CRAN incoming feasibility identified a new submission, a non-template MIT
   license file, an invalid placeholder ORCID, and a non-title-case Title.
2. `.github` was included in the source package.
3. The MIT license stub was not valid DCF, and the placeholder ORCID was
   invalid.
4. `devtools` and `maps` were declared but not imported or otherwise
   detected as used.
5. Static code analysis found unqualified functions and tidy-evaluation plot
   variables.
6. `occAssessReproducible.Rmd` had no vignette engine, and
   `vignettes/figure` looked like a leftover knitr directory.

The original testthat suite passed 17 expectations. These observations are
retained as historical baseline information; they do not describe the current
package metadata or check state.
