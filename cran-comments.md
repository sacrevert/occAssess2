## Test environments

- Local: Debian GNU/Linux 13, R 4.6.1 (R-devel at 2026-10-26)
- GitHub Actions: R release on Linux, Windows, and macOS; R-devel and
  oldrel-1 on Linux

## Initial check (unmodified master, 2026-10-26)

`R CMD check --as-cran` completed with 0 ERRORs, 0 WARNINGs, and 6 NOTEs:

1. CRAN incoming feasibility identified a new submission, a non-template MIT
   license file, an invalid ORCID placeholder, and a non-title-case Title.
2. `.github` was included in the source package.
3. The MIT license stub was not valid DCF, and the ORCID placeholder was
   invalid.
4. `devtools` and `maps` were declared but not imported or otherwise detected
   as used.
5. Static code analysis found unqualified functions and tidy-evaluation plot
   variables.
6. `occAssessReproducible.Rmd` had no vignette engine, and `vignettes/figure`
   looked like a leftover knitr directory.

The original testthat suite passed before changes (17 expectations).

## Final local check (2026-10-26)

Both direct `R CMD check --as-cran` and `rcmdcheck` completed on Debian 13
with R 4.6.1: 0 ERRORs, 0 WARNINGs, and 2 NOTEs. The only package issue in
either NOTE is the invalid ORCID placeholder described below; the incoming
check also reports the expected "New submission" status. The complete
testthat suite still passes (17 expectations).
