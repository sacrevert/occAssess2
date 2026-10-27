## Current check results (2026-10-27)

Local environment:

- Debian GNU/Linux 13
- R 4.6.1 (R-devel, 2026-06-24)

Results:

- `R CMD build`: passed, including all three vignettes.
- `R CMD check --as-cran`: 0 ERRORs, 0 WARNINGs, 1 NOTE.
- Complete testthat suite: 101 expectations passed.

The remaining NOTE is:

1. CRAN incoming feasibility reports this as a new submission.

The most recent completed GitHub Actions run on `master` before this change
tested R release on Linux, Windows, and macOS plus R-devel
and oldrel-1 on Linux; all five jobs passed:
https://github.com/sacrevert/occAssess2/actions/runs/40499503103

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
