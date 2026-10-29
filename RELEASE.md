# occAssess release checklist

The release-preparation pull request must be reviewed and merged before any of
these publication actions. Do not tag a branch commit that has not passed the
required checks.

- [ ] Review and merge the `Prepare occAssess 1.4.0 release` pull request into
      `master`.
- [ ] On the resulting `master` commit, confirm that `DESCRIPTION` contains
      `Package: occAssess` and `Version: 1.4.0`.
- [ ] Confirm that all required GitHub Actions jobs pass on that exact final
      commit: macOS release, Windows release, Ubuntu release, Ubuntu R-devel and
      Ubuntu oldrel-1.
- [ ] Review `NEWS.md`, `inst/CITATION`, `CITATION.cff`, and the paper DOI
      <https://doi.org/10.1002/ece3.8299>.
- [ ] Create an annotated tag named `v1.4.0` pointing to the validated release
      commit, then push that tag. Do not move the tag after publication.
- [ ] Publish a GitHub release titled `occAssess 1.4.0`, using the `NEWS.md`
      1.4.0 section as the basis of the release description.
- [ ] Verify the published tag and a clean tagged installation:
      `remotes::install_github("sacrevert/occAssess2", ref = "v1.4.0",
      dependencies = TRUE, build_vignettes = TRUE)`; then run
      `packageVersion("occAssess")`, `citation("occAssess")`, the README example,
      and `browseVignettes("occAssess")`.
- [ ] Prepare the initial CRAN submission as a separate, explicitly authorised
      task. Review the final CRAN policies, `cran-comments.md`, and checks again;
      do not assume that publishing the GitHub release submits to CRAN.
- [ ] After the stable GitHub release is published, consider advancing the
      development branch to `1.4.0.9000` in a separate commit, following R
      package development conventions. The release commit itself must remain
      `1.4.0`.
