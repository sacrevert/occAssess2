# Contributing

## Continuous integration

The full five-platform R CMD check workflow runs for pull requests, pushes to
`main` or `master`, and manual dispatches from the GitHub Actions page. When a
new commit starts the workflow for the same pull request or Git ref, GitHub
automatically cancels the superseded run. Runs for unrelated pull requests are
not affected.

For a trivial commit that does not need CI, include `[skip ci]` in the commit
message to prevent workflows triggered by `push` or `pull_request` from
starting. Use this sparingly: if branch protection requires the skipped status
check, the pull request may remain unmergeable until a later commit runs the
workflow successfully.

The workflow does not currently ignore README-only or `cran-comments.md`-only
changes. Branch-protection settings could not be verified with the repository's
GitHub Actions credentials, and a path-filtered required workflow can remain
pending and block a pull request. This is distinct from concurrency
cancellation: path or commit-message skipping prevents a run from starting,
whereas concurrency cancels an older run only after its replacement has
started.
