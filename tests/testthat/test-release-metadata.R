test_that("release version and installed citation are valid", {
  expect_equal(as.character(utils::packageVersion("occAssess")), "1.4.0")

  package_citation <- suppressWarnings(utils::citation("occAssess"))
  expect_s3_class(package_citation, "citation")
  formatted <- paste(format(package_citation), collapse = "\n")
  expect_match(formatted, "Boyd")
  expect_match(formatted, "10.1002/ece3.8299", fixed = TRUE)

  bibtex <- paste(utils::toBibtex(package_citation), collapse = "\n")
  expect_match(bibtex, "@Article", fixed = TRUE)
  expect_match(bibtex, "doi = {10.1002/ece3.8299}", fixed = TRUE)
})
