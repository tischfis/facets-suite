Sys.setenv("R_TESTS" = "") ## https://github.com/hadley/testthat/issues/86
suppressPackageStartupMessages(library(testthat))
suppressPackageStartupMessages(library(facets.suite))

test_check("facets.suite")
