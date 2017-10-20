test_that("run facets.suite", {
  samplefile <- system.file(package = "facets.suite", "example_data/samples.txt")
  maffile <- system.file(package = "facets.suite", "example_data/tcga.maf")
  outdir <- system.file(package = "facets.suite", "example_data/example_output")
  arg_line <- paste("-s", samplefile, "-m", maffile, "-o", outdir)
  cat("Rscript -e 'facets.suite::facets.suite()'", arg_line, "\n")
  facets.suite(arg_line)
})
