# Run the GGT tests from the repository root with:
#   Rscript tests/testthat.R
library(testthat)
test_dir(file.path("tests", "testthat"), reporter = "summary", stop_on_failure = TRUE)
