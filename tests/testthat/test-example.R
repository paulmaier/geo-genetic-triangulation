# End-to-end run of GGT.R on the bundled example (about 10-20 seconds).

example_dir <- file.path(ggt_root, "example", "example1")

test_that("the example gives the expected results", {
  out <- tempfile("ggt-")
  # fixed reference year so the birth-year cleaning does not depend on the clock;
  # LC_ALL=C checks that accented place names still match in a non-UTF-8 locale
  res <- run_ggt(c("-p", "example1", "-w", shQuote(example_dir), "-o", shQuote(out), "-y", "1770",
                   "--REFERENCE_YEAR", "2026"),
                 env = c(LC_ALL = "C"))
  expect_equal(res$status, 0L, info = paste(res$output, collapse = "\n"))
  log <- paste(res$output, collapse = "\n")
  expect_match(log, "200 ancestors geocoded to a whole province")
  expect_match(log, "655 of 1116 triangulated segments kept")
  expect_match(log, "87 triangulated locations in 26 hexagons")

  files <- c("_output_hex.pdf", "_output_points.pdf", "_output_hex.geojson",
             "_output_points.geojson", "_matches_geocoded_binned.txt",
             "_triangulated_locations.txt", "_run_parameters.txt")
  expect_true(all(file.exists(file.path(out, paste0("example1", files)))))

  hex <- sf::st_read(file.path(out, "example1_output_hex.geojson"), quiet = TRUE)
  expect_equal(nrow(hex), 72)
  expect_equal(sort(hex$Freq, decreasing = TRUE)[1:6], c(19, 13, 8, 8, 7, 5))
  expect_equal(sum(hex$Freq, na.rm = TRUE), 87)
  locs <- read.delim(file.path(out, "example1_triangulated_locations.txt"))
  expect_equal(nrow(locs), 87)
  expect_equal(length(unique(locs$segment)), 43)
})

test_that("province mode and weights run", {
  out <- tempfile("ggt-")
  res <- run_ggt(c("-p", "example1", "-w", shQuote(example_dir), "-o", shQuote(out), "-y", "1770",
                   "-H", "FALSE", "-u", "TRUE", "-W", "birth.year", "-b", "FALSE"))
  expect_equal(res$status, 0L, info = paste(res$output, collapse = "\n"))
  prov <- sf::st_read(file.path(out, "example1_output_provinces.geojson"), quiet = TRUE)
  expect_true(all(c("GID_1", "NAME_1", "Freq") %in% names(prov)))
  expect_gt(sum(!is.na(prov$Freq)), 10)
})

test_that("a run with no triangulated locations completes with a warning", {
  out <- tempfile("ggt-")
  res <- run_ggt(c("-p", "example1", "-w", shQuote(example_dir), "-o", shQuote(out), "-y", "1770",
                   "-M", "50", "-b", "FALSE"))
  log <- paste(res$output, collapse = "\n")
  expect_equal(res$status, 0L, info = log)
  expect_match(log, "no triangulated locations were found")
  expect_false(grepl("Inf|-Inf", log))
  expect_true(file.exists(file.path(out, "example1_output_points.geojson")))
})

test_that("bad arguments give a clear error", {
  res <- run_ggt(c("-p", "nothing", "-w", shQuote(tempdir())))
  expect_equal(res$status, 1L)
  expect_match(paste(res$output, collapse = "\n"), "cannot find triangulated segments file")
  res <- run_ggt(c("-p", "example1", "-w", shQuote(example_dir), "-W", "bogus"))
  expect_equal(res$status, 1L)
  expect_match(paste(res$output, collapse = "\n"), "WEIGHT_TYPE")
})

test_that("--help and --version work", {
  expect_match(paste(run_ggt("--help")$output, collapse = "\n"), "PREFIX")
  expect_match(paste(run_ggt("--version")$output, collapse = "\n"), "GGT 2")
})
