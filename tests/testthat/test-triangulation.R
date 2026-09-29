test_that("birth years are cleaned and imputed", {
  geo <- data.frame(Sample = as.character(1:6),
                    Relationship = c(1, 2, 3, 4, 2, NA),
                    BirthYear = c(1970, 1940, 1910, 1300, NA, NA))
  out <- filter_birth_years(geo, reference_year = 2000, max_age = 600, recent = 10)
  expect_equal(out$Sample, c("1", "2", "3", "5"))   # 1300 too old; 6 has nothing
  expect_equal(out$BirthYear[out$Sample == "5"], 1940)  # linear in generations
  expect_equal(attr(out, "n_no_year"), 1)
})

test_that("recent birth years are treated as unknown", {
  geo <- data.frame(Sample = as.character(1:4), Relationship = c(1, 2, 3, 1),
                    BirthYear = c(1970, 1940, 1910, 1995))
  out <- filter_birth_years(geo, reference_year = 2000)
  expect_equal(out$BirthYear[4], 1970)
})

ancestors <- data.frame(
  Sample = c("a", "a", "b", "b", "c", "d"),
  unit = c("H1", "H2", "H1", "H3", "H2", NA),
  BirthYear = c(1700, 1750, 1710, 1720, 1760, 1700),
  Long = c(1, 2, 1.1, 3, 2.1, 9), Lat = c(1, 2, 1.1, 3, 2.1, 9),
  stringsAsFactors = FALSE)

test_that("locations shared by segment members are found", {
  tri <- data.frame(Chromosome = c("1", "2", "3"), Start = 1, End = 2,
                    Members = c("a, b", "a, b, c", "a, d"))
  out <- triangulate_locations(tri, ancestors, min_common = 2)
  expect_equal(out$segment, c(1, 2, 2))
  expect_equal(out$unit, c("H1", "H1", "H2"))
  expect_equal(out$Sample, c("a", "a", "a"))   # first member's ancestor is reported
  expect_equal(out$n_members, c(2L, 2L, 2L))
  expect_equal(nrow(triangulate_locations(tri, ancestors, min_common = 3)), 0)
})

test_that("an empty result has the expected columns", {
  tri <- data.frame(Chromosome = "1", Start = 1, End = 2, Members = "x, y")
  out <- triangulate_locations(tri, ancestors)
  expect_equal(nrow(out), 0)
  expect_true(all(c("unit", "BirthYear", "Long", "Lat") %in% names(out)))
})

locs <- data.frame(unit = c("H1", "H1", "H2", "H3"), BirthYear = c(1700, 1740, 1600, 1760),
                   Long = c(1, 1, 2, 3), Lat = c(1, 1, 2, 3), pop.density = c(10, 10, 100, 1))

test_that("year filters apply", {
  out <- weight_and_filter(locs, year_target = 1770, max_year = 1745, years_diff_cutoff = 150)
  expect_equal(out$unit, c("H1", "H1"))
  expect_equal(out$weight, c(1, 1))
})

test_that("birth-year weights favour ancestors born close to the target year", {
  out <- weight_and_filter(locs, 1770, 2000, 1000, use_weights = TRUE, weight_type = "birth.year")
  expect_equal(out$weight, 1 - abs(locs$BirthYear - 1770) / 170)
  expect_true(length(unique(out$weight)) > 1)
})

test_that("population-density weights favour sparsely populated places", {
  out <- weight_and_filter(locs, 1770, 2000, 1000, use_weights = TRUE, weight_type = "pop.density")
  expect_equal(out$weight, 1 / (locs$pop.density / 100))
  out2 <- weight_and_filter(locs, 1770, 2000, 1000, TRUE, "pop.density", weight_coefficient = 2)
  expect_equal(out2$weight, 1 / sqrt(locs$pop.density / 100))
})

test_that("weights are summed per unit and points are counted", {
  x <- weight_and_filter(locs, 1770, 2000, 1000)
  agg <- aggregate_units(x)
  expect_equal(agg$Freq[match(c("H1", "H2", "H3"), agg$unit)], c(2, 1, 1))
  pts <- aggregate_points(x)
  expect_equal(pts$Freq, c(2L, 1L, 1L))
})
