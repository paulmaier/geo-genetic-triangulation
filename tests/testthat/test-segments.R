test_that("members are split and trimmed", {
  expect_equal(split_members(c("205, 3297", "1,2 ,3", "")),
               list(c("205", "3297"), c("1", "2", "3"), character()))
})

test_that("genetic positions use the nearest map position at or below", {
  map <- data.frame(chrom = c("1", "1", "1", "2"), pos = c(100, 200, 300, 100),
                    cm = c(1, 2, 3, 10))
  tri <- data.frame(Chromosome = c("1", "1", "2"), Start = c(150, 50, 100), End = c(350, 250, 500),
                    Members = "a, b")
  out <- add_cm_positions(tri, map)
  expect_equal(out$cM.start, c(1, 1, 10))   # 50 is before the first position
  expect_equal(out$cM.end, c(3, 2, 10))
  expect_equal(out$cM, c(2, 1, 0))
})

test_that("segments on chromosomes missing from the map are dropped with a warning", {
  map <- data.frame(chrom = "1", pos = c(1, 1000), cm = c(0, 5))
  tri <- data.frame(Chromosome = c("1", "X"), Start = 1, End = 500, Members = "a, b")
  expect_warning(out <- add_cm_positions(tri, map), "absent from the genetic map")
  expect_equal(nrow(out), 1)
})

test_that("chr prefixes are ignored", {
  expect_equal(normalize_chrom(c("chr1", "CHRX", "7")), c("1", "X", "7"))
})

dedupe_input <- function() {
  data.frame(Chromosome = c("1", "1", "1", "1", "2"),
             cM.start = c(10, 10.5, 10.9, 30, 10), cM.end = c(20, 20.4, 20.8, 40, 20),
             Members = c("a, b", "b, c", "c, d", "a, b", "a, b"),
             MemberCount = 3L, DistantMemberCount = 2L, stringsAsFactors = FALSE)
}

test_that("overlapping segments sharing members are merged", {
  out <- dedupe_segments(dedupe_input(), threshold = 1)
  # rows 1-2 merge (same region, share b); row 3 is within 1 cM of row 2 and
  # shares c, and is gathered by row 2 as in GGT 1.x
  expect_equal(nrow(out), 4)
  expect_equal(out$Members, c("a, b, c", "c, d", "a, b", "a, b"))
  expect_equal(out$DistantMemberCount, c(3L, 2L, 2L, 2L))
  expect_equal(out$MemberCount, c(4L, 3L, 3L, 3L))
})

test_that("segments without shared members are not merged", {
  x <- dedupe_input()[1:2, ]
  x$Members <- c("a, b", "c, d")
  expect_equal(nrow(dedupe_segments(x, 1)), 2)
})
