test_that("hexagon centers match sp::spsample(type = 'hexagonal')", {
  skip_if_not_installed("sp")
  bb <- matrix(c(0, 0, 1000, 800), 2, dimnames = list(c("x", "y"), c("min", "max")))
  poly <- sp::SpatialPolygons(list(sp::Polygons(list(sp::Polygon(
    cbind(c(0, 1000, 1000, 0, 0), c(0, 0, 800, 800, 0)))), "a")))
  ref <- sp::coordinates(sp::spsample(poly, type = "hexagonal", cellsize = 100,
                                      offset = c(0.5, 0.5)))
  ours <- hex_centers(0, 1000, 0, 800, 100)
  expect_equal(unname(ours), unname(ref))
})

test_that("hexagons are pointy-topped with the right area", {
  hex <- hex_polygons(cbind(x = 0, y = 0), dx = 100, crs = sf::NA_crs_)
  expect_equal(as.numeric(sf::st_area(hex)), sqrt(3) / 2 * 100^2)
  bb <- sf::st_bbox(hex)
  expect_equal(unname(bb[["xmax"]] - bb[["xmin"]]), 100)
})

test_that("the hexagonal grid is clipped to the mask and numbered", {
  mask <- square(5, 45, 4)
  grid <- make_hex_grid(mask, cellsize = 100000, lat_range = c(44, 50))
  expect_equal(grid$HexName, seq_len(nrow(grid)))
  expect_true(all(sf::st_is(grid, c("POLYGON", "MULTIPOLYGON"))))
  # hexagons lie within the mask; as in GGT 1.x only hexagons whose centres
  # fall inside the mask's bounding box are made, so the mask's outer edges
  # are not covered (a large share for a small square, ~0.1% for Europe)
  area <- function(x) as.numeric(sum(sf::st_area(sf::st_transform(x, 3035))))  # planar, metres
  outside <- suppressMessages(sf::st_difference(sf::st_union(grid), sf::st_geometry(mask)))
  expect_lt(area(outside), 1e-3 * area(mask))
  expect_gt(area(grid), 0.7 * area(mask))
  expect_equal(sf::st_crs(grid)$epsg, 4326L)
})

test_that("ancestors are assigned to the unit containing them", {
  units <- rbind(square(0, 0, 1, unit = "A"), square(1, 0, 1, unit = "B"))
  expect_equal(assign_units(c(0.5, 1.5, 5, NA), c(0.5, 0.5, 5, 1), units, "unit"),
               c("A", "B", NA, NA))
})

test_that("province-level locations are placed at one random point per province", {
  provinces <- rbind(square(0, 0, 1, NAME_0 = "Land", NAME_1 = "North"),
                     square(0, 5, 1, NAME_0 = "Land", NAME_1 = "South"))
  geo <- data.frame(Sample = c("1", "2", "3", "4"),
                    BirthPlaceGoogle = c("North, Land", "North, Land", "Land", NA),
                    Long = c(0.5, 0.5, 9, 3), Lat = c(0.5, 0.5, 9, 3), stringsAsFactors = FALSE)
  out <- place_general_locations(geo, provinces, "random", seed = 1)
  expect_equal(attr(out, "n_general"), 2)
  expect_equal(out$Long[1], out$Long[2])
  expect_true(out$Long[1] > 0 && out$Long[1] < 1 && out$Lat[1] > 0 && out$Lat[1] < 1)
  expect_false(isTRUE(all.equal(out$Long[1], 0.5)))
  expect_equal(plain(out)[3:4, ], geo[3:4, ])
  expect_equal(plain(place_general_locations(geo, provinces, "geocoded")), geo)
  expect_equal(nrow(place_general_locations(geo, provinces, "drop")), 2)
  expect_equal(place_general_locations(geo, provinces, "random", seed = 1)$Long, out$Long)
})

test_that("province names match regardless of the locale's encoding", {
  provinces <- square(0, 0, 1, NAME_0 = "Schweiz", NAME_1 = enc2utf8("Zürich"))
  geo <- data.frame(Sample = "1", BirthPlaceGoogle = enc2utf8("Zürich, Schweiz"),
                    Long = 0.5, Lat = 0.5, stringsAsFactors = FALSE)
  expect_equal(attr(place_general_locations(geo, provinces, "geocoded"), "n_general"), 1)
})

test_that("Google-style country names are recognised as country-only", {
  provinces <- square(0, 0, 1, NAME_0 = "Czech Republic", NAME_1 = "Praha")
  cn <- country_names(provinces)
  expect_true(all(c("Czech Republic", "Czechia", enc2utf8("C\u00f4te d'Ivoire")) %in% cn))
})
