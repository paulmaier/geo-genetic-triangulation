# Geo-Genetic Triangulation (GGT) -- core functions.
#
# Sourced by GGT.R (the command-line interface) and by the tests. Each step of
# the method is a separate function so it can be tested and reused:
#
#   read_triangulated() / read_geocoded() / read_cm_map()   input
#   add_cm_positions() / dedupe_segments()                   segment cleanup
#   make_hex_grid() / province_units()                       spatial units
#   place_general_locations() / assign_units()               ancestor locations
#   filter_birth_years()                                     ancestor years
#   triangulate_locations()                                  the GGT step
#   weight_and_filter() / aggregate_units()                  summaries
#
# Spatial operations are planar (sf_use_s2(FALSE)), as in the original
# sp/rgeos implementation, so hexagon edges are straight lines in the stated
# coordinate reference system.

GGT_VERSION <- "2.0.0"

GGT_REQUIRED <- c("optparse", "data.table", "sf", "ggplot2", "maps")

# Common English and Google Geocoding API names of countries and territories.
# Together with the NAME_0 column of the province layer, an ancestor whose
# birthplace is exactly one of these names is treated as geocoded to a whole
# country and is not used.
COUNTRY_ALIASES <- c(
  "Afghanistan", "\u00c5land Islands", "American Samoa", "Antigua and Barbuda",
  "Antigua & Barbuda", "Aruba", "Bahamas", "The Bahamas", "Bermuda", "Bonaire",
  "Bosnia and Herzegovina", "Bosnia & Herzegovina", "British Virgin Islands", "Brunei",
  "Burma", "Myanmar", "Myanmar (Burma)", "Cabo Verde", "Cape Verde", "Cayman Islands",
  "Central African Republic", "Channel Islands", "Congo", "Republic of the Congo",
  "Congo - Brazzaville", "Democratic Republic of the Congo", "Congo - Kinshasa",
  "C\u00f4te d'Ivoire", "Cote d'Ivoire", "Ivory Coast", "Cura\u00e7ao", "Czechia",
  "Czech Republic", "East Timor", "Timor-Leste", "England", "Scotland", "Wales",
  "Northern Ireland", "Great Britain", "Britain", "United Kingdom", "UK", "Eswatini",
  "Swaziland", "Falkland Islands", "Faroe Islands", "French Guiana", "French Polynesia",
  "Gambia", "The Gambia", "Gibraltar", "Greenland", "Guadeloupe", "Guam", "Guernsey",
  "Holland", "Hong Kong", "Isle of Man", "Jersey", "Kosovo", "Laos", "Macao", "Macau",
  "Macedonia", "North Macedonia", "Martinique", "Mayotte", "Micronesia", "Moldova",
  "Montserrat", "New Caledonia", "North Korea", "South Korea", "Korea",
  "Northern Mariana Islands", "Palestine", "Palestinian Territories", "Puerto Rico",
  "R\u00e9union", "Reunion", "Russia", "Russian Federation", "Saint Barth\u00e9lemy",
  "Saint Kitts and Nevis", "St Kitts & Nevis", "Saint Lucia", "St Lucia",
  "Saint Martin", "Sint Maarten", "Saint Pierre and Miquelon",
  "Saint Vincent and the Grenadines", "St Vincent & Grenadines",
  "S\u00e3o Tom\u00e9 and Pr\u00edncipe", "S\u00e3o Tom\u00e9 & Pr\u00edncipe",
  "Sao Tome and Principe", "Syria", "Taiwan", "Tanzania", "Trinidad and Tobago",
  "Trinidad & Tobago", "Turkey", "T\u00fcrkiye", "Turks and Caicos Islands",
  "Turks & Caicos Islands", "United States", "United States of America", "USA",
  "US Virgin Islands", "Vatican City", "Holy See", "Vietnam", "Wallis and Futuna",
  "Western Sahara")

#' Names that mark an ancestor as geocoded only to a whole country.
country_names <- function(provinces) {
  unique(enc2utf8(c(as.character(provinces$NAME_0), COUNTRY_ALIASES)))
}
GGT_OPTIONAL <- c("terra")

ggt_message <- function(...) {
  if (!isTRUE(getOption("ggt.quiet"))) message(...)
}

ggt_missing_packages <- function(packages = GGT_REQUIRED) {
  packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
}

# Split a "Members" field ("205, 3297") into a character vector of IDs.
split_members <- function(x) {
  lapply(strsplit(as.character(x), ",", fixed = TRUE), function(m) {
    m <- trimws(m)
    m[nzchar(m)]
  })
}

normalize_chrom <- function(x) {
  sub("^chr", "", as.character(x), ignore.case = TRUE)
}


# --------------------------------------------------------------------------- #
# Input
# --------------------------------------------------------------------------- #

check_columns <- function(df, required, path) {
  missing <- setdiff(required, names(df))
  if (length(missing)) {
    stop(sprintf("'%s' is missing required column(s): %s", path,
                 paste(missing, collapse = ", ")), call. = FALSE)
  }
}

#' Triangulated segments: one row per segment, tab-delimited, with at least
#' Chromosome, Start, End and Members (comma-separated sample IDs).
read_triangulated <- function(path) {
  tri <- data.table::fread(path, sep = "\t", colClasses = list(character = "Members"),
                           encoding = "UTF-8",
                           data.table = FALSE, showProgress = FALSE)
  check_columns(tri, c("Chromosome", "Start", "End", "Members"), path)
  tri$Chromosome <- normalize_chrom(tri$Chromosome)
  tri$Start <- as.numeric(tri$Start)
  tri$End <- as.numeric(tri$End)
  tri
}

#' Geocoded ancestors: one row per ancestor of each matching sample, with
#' Sample, Relationship, BirthYear, BirthPlaceGoogle, Long, Lat.
read_geocoded <- function(path) {
  geo <- data.table::fread(path, sep = "\t", quote = "\"", na.strings = c("NA", ""),
                           encoding = "UTF-8",
                           colClasses = list(character = c("Sample", "BirthPlaceGoogle")),
                           data.table = FALSE, showProgress = FALSE)
  check_columns(geo, c("Sample", "Relationship", "BirthYear", "BirthPlaceGoogle",
                       "Long", "Lat"), path)
  geo$BirthPlaceGoogle <- enc2utf8(geo$BirthPlaceGoogle)
  for (col in c("Relationship", "BirthYear", "Long", "Lat")) {
    geo[[col]] <- suppressWarnings(as.numeric(geo[[col]]))
  }
  geo
}

#' Genetic map: three whitespace-delimited columns without a header
#' (chromosome, bp position, cM). May be gzipped.
read_cm_map <- function(path) {
  if (grepl("\\.gz$", path)) {
    tmp <- gunzip_to_temp(path)
    on.exit(unlink(tmp))
    map <- data.table::fread(tmp, header = FALSE, showProgress = FALSE)
  } else {
    map <- data.table::fread(path, header = FALSE, showProgress = FALSE)
  }
  if (ncol(map) < 3) stop(sprintf("'%s' should have 3 columns: chromosome, bp, cM", path),
                          call. = FALSE)
  map <- data.frame(chrom = normalize_chrom(map[[1]]), pos = as.numeric(map[[2]]),
                    cm = as.numeric(map[[3]]), stringsAsFactors = FALSE)
  map[order(map$chrom, map$pos), ]
}

# Decompress a .gz file to a temporary file (fread needs R.utils to do this itself).
gunzip_to_temp <- function(path) {
  tmp <- tempfile(fileext = ".txt")
  src <- gzfile(path, "rb")
  dst <- file(tmp, "wb")
  on.exit({ close(src); close(dst) })
  repeat {
    chunk <- readBin(src, "raw", 2^24)
    if (!length(chunk)) break
    writeBin(chunk, dst)
  }
  tmp
}

read_known_locations <- function(path) {
  known <- utils::read.delim(path, stringsAsFactors = FALSE, encoding = "UTF-8")
  names(known) <- tolower(names(known))
  check_columns(known, c("long", "lat"), path)
  if (!"count" %in% names(known)) known$count <- 1
  known
}

read_country_list <- function(path) {
  x <- enc2utf8(trimws(readLines(path, warn = FALSE, encoding = "UTF-8")))
  x[nzchar(x) & !startsWith(x, "#")]
}


# --------------------------------------------------------------------------- #
# Segments
# --------------------------------------------------------------------------- #

#' Look up the genetic position (cM) of each segment's start and end.
#'
#' Each position takes the cM value of the nearest map position at or below it
#' (positions before the first map entry take the first entry's value).
#' Segment length in cM is recomputed as end - start.
add_cm_positions <- function(tri, map) {
  tri$cM.start <- NA_real_
  tri$cM.end <- NA_real_
  for (chrom in unique(tri$Chromosome)) {
    rows <- which(tri$Chromosome == chrom)
    m <- map[map$chrom == chrom, ]
    if (!nrow(m)) next
    lookup <- function(pos) m$cm[pmax(findInterval(pos, m$pos), 1L)]
    tri$cM.start[rows] <- lookup(tri$Start[rows])
    tri$cM.end[rows] <- lookup(tri$End[rows])
  }
  unmapped <- is.na(tri$cM.start) | is.na(tri$cM.end)
  if (any(unmapped)) {
    warning(sprintf("%d segment(s) on chromosomes absent from the genetic map (%s) were removed",
                    sum(unmapped), paste(unique(tri$Chromosome[unmapped]), collapse = ", ")),
            call. = FALSE)
    tri <- tri[!unmapped, , drop = FALSE]
  }
  tri$cM <- tri$cM.end - tri$cM.start
  tri
}

#' Merge overlapping triangulated segments.
#'
#' Each segment in turn gathers the not-yet-merged segments on its chromosome
#' whose cM start and end both lie within `threshold` cM of its own and that
#' share at least one member with it (itself included). Each group becomes one
#' segment: the group's first row, with the members of the whole group pooled.
#' A segment that is already in a group can still gather further segments, as
#' in GGT 1.x.
dedupe_segments <- function(tri, threshold) {
  n <- nrow(tri)
  if (n == 0) return(tri)
  members <- split_members(tri$Members)
  grouped <- logical(n)
  groups <- list()
  for (i in seq_len(n)) {
    cand <- which(!grouped &
                    tri$Chromosome == tri$Chromosome[i] &
                    abs(tri$cM.start - tri$cM.start[i]) < threshold &
                    abs(tri$cM.end - tri$cM.end[i]) < threshold)
    cand <- cand[vapply(members[cand], function(m) any(m %in% members[[i]]), logical(1))]
    if (length(cand)) {
      groups[[length(groups) + 1]] <- cand
      grouped[cand] <- TRUE
    }
  }
  out <- tri[vapply(groups, `[`, integer(1), 1), , drop = FALSE]
  pooled <- lapply(groups, function(g) unique(unlist(members[g])))
  out$Members <- vapply(pooled, paste, character(1), collapse = ", ")
  if ("DistantMemberCount" %in% names(out)) out$DistantMemberCount <- lengths(pooled)
  if ("MemberCount" %in% names(out)) out$MemberCount <- lengths(pooled) + 1L
  rownames(out) <- NULL
  out
}


# --------------------------------------------------------------------------- #
# Spatial units
# --------------------------------------------------------------------------- #

#' Hexagon centers on a triangular lattice covering a bounding box.
#'
#' Port of sp::spsample(type = "hexagonal", offset = c(0.5, 0.5)) so the grid
#' is identical to the original implementation.
hex_centers <- function(xmin, xmax, ymin, ymax, dx, offset = c(0.5, 0.5)) {
  dy <- sqrt(3) * dx / 2
  x <- seq(xmin, xmax - dx / 2, dx)
  y <- seq(ymin, ymax, dy)
  y <- rep(y, each = length(x))
  x <- rep(c(x, x + dx / 2), length.out = length(y))
  x <- x + (xmax - max(x)) / 2 + offset[1] * dx
  y <- y + (ymax - max(y)) / 2 + offset[2] * dy
  keep <- x >= xmin & x <= xmax & y >= ymin & y <= ymax
  cbind(x = x[keep], y = y[keep])
}

#' Pointy-topped hexagons of width `dx` around each center (sp's
#' HexPoints2SpatialPolygons).
hex_polygons <- function(centers, dx, crs) {
  dy <- dx / sqrt(3)
  xo <- c(-dx / 2, 0, dx / 2, dx / 2, 0, -dx / 2, -dx / 2)
  yo <- c(dy / 2, dy, dy / 2, -dy / 2, -dy, -dy / 2, dy / 2)
  polys <- lapply(seq_len(nrow(centers)), function(i) {
    sf::st_polygon(list(cbind(centers[i, 1] + xo, centers[i, 2] + yo)))
  })
  sf::st_sfc(polys, crs = crs)
}

#' Hexagonal grid clipped to the analysis mask.
#'
#' The mask is simplified (tolerance in degrees), projected to an Albers
#' equal-area projection with standard parallels `lat_range`, covered with
#' hexagons of width `cellsize` metres, clipped, and returned in WGS84 with a
#' HexName column numbering the hexagons 1..n.
make_hex_grid <- function(mask, cellsize, lat_range = NULL, simplify_tol = 0.05) {
  mask <- suppressMessages(sf::st_union(sf::st_geometry(mask)))
  if (is.null(lat_range)) {
    bb <- sf::st_bbox(sf::st_transform(mask, 4326))
    lat_range <- c(bb[["ymin"]], bb[["ymax"]])
  }
  aea <- sprintf("+proj=aea +lat_1=%s +lat_2=%s", lat_range[1], lat_range[2])
  mask_proj <- suppressWarnings(suppressMessages(
    sf::st_simplify(mask, preserveTopology = FALSE, dTolerance = simplify_tol)))
  mask_proj <- sf::st_transform(mask_proj, aea)
  bb <- sf::st_bbox(mask_proj)
  centers <- hex_centers(bb[["xmin"]] - 500, bb[["xmax"]] + 500,
                         bb[["ymin"]] - 500, bb[["ymax"]] + 500, cellsize)
  hexes <- hex_polygons(centers, cellsize, sf::st_crs(mask_proj))
  clipped <- suppressWarnings(suppressMessages(sf::st_intersection(hexes, mask_proj)))
  hit <- attr(clipped, "idx")[, 1]
  if (any(sf::st_is(clipped, "GEOMETRYCOLLECTION"))) {
    clipped <- sf::st_collection_extract(clipped, "POLYGON", warn = FALSE)
  }
  keep <- as.numeric(sf::st_area(clipped)) > 0
  grid <- sf::st_sf(HexName = seq_len(sum(keep)),
                    geometry = sf::st_transform(clipped[keep], 4326))
  pts <- sf::st_sfc(lapply(seq_len(nrow(centers)), function(i) sf::st_point(centers[i, ])),
                    crs = sf::st_crs(mask_proj))
  inside <- lengths(sf::st_intersects(pts, mask_proj)) > 0
  attr(grid, "debug") <- list(mask = mask_proj, hexes = hexes, centers = pts[inside])
  grid
}

#' Province polygons as spatial units (used with --HEX FALSE).
province_units <- function(provinces, countries = NULL) {
  units <- provinces
  if (!is.null(countries)) units <- units[enc2utf8(as.character(units$NAME_0)) %in% countries, ]
  id <- if ("GID_1" %in% names(units)) units$GID_1 else paste(units$NAME_0, units$NAME_1)
  units$unit <- as.character(id)
  units
}

#' Assign each ancestor to the spatial unit containing its coordinates.
assign_units <- function(long, lat, units, id_col) {
  out <- rep(NA_character_, length(long))
  ok <- !is.na(long) & !is.na(lat)
  if (!any(ok)) return(out)
  pts <- sf::st_as_sf(data.frame(long = long[ok], lat = lat[ok]), coords = c("long", "lat"),
                      crs = 4326)
  hits <- suppressMessages(sf::st_intersects(pts, sf::st_transform(units, 4326)))
  first <- vapply(hits, function(h) if (length(h)) h[1] else NA_integer_, integer(1))
  out[ok] <- as.character(units[[id_col]][first])
  out
}

#' Handle ancestors geocoded only to a whole province ("Bavaria, Germany").
#'
#' mode = "geocoded": keep the geocoder's coordinates (usually the province
#'   centroid); what GGT 1.x did when assigning hexagons.
#' mode = "random": move them to a random point inside the province (one point
#'   per province, reproducible via `seed`).
#' mode = "drop": remove them.
place_general_locations <- function(geo, provinces, mode = "geocoded", seed = 1) {
  full_names <- enc2utf8(paste0(provinces$NAME_1, ", ", provinces$NAME_0))
  general <- !is.na(geo$BirthPlaceGoogle) & geo$BirthPlaceGoogle %in% full_names
  attr(geo, "n_general") <- sum(general)
  if (!any(general) || mode == "geocoded") return(geo)
  if (mode == "drop") return(geo[!general, , drop = FALSE])
  for (place in unique(geo$BirthPlaceGoogle[general])) {
    shape <- sf::st_geometry(provinces[full_names == place, ])[1]
    pt <- random_point(shape, seed)
    rows <- which(general & geo$BirthPlaceGoogle == place)
    geo$Long[rows] <- pt[1]
    geo$Lat[rows] <- pt[2]
  }
  geo
}

random_point <- function(shape, seed, tries = 1000) {
  set.seed(seed)
  bb <- sf::st_bbox(shape)
  for (i in seq_len(tries)) {
    p <- c(stats::runif(1, bb[["xmin"]], bb[["xmax"]]), stats::runif(1, bb[["ymin"]], bb[["ymax"]]))
    pt <- sf::st_sfc(sf::st_point(p), crs = sf::st_crs(shape))
    if (lengths(suppressMessages(sf::st_intersects(pt, shape))) > 0) return(p)
  }
  as.numeric(sf::st_coordinates(suppressWarnings(sf::st_point_on_surface(shape))))[1:2]
}


# --------------------------------------------------------------------------- #
# Ancestor birth years
# --------------------------------------------------------------------------- #

#' Clean ancestor birth years.
#'
#' Birth years within 10 years of `reference_year` are treated as unknown,
#' ancestors born more than `max_age` years before it are removed, and missing
#' years are imputed from a linear model of age on generations removed
#' (Relationship). Ancestors with neither a year nor a relationship are removed.
filter_birth_years <- function(geo, reference_year, max_age = 600, recent = 10) {
  y <- geo$BirthYear
  y[!is.na(y) & y > reference_year - recent] <- NA
  geo$BirthYear <- y
  geo <- geo[!(!is.na(y) & y < reference_year - max_age), , drop = FALSE]
  age <- reference_year - geo$BirthYear
  known <- !is.na(age) & !is.na(geo$Relationship)
  model <- NULL
  if (sum(known) >= 2 && length(unique(geo$Relationship[known])) >= 2) {
    model <- stats::lm(age ~ Relationship, data = data.frame(age = age, Relationship = geo$Relationship)[known, ])
    impute <- is.na(geo$BirthYear) & !is.na(geo$Relationship)
    geo$BirthYear[impute] <- reference_year -
      (stats::coef(model)[1] + stats::coef(model)[2] * geo$Relationship[impute])
  } else {
    warning("too few ancestors with both a birth year and a relationship to impute birth years",
            call. = FALSE)
  }
  n_removed <- sum(is.na(geo$BirthYear))
  geo <- geo[!is.na(geo$BirthYear), , drop = FALSE]
  attr(geo, "age_model") <- model
  attr(geo, "n_no_year") <- n_removed
  geo
}


# --------------------------------------------------------------------------- #
# Triangulation of locations
# --------------------------------------------------------------------------- #

#' Find ancestor locations shared by the members of each triangulated segment.
#'
#' For each segment, every member's ancestors are reduced to one entry per
#' spatial unit (the first listed ancestor in that unit). A unit counts as a
#' triangulated location for the segment when at least `min_common` members
#' have an ancestor there. Returns one row per segment x location.
triangulate_locations <- function(tri, ancestors, min_common = 2) {
  anc <- ancestors[!is.na(ancestors$unit), , drop = FALSE]
  # first ancestor per (sample, unit), in input order
  anc <- anc[!duplicated(paste(anc$Sample, anc$unit, sep = "\r")), , drop = FALSE]
  by_sample <- split(seq_len(nrow(anc)), anc$Sample)
  has_ancestors <- unique(ancestors$Sample)
  members <- split_members(tri$Members)
  rows <- vector("list", nrow(tri))
  for (i in seq_len(nrow(tri))) {
    m <- members[[i]]
    m <- m[m %in% has_ancestors]
    if (length(m) < 2) next
    idx <- unlist(by_sample[m], use.names = FALSE)
    if (!length(idx)) next
    hits <- anc[idx, , drop = FALSE]
    n_members <- table(hits$unit)
    first <- hits[!duplicated(hits$unit), , drop = FALSE]
    first <- first[as.integer(n_members[first$unit]) >= min_common, , drop = FALSE]
    if (!nrow(first)) next
    rows[[i]] <- data.frame(
      segment = i, Chromosome = tri$Chromosome[i], Start = tri$Start[i], End = tri$End[i],
      unit = first$unit, n_members = as.integer(n_members[first$unit]),
      Sample = first$Sample, BirthYear = first$BirthYear,
      Long = first$Long, Lat = first$Lat,
      pop.density = if ("pop.density" %in% names(first)) first$pop.density else NA_real_,
      stringsAsFactors = FALSE)
  }
  out <- do.call(rbind, rows)
  if (is.null(out)) {
    out <- data.frame(segment = integer(), Chromosome = character(), Start = numeric(),
                      End = numeric(), unit = character(), n_members = integer(),
                      Sample = character(), BirthYear = numeric(), Long = numeric(),
                      Lat = numeric(), pop.density = numeric(), stringsAsFactors = FALSE)
  }
  out
}

#' Apply year filters and weights to triangulated locations.
#'
#' Keeps locations whose ancestor was born less than `years_diff_cutoff` years
#' from `year_target` and before `max_year`. With `use_weights`, each location
#' is weighted by closeness in birth year to the target ("birth.year") or by
#' inverse relative population density ("pop.density").
weight_and_filter <- function(locs, year_target, max_year, years_diff_cutoff,
                              use_weights = FALSE, weight_type = "pop.density",
                              weight_coefficient = 1) {
  locs$years.diff <- abs(locs$BirthYear - year_target)
  locs <- locs[!is.na(locs$years.diff), , drop = FALSE]
  locs$weight <- rep(1, nrow(locs))
  if (use_weights && nrow(locs)) {
    if (weight_type == "birth.year") {
      w <- 1 - locs$years.diff / max(locs$years.diff)
      locs$weight <- w^(1 / weight_coefficient)
    } else {
      d <- locs$pop.density
      if (all(is.na(d))) stop("population density weights requested but no densities available",
                              call. = FALSE)
      d[!is.na(d) & d <= 0] <- min(d[!is.na(d) & d > 0])
      rel <- d / max(d, na.rm = TRUE)
      locs$weight <- 1 / rel^(1 / weight_coefficient)
      if (any(is.na(rel))) {
        warning(sprintf("%d location(s) without population density were given weight 0",
                        sum(is.na(rel))), call. = FALSE)
        locs$weight[is.na(rel)] <- 0
      }
    }
  }
  locs[locs$years.diff < years_diff_cutoff & locs$BirthYear < max_year, , drop = FALSE]
}

#' Sum of location weights per spatial unit.
aggregate_units <- function(locs) {
  if (!nrow(locs)) return(data.frame(unit = character(), Freq = numeric()))
  agg <- stats::aggregate(weight ~ unit, data = locs, FUN = sum)
  names(agg) <- c("unit", "Freq")
  agg
}

#' Number of triangulated locations at each distinct coordinate.
aggregate_points <- function(locs) {
  if (!nrow(locs)) return(data.frame(long = numeric(), lat = numeric(), Freq = integer()))
  key <- paste(locs$Long, locs$Lat)
  first <- !duplicated(key)
  data.frame(long = locs$Long[first], lat = locs$Lat[first],
             Freq = as.integer(table(key)[key[first]]))
}


# --------------------------------------------------------------------------- #
# Output
# --------------------------------------------------------------------------- #

# sf polygons -> data frame for ggplot2::geom_polygon (holes via subgroup)
polygon_df <- function(x, value = NULL) {
  geom <- sf::st_geometry(sf::st_cast(sf::st_transform(x, 4326), "MULTIPOLYGON", warn = FALSE))
  xy <- as.data.frame(sf::st_coordinates(geom))
  df <- data.frame(long = xy$X, lat = xy$Y,
                   group = paste(xy$L3, xy$L2, sep = "."),
                   subgroup = paste(xy$L3, xy$L2, xy$L1, sep = "."))
  if (!is.null(value)) df$value <- value[xy$L3]
  df
}

plot_map <- function(units_df, points = NULL, known = NULL, fill = TRUE, title, subtitle,
                     xlim, ylim, legend) {
  world <- ggplot2::map_data("world")
  g <- ggplot2::ggplot(world, ggplot2::aes(x = long, y = lat, group = group)) +
    ggplot2::geom_polygon(fill = "white", colour = "black")
  if (fill) {
    g <- g + ggplot2::geom_polygon(data = units_df, colour = "black", linewidth = 0.25,
                                   ggplot2::aes(fill = value, subgroup = subgroup)) +
      ggplot2::scale_fill_viridis_c(option = "plasma", na.value = "grey50")
  } else {
    g <- g + ggplot2::geom_polygon(data = units_df, colour = "black", linewidth = 0.25,
                                   fill = "gray60", ggplot2::aes(subgroup = subgroup))
  }
  g <- g + ggplot2::geom_polygon(fill = NA, colour = "black")
  if (!is.null(known) && nrow(known)) {
    g <- g + ggplot2::geom_point(data = known, inherit.aes = FALSE, shape = 21, fill = "white",
                                 ggplot2::aes(x = long, y = lat, size = count)) +
      ggplot2::scale_size_continuous(name = legend$known)
  }
  if (!is.null(points) && nrow(points)) {
    g <- g + ggplot2::geom_point(data = points, inherit.aes = FALSE,
                                 ggplot2::aes(x = long, y = lat, colour = Freq)) +
      ggplot2::scale_colour_viridis_c(option = "plasma")
  }
  labels <- ggplot2::labs(x = "Longitude", y = "Latitude", title = title, subtitle = subtitle)
  if (fill) labels$fill <- legend$fill
  if (!is.null(points) && nrow(points)) labels$colour <- legend$colour
  g + ggplot2::coord_fixed(ratio = 1.5, xlim = xlim, ylim = ylim) + labels +
    ggplot2::theme(panel.background = ggplot2::element_rect(fill = "#d4f5ff"),
                   plot.title = ggplot2::element_text(size = 24, hjust = 0.5, vjust = 1),
                   plot.subtitle = ggplot2::element_text(size = 18, hjust = 0.5))
}

write_geojson <- function(x, path) {
  if (file.exists(path)) file.remove(path)
  sf::st_write(x, path, driver = "GeoJSON", quiet = TRUE)
}

write_tsv <- function(x, path) {
  utils::write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
}
