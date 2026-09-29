#!/usr/bin/env Rscript
#
# GGT.R -- Geo-Genetic Triangulation.
#
# Locates probable ancestor locations of a subject from (1) IBD segments that
# the subject shares with three or more genetic matches ("triangulated"
# segments) and (2) geocoded birthplaces from those matches' family trees.
# See README.md, or run `Rscript GGT.R --help`.
#
# Begg et al. (2023) Genomic analyses of hair from Ludwig van Beethoven.
# Current Biology 33, 1431-1447. https://doi.org/10.1016/j.cub.2023.02.041

# ---------------------------------------------------------------------------
# Locate this script's directory (bundled data and R/ggt.R live alongside it)
# ---------------------------------------------------------------------------
script_dir <- local({
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[1])))
  else getwd()
})
source(file.path(script_dir, "R", "ggt.R"))

# ---------------------------------------------------------------------------
# Packages
# ---------------------------------------------------------------------------
missing <- ggt_missing_packages(GGT_REQUIRED)
if (length(missing)) {
  if ("--install-missing" %in% commandArgs(trailingOnly = TRUE)) {
    lib <- Sys.getenv("R_LIBS_USER")
    dir.create(lib, showWarnings = FALSE, recursive = TRUE)
    .libPaths(c(lib, .libPaths()))
    utils::install.packages(missing, lib = lib, repos = "https://cloud.r-project.org")
    missing <- ggt_missing_packages(GGT_REQUIRED)
  }
  if (length(missing)) {
    stop("Missing R packages: ", paste(missing, collapse = ", "),
         "\nInstall them with:\n  install.packages(c(",
         paste0('"', missing, '"', collapse = ", "),
         "))\nor rerun with --install-missing.", call. = FALSE)
  }
}
suppressPackageStartupMessages(library(optparse))

# ---------------------------------------------------------------------------
# Options (names and short flags are unchanged from GGT 1.x)
# ---------------------------------------------------------------------------
opt_list <- list(
  make_option(c("-p", "--PREFIX"), type = "character",
              help = "Sample name; also the prefix of input and output files."),
  make_option(c("-w", "--WORK_DIR"), type = "character", default = ".",
              help = "Directory holding the input files [default: current directory]."),
  make_option(c("-g", "--GGT_DIR"), type = "character", default = script_dir,
              help = "GGT directory holding gis/ and map/ [default: this script's directory]."),
  make_option("--TRIANGULATED", type = "character",
              help = "Triangulated segments file [default: WORK_DIR/PREFIX_matches_triangulated.txt]."),
  make_option("--GEOCODED", type = "character",
              help = "Geocoded ancestors file [default: WORK_DIR/PREFIX_matches_geocoded.txt]."),
  make_option(c("-o", "--OUT_DIR"), type = "character",
              help = "Output directory [default: WORK_DIR/output]."),
  make_option(c("-c", "--CENTIMORGAN_MAP"), type = "character", default = "genetic_map_GRCh37.txt.gz",
              help = "Genetic map (chromosome, bp, cM; no header; may be gzipped), a path or a file in GGT_DIR/map [default: %default]."),
  make_option(c("-d", "--POP_DENSITY"), type = "character", default = "pop_density.tif",
              help = "Population density raster in long/lat, a path or a file in GGT_DIR/gis [default: %default]. Needs the terra package; used only for pop.density weights and the binned output."),
  make_option(c("-P", "--PROVINCE_LAYER"), type = "character", default = "provinces",
              help = "Province polygons with NAME_0 (country) and NAME_1 (province) columns, a path or a layer in GGT_DIR/gis [default: %default]."),
  make_option(c("-m", "--MASK_LAYER"), type = "character", default = "central_europe_mask",
              help = "Polygon(s) delimiting the analysis area, a path or a layer in GGT_DIR/gis [default: %default]."),
  make_option(c("-C", "--COUNTRIES_LIST"), type = "character", default = "country_include.txt",
              help = "Countries to map when HEX is FALSE, one per line, a path or a file in WORK_DIR [default: %default]."),
  make_option(c("-k", "--KNOWN_LOCATIONS"), type = "character", default = "known_locations.txt",
              help = "Optional known ancestor locations to plot (tab-delimited: long, lat, count), a path or a file in WORK_DIR; NULL to skip [default: %default]."),
  make_option(c("-M", "--MIN_COMMON_PLACES"), type = "integer", default = 2L,
              help = "Minimum number of segment members with an ancestor in a place [default: %default]."),
  make_option(c("-t", "--CM_THRESHOLD"), type = "double", default = 2,
              help = "Keep triangulated segments longer than this many cM [default: %default]."),
  make_option(c("-D", "--DEDUP"), type = "logical", default = TRUE,
              help = "Merge overlapping triangulated segments that share members [default: %default]."),
  make_option(c("-s", "--SEG_DUPE_THRESHOLD"), type = "double", default = 1,
              help = "Segments are merged when both their start and end lie within this many cM [default: %default]."),
  make_option(c("-u", "--USE_WEIGHTS"), type = "logical", default = FALSE,
              help = "Weight triangulated locations (see WEIGHT_TYPE) [default: %default]."),
  make_option(c("-W", "--WEIGHT_TYPE"), type = "character", default = "pop.density",
              help = "'pop.density' (inverse relative population density) or 'birth.year' (closeness to YEAR_TARGET) [default: %default]."),
  make_option(c("-f", "--WEIGHT_COEFFICIENT"), type = "double", default = 1,
              help = "Weights are raised to the power 1/WEIGHT_COEFFICIENT [default: %default]."),
  make_option(c("-y", "--YEAR_TARGET"), type = "character",
              help = "Birth year of the subject [default: the reference year]."),
  make_option(c("-Y", "--MAX_YEAR"), type = "character",
              help = "Only use ancestors born before this year [default: YEAR_TARGET - 30]."),
  make_option(c("-Z", "--YEARS_DIFF_CUTOFF"), type = "double", default = 1000,
              help = "Only use ancestors born within this many years of YEAR_TARGET [default: %default]."),
  make_option(c("-H", "--HEX"), type = "logical", default = TRUE,
              help = "Group locations in a hexagonal grid (TRUE) or by province (FALSE) [default: %default]."),
  make_option(c("-I", "--HEX_SIZE"), type = "double", default = 200000,
              help = "Hexagon width in metres [default: %default]."),
  make_option(c("-L", "--HEX_Y_RANGE"), type = "character", default = "29.5,42.5",
              help = "Standard parallels (lat1,lat2) of the Albers equal-area projection used to build the grid; NULL for the mask's latitude range [default: %default]."),
  make_option(c("-l", "--PLOT_X_RANGE"), type = "character", default = "-10,30",
              help = "Longitude range of the maps [default: %default]."),
  make_option(c("-A", "--PLOT_Y_RANGE"), type = "character", default = "37,59",
              help = "Latitude range of the maps [default: %default]."),
  make_option("--GENERAL_LOCATIONS", type = "character", default = "geocoded",
              help = "Ancestors geocoded only to a whole province: 'geocoded' (keep the geocoder's point, usually the province centroid, as in GGT 1.x), 'random' (a random point inside the province) or 'drop' [default: %default]."),
  make_option("--REFERENCE_YEAR", type = "integer",
              help = "'Present' year for birth-year filtering and imputation; set it to reproduce earlier runs [default: current year]."),
  make_option("--SEED", type = "integer", default = 1L,
              help = "Random seed [default: %default]."),
  make_option(c("-b", "--DEBUG"), type = "logical", default = TRUE,
              help = "Write diagnostic plots [default: %default]."),
  make_option(c("-q", "--QUIET"), action = "store_true", default = FALSE,
              help = "Only print warnings and errors."),
  make_option("--install-missing", action = "store_true", default = FALSE, dest = "INSTALL",
              help = "Install missing R packages into the user library."),
  make_option("--version", action = "store_true", default = FALSE, dest = "VERSION",
              help = "Print the version and exit.")
)

parser <- OptionParser(
  prog = "GGT.R",
  usage = "%prog -p PREFIX [-w WORK_DIR] [options]",
  description = "Geo-Genetic Triangulation: probable ancestor locations from triangulated IBD segments and geocoded family trees.",
  epilog = "Example:\n  Rscript GGT.R -p example1 -w example/example1 -y 1770\n\nCite: Begg et al. (2023) Current Biology 33, 1431-1447. https://doi.org/10.1016/j.cub.2023.02.041",
  option_list = opt_list)

if (!length(commandArgs(trailingOnly = TRUE))) {
  print_help(parser)
  quit(status = 0)
}
opt <- parse_args(parser)
if (opt$VERSION) {
  cat("GGT", GGT_VERSION, "\n")
  quit(status = 0)
}
options(ggt.quiet = opt$QUIET)

fail <- function(...) {
  message("ERROR: ", ...)
  quit(status = 1)
}
is_null <- function(x) is.null(x) || is.na(x) || toupper(x) %in% c("NULL", "NA", "")
num_pair <- function(x, name) {
  v <- suppressWarnings(as.numeric(strsplit(x, ",")[[1]]))
  if (length(v) != 2 || anyNA(v)) fail(name, " should be two comma-separated numbers, e.g. '29.5,42.5'")
  v
}
opt_year <- function(x, name) {
  if (is_null(x)) return(NULL)
  v <- suppressWarnings(as.numeric(x))
  if (is.na(v)) fail(name, " should be a year")
  v
}
# a path as given, or a file inside `dir` (optionally trying extra extensions)
resolve <- function(x, dir, ext = "") {
  for (cand in c(x, file.path(dir, x))) {
    for (e in ext) if (file.exists(paste0(cand, e))) return(paste0(cand, e))
  }
  NA_character_
}

if (is_null(opt$PREFIX) && (is.null(opt$TRIANGULATED) || is.null(opt$GEOCODED))) {
  fail("give -p PREFIX (or both --TRIANGULATED and --GEOCODED); see --help")
}
if (!dir.exists(opt$WORK_DIR)) fail("cannot find working directory: ", opt$WORK_DIR)
if (!dir.exists(opt$GGT_DIR)) fail("cannot find GGT directory: ", opt$GGT_DIR)
if (!opt$WEIGHT_TYPE %in% c("pop.density", "birth.year")) {
  fail("WEIGHT_TYPE should be 'pop.density' or 'birth.year'")
}
if (!opt$GENERAL_LOCATIONS %in% c("random", "geocoded", "drop")) {
  fail("GENERAL_LOCATIONS should be 'random', 'geocoded' or 'drop'")
}

prefix <- if (is_null(opt$PREFIX)) sub("_matches_triangulated.*$", "", basename(opt$TRIANGULATED)) else opt$PREFIX
gis_dir <- file.path(opt$GGT_DIR, "gis")
map_dir <- file.path(opt$GGT_DIR, "map")
tri_file <- if (is.null(opt$TRIANGULATED)) file.path(opt$WORK_DIR, paste0(prefix, "_matches_triangulated.txt")) else opt$TRIANGULATED
geo_file <- if (is.null(opt$GEOCODED)) file.path(opt$WORK_DIR, paste0(prefix, "_matches_geocoded.txt")) else opt$GEOCODED
map_file <- resolve(opt$CENTIMORGAN_MAP, map_dir, c("", ".gz"))
if (is.na(map_file)) map_file <- resolve(sub("\\.gz$", "", opt$CENTIMORGAN_MAP), map_dir, c("", ".gz"))
province_file <- resolve(opt$PROVINCE_LAYER, gis_dir, c("", ".shp", ".gpkg"))
mask_file <- resolve(opt$MASK_LAYER, gis_dir, c("", ".shp", ".gpkg"))
density_file <- if (is_null(opt$POP_DENSITY)) NA_character_ else resolve(opt$POP_DENSITY, gis_dir)
known_file <- if (is_null(opt$KNOWN_LOCATIONS)) NA_character_ else resolve(opt$KNOWN_LOCATIONS, opt$WORK_DIR)
countries_file <- if (is_null(opt$COUNTRIES_LIST)) NA_character_ else resolve(opt$COUNTRIES_LIST, opt$WORK_DIR)
out_dir <- if (is.null(opt$OUT_DIR)) file.path(opt$WORK_DIR, "output") else opt$OUT_DIR

for (f in list(c(tri_file, "triangulated segments file"), c(geo_file, "geocoded ancestors file"))) {
  if (!file.exists(f[1])) fail("cannot find ", f[2], ": ", f[1])
}
if (is.na(map_file)) fail("cannot find genetic map: ", opt$CENTIMORGAN_MAP)
if (is.na(province_file)) fail("cannot find province layer: ", opt$PROVINCE_LAYER)
if (opt$HEX && is.na(mask_file)) fail("cannot find mask layer: ", opt$MASK_LAYER)
if (!opt$HEX && is.na(countries_file)) fail("HEX = FALSE needs a countries list: ", opt$COUNTRIES_LIST)
if (!is_null(opt$KNOWN_LOCATIONS) && is.na(known_file) && opt$KNOWN_LOCATIONS != "known_locations.txt") {
  fail("cannot find known ancestor locations file: ", opt$KNOWN_LOCATIONS, " (use -k NULL to skip)")
}

use_density <- !is.na(density_file) && requireNamespace("terra", quietly = TRUE)
if (opt$USE_WEIGHTS && opt$WEIGHT_TYPE == "pop.density" && !use_density) {
  fail("pop.density weights need the terra package and a population density raster")
}

reference_year <- if (is.null(opt$REFERENCE_YEAR)) as.integer(format(Sys.Date(), "%Y")) else opt$REFERENCE_YEAR
year_target <- opt_year(opt$YEAR_TARGET, "YEAR_TARGET")
if (is.null(year_target)) year_target <- reference_year
max_year <- opt_year(opt$MAX_YEAR, "MAX_YEAR")
if (is.null(max_year)) max_year <- year_target - 30
plot_x <- num_pair(opt$PLOT_X_RANGE, "PLOT_X_RANGE")
plot_y <- num_pair(opt$PLOT_Y_RANGE, "PLOT_Y_RANGE")
hex_y <- if (is_null(opt$HEX_Y_RANGE)) NULL else num_pair(opt$HEX_Y_RANGE, "HEX_Y_RANGE")

# ---------------------------------------------------------------------------
# Run
# ---------------------------------------------------------------------------
suppressMessages(sf::sf_use_s2(FALSE))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
out <- function(name) file.path(out_dir, paste0(prefix, name))
ggt_message("GGT ", GGT_VERSION, ": processing sample '", prefix, "'")

ggt_message("Reading triangulated segments")
tri <- read_triangulated(tri_file)
n_input_segments <- nrow(tri)
ggt_message("Reading geocoded ancestors")
geo <- read_geocoded(geo_file)

ggt_message("Reading GIS layers")
provinces <- sf::st_read(province_file, quiet = TRUE)
check_columns(provinces, c("NAME_0", "NAME_1"), province_file)
countries <- if (!is.na(countries_file)) read_country_list(countries_file) else NULL

if (opt$HEX) {
  ggt_message("Making hexagonal grid")
  mask <- sf::st_read(mask_file, quiet = TRUE)
  units <- make_hex_grid(mask, opt$HEX_SIZE, hex_y)
  unit_col <- "HexName"
  units$unit <- as.character(units$HexName)
  if (opt$DEBUG) {
    dbg <- attr(units, "debug")
    grDevices::pdf(file.path(out_dir, "DEBUG_hex_extent.pdf"), width = 8.5, height = 8.5)
    plot(dbg$mask, col = "grey50", bg = "lightblue", axes = TRUE)
    plot(dbg$centers, col = "black", pch = 20, cex = 0.5, add = TRUE)
    plot(sf::st_transform(units, sf::st_crs(dbg$mask))$geometry, border = "orange", add = TRUE)
    invisible(grDevices::dev.off())
  }
} else {
  units <- province_units(provinces, countries)
  unit_col <- "unit"
}

ggt_message("Placing ancestors in ", if (opt$HEX) "hexagons" else "provinces")
geo <- place_general_locations(geo, provinces, opt$GENERAL_LOCATIONS, opt$SEED)
if (attr(geo, "n_general") > 0) {
  ggt_message("  ", attr(geo, "n_general"), " ancestors geocoded to a whole province (",
              opt$GENERAL_LOCATIONS, ")")
}
geo$unit <- assign_units(geo$Long, geo$Lat, units, "unit")
geo$NAME_1 <- if (opt$HEX) geo$unit else as.character(units$NAME_1[match(geo$unit, units$unit)])
geo$pop.density <- NA_real_
if (use_density) {
  ok <- !is.na(geo$Long) & !is.na(geo$Lat)
  r <- terra::rast(density_file)
  geo$pop.density[ok] <- terra::extract(r, cbind(geo$Long[ok], geo$Lat[ok]))[, 1]
}

ggt_message("Filtering ancestors by birth year (reference year ", reference_year, ")")
geo <- filter_birth_years(geo, reference_year)
if (attr(geo, "n_no_year") > 0) {
  ggt_message("  removed ", attr(geo, "n_no_year"), " ancestors with no birth year or relationship")
}
if (opt$DEBUG && !is.null(attr(geo, "age_model"))) {
  grDevices::pdf(file.path(out_dir, "DEBUG_years_generations.pdf"), width = 8.5, height = 8.5)
  plot(geo$Relationship, reference_year - geo$BirthYear,
       xlab = "Relationship (generations removed)", ylab = "Birth year (years before present)",
       main = "Ancestor birth year vs. generations removed")
  graphics::abline(attr(geo, "age_model"), col = "orange")
  invisible(grDevices::dev.off())
}
binned <- geo[, setdiff(names(geo), "unit")]
write_tsv(binned, out("_matches_geocoded_binned.txt"))

ggt_message("Converting segment positions to cM")
tri <- add_cm_positions(tri, read_cm_map(map_file))
tri <- tri[tri$cM > opt$CM_THRESHOLD, , drop = FALSE]
if (opt$DEDUP) {
  ggt_message("Merging overlapping segments")
  tri <- dedupe_segments(tri, opt$SEG_DUPE_THRESHOLD)
}
ggt_message("  ", nrow(tri), " of ", n_input_segments, " triangulated segments kept")

# ancestors usable for triangulation: known sample, not geocoded to a whole country
countries_only <- country_names(provinces)
ancestors <- geo[!is.na(geo$Sample) &
                   (is.na(geo$BirthPlaceGoogle) | !geo$BirthPlaceGoogle %in% countries_only), ]

ggt_message("Triangulating ancestor locations")
locs <- triangulate_locations(tri, ancestors, opt$MIN_COMMON_PLACES)
locs <- weight_and_filter(locs, year_target, max_year, opt$YEARS_DIFF_CUTOFF,
                          opt$USE_WEIGHTS, opt$WEIGHT_TYPE, opt$WEIGHT_COEFFICIENT)
unit_freq <- aggregate_units(locs)
units$Freq <- unit_freq$Freq[match(units$unit, unit_freq$unit)]
locs <- locs[locs$unit %in% units$unit, , drop = FALSE]
pts <- aggregate_points(locs)
pts <- pts[order(pts$long, pts$lat), , drop = FALSE]
ggt_message("  ", nrow(locs), " triangulated locations in ", sum(!is.na(units$Freq)),
            if (opt$HEX) " hexagons" else " provinces")

if (!nrow(locs)) {
  warning("no triangulated locations were found; check the inputs, the year filters (-y, -Y, -Z) ",
          "and --MIN_COMMON_PLACES", call. = FALSE)
}

loc_out <- locs
names(loc_out)[names(loc_out) == "unit"] <- if (opt$HEX) "HexName" else "GID_1"
write_tsv(loc_out, out("_triangulated_locations.txt"))

ggt_message("Plotting")
known <- if (!is.na(known_file)) read_known_locations(known_file) else NULL
subtitle <- sprintf("%d triangulated segments; %d triangulated locations", nrow(tri), sum(pts$Freq))
legend <- list(known = paste0(prefix, "'s\nKnown\nAncestor\nLocations"),
               fill = paste0("GGT-Inferred\nAncestor\nLocations\nPer ", if (opt$HEX) "Hexagon" else "Province"),
               colour = "GGT-Inferred\nAncestor\nLocations\nPer Point")
title <- paste0(prefix, " - Triangulated Ancestor Locations")
units_df <- polygon_df(units, units$Freq)
grDevices::pdf(out(if (opt$HEX) "_output_hex.pdf" else "_output_provinces.pdf"),
               width = 10, height = 7, useDingbats = FALSE)
print(plot_map(units_df, known = known, title = title, subtitle = subtitle,
               xlim = plot_x, ylim = plot_y, legend = legend))
invisible(grDevices::dev.off())
grDevices::pdf(out("_output_points.pdf"), width = 10, height = 7, useDingbats = FALSE)
print(plot_map(units_df, points = pts, known = known, fill = FALSE, title = title,
               subtitle = subtitle, xlim = plot_x, ylim = plot_y, legend = legend))
invisible(grDevices::dev.off())

ggt_message("Writing GeoJSON")
if (opt$HEX) {
  write_geojson(units[, c("HexName", "Freq")], out("_output_hex.geojson"))
} else {
  write_geojson(units[, intersect(c("GID_1", "NAME_0", "NAME_1", "Freq"), names(units))],
                out("_output_provinces.geojson"))
}
pts_sf <- if (nrow(pts)) {
  sf::st_as_sf(pts, coords = c("long", "lat"), crs = 4326, remove = FALSE)
} else {
  sf::st_sf(pts, geometry = sf::st_sfc(crs = 4326))
}
write_geojson(pts_sf, out("_output_points.geojson"))

params <- c(GGT_VERSION = GGT_VERSION, R = R.version.string,
            sf = as.character(utils::packageVersion("sf")),
            date = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
            triangulated = tri_file, geocoded = geo_file, genetic_map = map_file,
            provinces = province_file, mask = if (opt$HEX) mask_file else NA,
            reference_year = reference_year, year_target = year_target, max_year = max_year,
            unlist(opt[setdiff(names(opt), c("help", "INSTALL", "VERSION", "QUIET", "YEAR_TARGET", "MAX_YEAR"))]))
write_tsv(data.frame(parameter = names(params), value = unname(params)), out("_run_parameters.txt"))
ggt_message("Done. Results are in ", normalizePath(out_dir))
