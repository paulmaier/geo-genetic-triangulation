#!/usr/bin/env Rscript
#
# Build GGT's bundled boundary layers from geoBoundaries CGAZ (CC BY 4.0).
#
#   gis/provinces.gpkg            world first-level divisions: GID_1, NAME_0, NAME_1
#   gis/central_europe_mask.gpkg  outline of the 12 countries analysed in Begg et al. (2023)
#
# Usage (from the repository root; downloads about 300 MB):
#   Rscript data-raw/build_gis.R [download_dir]
#
# Source: Runfola D. et al. (2020) geoBoundaries: A global database of political
# administrative boundaries. PLoS ONE 15(4): e0231866.
# https://www.geoboundaries.org

suppressPackageStartupMessages(library(sf))
suppressMessages(sf_use_s2(FALSE))

args <- commandArgs(trailingOnly = TRUE)
dl <- if (length(args)) args[1] else file.path(tempdir(), "geoboundaries")
dir.create(dl, showWarnings = FALSE, recursive = TRUE)
base <- "https://github.com/wmgeolab/geoBoundaries/raw/main/releaseData/CGAZ"

fetch <- function(name) {
  path <- file.path(dl, name)
  if (!file.exists(path)) {
    message("Downloading ", name)
    options(timeout = max(3600, getOption("timeout")))
    utils::download.file(file.path(base, name), path, mode = "wb", quiet = TRUE)
  }
  path
}

# Countries of the Central Europe analysis area (ISO 3166-1 alpha-3)
MASK_COUNTRIES <- c("AUT", "BEL", "CZE", "DNK", "FRA", "DEU", "HUN", "LUX", "NLD", "POL",
                    "CHE", "SVK")
# Simplification tolerance (degrees) for the bundled layers. GGT simplifies the
# mask by 0.05 degrees before building the hexagonal grid, so this is far finer
# than anything the analysis uses.
TOLERANCE <- 0.005

# Keep the parts of a country within MAX_DISTANCE_KM of its largest part. This
# drops distant territories (e.g. French overseas departments, the Faroe Islands
# in Denmark) and keeps near-shore islands (e.g. Corsica, Bornholm).
MAX_DISTANCE_KM <- 500
home_territory <- function(geom) {
  parts <- st_cast(st_make_valid(geom), "POLYGON")
  parts_m <- st_transform(parts, 3035)
  main <- which.max(as.numeric(st_area(parts_m)))
  near <- as.numeric(st_distance(parts_m, parts_m[main])) <= MAX_DISTANCE_KM * 1000
  st_union(parts[near])
}

adm0 <- st_read(fetch("geoBoundariesCGAZ_ADM0.gpkg"), quiet = TRUE)
adm1 <- st_read(fetch("geoBoundariesCGAZ_ADM1.gpkg"), quiet = TRUE)

country <- setNames(as.character(adm0$shapeName), adm0$shapeGroup)
provinces <- st_sf(
  GID_1 = as.character(adm1$shapeID),
  GID_0 = as.character(adm1$shapeGroup),
  NAME_0 = unname(country[as.character(adm1$shapeGroup)]),
  NAME_1 = as.character(adm1$shapeName),
  geometry = st_geometry(adm1))
provinces <- provinces[!is.na(provinces$NAME_0), ]
provinces <- suppressWarnings(st_simplify(st_make_valid(provinces), dTolerance = TOLERANCE,
                                          preserveTopology = TRUE))
provinces <- provinces[!st_is_empty(provinces), ]

mask <- adm0[adm0$shapeGroup %in% MASK_COUNTRIES, ]
stopifnot(nrow(mask) == length(MASK_COUNTRIES))
mask <- do.call(c, lapply(seq_len(nrow(mask)), function(i) home_territory(st_geometry(mask)[i])))
mask <- st_make_valid(st_union(st_make_valid(mask)))
mask <- st_sf(name = "Central Europe",
              geometry = suppressWarnings(st_simplify(mask, dTolerance = TOLERANCE,
                                                      preserveTopology = TRUE)))

provinces <- st_transform(provinces, 4326)
mask <- st_transform(mask, 4326)
dir.create("gis", showWarnings = FALSE)
st_write(provinces, "gis/provinces.gpkg", quiet = TRUE, delete_dsn = TRUE)
st_write(mask, "gis/central_europe_mask.gpkg", quiet = TRUE, delete_dsn = TRUE)
message(sprintf("Wrote gis/provinces.gpkg (%d divisions, %d countries) and gis/central_europe_mask.gpkg",
                nrow(provinces), length(unique(provinces$NAME_0))))
