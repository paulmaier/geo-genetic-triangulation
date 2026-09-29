# Changelog

## 2.0.0 (2026)

A port to current R. With GGT 1.x's boundary files, the example dataset gives the same
hexagon grid, hexagon counts, triangulated locations and binned ancestor table as GGT 1.x.

### Data

- The boundary layers now come from geoBoundaries CGAZ (CC BY 4.0) instead of GADM 3.6:
  `gis/provinces.gpkg` and `gis/central_europe_mask.gpkg`. They are built by
  `data-raw/build_gis.R`.
- The genetic map is now the public HapMap Phase II map lifted to GRCh37, as distributed
  with Beagle (`map/genetic_map_GRCh37.txt.gz`, built by `data-raw/build_map.R`). It
  replaces an earlier GRCh37 map of unrecorded origin, which correlates with it at
  r > 0.999 on every chromosome.
- Together, these give 87 triangulated locations on the example instead of 97, in 26
  hexagons instead of 27, with the same peak hexagon; 655 rather than 665 segments pass
  the 2 cM threshold and merging. For comparison, shifting the hexagon lattice by a
  fraction of a cell gives 83–97 locations.
- 200 rather than 208 ancestors are recognised as geocoded to a whole province, because
  geoBoundaries and GADM spell some province names differently. With the default
  `--GENERAL_LOCATIONS geocoded` this affects only the count shown in the log.
- Ancestors geocoded to a whole country are recognised using the layer's country names
  plus a built-in list of common English and Google Geocoding API names. For example,
  "Jersey", "Czechia" and "Czech Republic" are all recognised.

### Compatibility

- Uses `sf` (and optionally `terra`) in place of `rgdal`, `rgeos`, `raster` and `sp`.
  `rgdal` and `rgeos` were removed from CRAN in 2023, so GGT 1.x can no longer be
  installed on current R.
- The hexagon lattice is a direct port of `sp::spsample(type = "hexagonal")`, so the grid
  is identical to 1.x. Newer GEOS versions simplify the mask outline slightly
  differently, which moves hexagon edges along coasts and borders slightly (about 0.03%
  of the area); on the example, hexagon counts are unaffected.
- `ggmap`, `geosphere` and `reshape2` (unused) and `dplyr`, `broom` and `viridis` are no
  longer needed.
- Missing packages are reported with the command to install them, rather than installed
  automatically; `--install-missing` restores the old behavior.
- The genetic map is stored gzipped.
- Place names are compared as UTF-8, so names with accents (Zürich, Baden-Württemberg)
  are matched correctly even when R runs in a non-UTF-8 locale.
- Runs in about 10 seconds instead of 60 on the example.

### Command line

- `GGT.R` has moved from `scripts/` to the repository root, and finds `gis/` and `map/`
  relative to itself, so `-g` is optional. `-w` defaults to the current directory.
- All 1.x options are unchanged. New options:
  - `--TRIANGULATED` and `--GEOCODED` give the input files directly.
  - `-o/--OUT_DIR` sets the output directory.
  - `--REFERENCE_YEAR` fixes the "present" year used in birth-year cleaning, which 1.x
    took from the system clock, for reproducible runs.
  - `--GENERAL_LOCATIONS` sets how province-level birthplaces are handled.
  - `--SEED`, `-q/--QUIET`, `--version` and `--install-missing`.
- Running without arguments prints the help. Bad arguments give a clear error and exit
  status 1.
- `-C/--COUNTRIES_LIST` is needed only for province mode (`-H FALSE`).
- Input and support files can be given as paths as well as names inside `WORK_DIR`,
  `gis/` or `map/`. Any polygon format readable by `sf` can be used for the mask and
  province layers.

### New outputs

- `PREFIX_triangulated_locations.txt`: one row per segment × triangulated location, with
  the ancestor, birth year and coordinates behind it.
- `PREFIX_run_parameters.txt`: settings, input files and software versions.

### Fixes (these can change results)

- **Weights.** Weights (`-u TRUE`) had no effect: `ifelse()` returned only the first
  location's weight, so every location got the same weight. Birth-year and
  population-density weights now vary by location as described. Locations in cells with
  zero or missing population density no longer produce infinite or missing totals.
- **Province mode (`-H FALSE`).** Provinces were matched by name alone, so provinces with
  the same name in different countries (e.g. Limburg in Belgium and the Netherlands) were
  pooled. They are now matched by their ID (`GID_1`).
- **Province-level birthplaces.**
  - With `--GENERAL_LOCATIONS random`, the random point is now used to choose the
    hexagon. 1.x chose the hexagon from the geocoder's point and used the random point
    only in the points map, population density and binned table.
  - The default, `geocoded`, keeps 1.x's hexagon assignment and uses the same point
    everywhere.
  - Random points are drawn inside the named province of the named country. 1.x could
    draw them in a same-named province elsewhere.
- **Robustness.**
  - Segments on chromosomes missing from the genetic map (e.g. X) are removed with a
    warning. Previously they became rows of missing values.
  - Positions before the first map position no longer produce missing cM values.
  - Datasets with no triangulated locations, or too few birth years to fit the birth-year
    model, no longer stop with an error.

## 1.0 (2023)

Version published with Begg et al. (2023), archived on
[Dryad](https://doi.org/10.5061/dryad.k0p2ngfc4).
