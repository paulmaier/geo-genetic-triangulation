# GGT: Geo-Genetic Triangulation

Geo-Genetic Triangulation (GGT) estimates where a person's autosomal ancestors lived. It
combines two sources:

1. **Genetics:** DNA segments that the person shares identical-by-descent (IBD) with
   several genetic matches.
2. **Geography:** the birthplaces recorded in those matches' family trees.

GGT was developed for the genomic analysis of Ludwig van Beethoven's hair
([Begg et al. 2023](https://doi.org/10.1016/j.cub.2023.02.041)). There it traced
Beethoven's likely ancestral origins to the Rhine valley and present-day North
Rhine-Westphalia.

From the paper:

> GGT ensures that only locations with a high likelihood of being ancestral to [a subject]
> are selected. Broadly, there are three steps: (1) segment matches between the subject
> and other genetic testers are identified; (2) genetic triangulation between three or more
> individuals identifies segments inherited identical-by-descent (IBD) from a recent common
> ancestor; (3) each triangulated segment is screened for coinciding ancestor locations,
> which increases the confidence of recent shared ancestry. When an IBD segment is
> inherited by two individuals, only one of many ancestral lines is typically shared
> between the two. The GGT method substantially reduces the noise of irrelevant ancestor
> locations.

`GGT.R` performs step 3. Its inputs are the triangulated segments from step 2 and the
geocoded ancestors of the matches.

## Contents

- [Requirements](#requirements)
- [Quick start](#quick-start)
- [Input files](#input-files)
- [How it works](#how-it-works)
- [Options](#options)
- [Output files](#output-files)
- [Bundled data](#bundled-data)
- [Changes from version 1](#changes-from-version-1)
- [Testing](#testing)
- [Citation](#citation)
- [License](#license)

## Requirements

- [R](https://www.r-project.org/) 4.1 or newer.
- R packages `optparse`, `data.table`, `sf`, `ggplot2` and `maps`.
- Optionally, `terra`, for population-density weights and the population-density column
  of the binned output.

```r
install.packages(c("optparse", "data.table", "sf", "ggplot2", "maps", "terra"))
```

On macOS and Windows, `sf` and `terra` install as ready-made binaries. On Linux they need
the GDAL, GEOS and PROJ system libraries; see the
[sf installation notes](https://r-spatial.github.io/sf/#installing). Alternatively, run
`Rscript GGT.R --install-missing ...` once to install any missing packages into your
user library.

GGT 1.x needed `rgdal`, `rgeos` and `raster`. `rgdal` and `rgeos` were removed from CRAN
in 2023, so 1.x no longer installs on current R. This version uses `sf` instead.

## Quick start

```sh
git clone https://github.com/paulmaier/geo-genetic-triangulation.git
cd geo-genetic-triangulation
Rscript GGT.R -p example1 -w example/example1 -y 1770
```

This runs GGT on the example dataset, which resembles the paper's analysis (for a subject
born in 1770). It takes about 15 seconds. The results are written to
`example/example1/output/`, and the main map is `example1_output_hex.pdf`.

Run `Rscript GGT.R --help` for all options. On macOS or Linux, you can also run
`./GGT.R ...` directly, or put the folder on your `PATH` and run `GGT.R ...` from
anywhere. The script finds its bundled data relative to its own location.

## Input files

By default GGT looks in the working directory (`-w`) for two files named after the sample
(`-p`). Paths can also be given directly with `--TRIANGULATED` and `--GEOCODED`.

### `PREFIX_matches_triangulated.txt` (required)

One row per triangulated segment, tab-delimited, with a header:

| Column | Description |
|--------|-------------|
| `Chromosome` | 1–22 (a `chr` prefix is ignored). Segments on chromosomes missing from the genetic map (e.g. X) are removed with a warning. |
| `Start`, `End` | Segment start and end, GRCh37 base-pair positions |
| `Members` | Comma-separated IDs of the matches who share the segment with the subject |
| `MemberCount`, `DistantMemberCount` | Optional. Updated when segments are merged. |
| `Mbp`, `cM` | Optional and not used; cM lengths are recalculated from the genetic map. |

### `PREFIX_matches_geocoded.txt` (required)

One row per ancestor of each match, tab-delimited, with a header:

| Column | Description |
|--------|-------------|
| `Sample` | Match ID, as used in `Members` |
| `Relationship` | Generations between the match and this ancestor |
| `BirthYear` | Ancestor's birth year; may be `NA` |
| `BirthPlaceGoogle` | Formatted birthplace returned by the [Google Geocoding API](https://developers.google.com/maps/documentation/geocoding), e.g. `Hessen, Germany`; may be `NA` |
| `Long`, `Lat` | WGS84 coordinates of the birthplace |

The example's geocoded data are not real.

### Other files (optional)

- `known_locations.txt` (`-k`): the subject's documented ancestor locations, plotted for
  comparison. Tab-delimited with columns `long`, `lat` and `count`. Use `-k NULL` to skip.
- `country_include.txt` (`-C`): countries to map when grouping by province (`-H FALSE`),
  one per line.

## How it works

1. **Segment lengths.** Each segment's start and end are converted to genetic positions
   (cM) using the bundled HapMap GRCh37 genetic map. Each position takes the value of the nearest
   map position at or below it. Segments of 2 cM or less (`-t`) are removed.
2. **Merging.** Segments on the same chromosome that share at least one member, and whose
   start and end both lie within 1 cM of each other (`-s`), are merged into one segment
   with the members pooled (`-D FALSE` turns this off).
3. **Spatial units.** The analysis area (`-m`, by default 12 Central European countries)
   is covered with a grid of hexagons 200 km wide (`-I`), in an Albers equal-area
   projection (`-L`). Alternatively, provinces can be the units (`-H FALSE`).
4. **Ancestors.**
   - Ancestors geocoded only to a country (e.g. just `Germany`) are not used.
   - Ancestors geocoded to a whole province are kept at the geocoder's point, which is
     usually the province centroid (`--GENERAL_LOCATIONS`).
   - Birth years within 10 years of the present are treated as unknown, and ancestors
     born more than 600 years ago are removed.
   - Missing birth years are estimated from a linear regression of age on generations
     removed.
5. **Triangulation.** For each segment, each member's ancestors are reduced to one per
   spatial unit. A unit is a *triangulated location* for that segment when at least 2
   members (`-M`) have an ancestor there.
6. **Filtering and weighting.** Only ancestors born before `YEAR_TARGET − 30` (`-Y`) and
   within 1000 years of `YEAR_TARGET` (`-Z`) are counted. Optionally (`-u TRUE`),
   locations are weighted by closeness in birth year to the subject, or by inverse
   population density (`-W`).
7. **Summaries.** Triangulated locations are counted per hexagon (or province) and per
   birthplace, then mapped and exported.

"Present" means the current calendar year unless you give `--REFERENCE_YEAR`. Set it to
reproduce an earlier run exactly.

A small limitation, inherited from GGT 1.x: hexagons are placed only where their centers
fall inside the analysis area's bounding box. So the outermost edges of the area can go
uncovered. For the bundled Central Europe mask this is 0.1% of the area. For a small or
rectangular custom mask it can be much more.

## Options

| Short | Long | Default | Description |
|-------|------|---------|-------------|
| `-p` | `--PREFIX` | | Sample name: the prefix of the input and output files |
| `-w` | `--WORK_DIR` | `.` | Directory with the input files |
| `-o` | `--OUT_DIR` | `WORK_DIR/output` | Output directory |
| | `--TRIANGULATED`, `--GEOCODED` | | Input files, if not named after the prefix |
| `-y` | `--YEAR_TARGET` | reference year | The subject's birth year |
| `-Y` | `--MAX_YEAR` | `YEAR_TARGET − 30` | Only use ancestors born before this year |
| `-Z` | `--YEARS_DIFF_CUTOFF` | `1000` | Only use ancestors born within this many years of `YEAR_TARGET` |
| `-M` | `--MIN_COMMON_PLACES` | `2` | Members needed with an ancestor in a unit |
| `-t` | `--CM_THRESHOLD` | `2` | Keep segments longer than this (cM) |
| `-D` | `--DEDUP` | `TRUE` | Merge overlapping segments |
| `-s` | `--SEG_DUPE_THRESHOLD` | `1` | Merge distance (cM) for segment starts and ends |
| `-H` | `--HEX` | `TRUE` | Hexagons (`TRUE`) or provinces (`FALSE`) as units |
| `-I` | `--HEX_SIZE` | `200000` | Hexagon width (m) |
| `-L` | `--HEX_Y_RANGE` | `29.5,42.5` | Standard parallels of the grid projection; `NULL` uses the mask's latitude range |
| `-u` | `--USE_WEIGHTS` | `FALSE` | Weight locations |
| `-W` | `--WEIGHT_TYPE` | `pop.density` | `pop.density` or `birth.year` |
| `-f` | `--WEIGHT_COEFFICIENT` | `1` | Weights are raised to the power 1/f |
| | `--GENERAL_LOCATIONS` | `geocoded` | Province-level birthplaces: `geocoded`, `random` (a random point in the province, seeded by `--SEED`) or `drop` |
| | `--REFERENCE_YEAR` | current year | "Present" year for birth-year cleaning |
| `-l`, `-A` | `--PLOT_X_RANGE`, `--PLOT_Y_RANGE` | `-10,30`, `37,59` | Map extent (longitude, latitude) |
| `-m` | `--MASK_LAYER` | Central Europe | Analysis area (any polygon file readable by `sf`) |
| `-P` | `--PROVINCE_LAYER` | `provinces` | Province polygons with `NAME_0` and `NAME_1` columns |
| `-c` | `--CENTIMORGAN_MAP` | `genetic_map_GRCh37.txt.gz` | Genetic map: chromosome, bp, cM; no header |
| `-d` | `--POP_DENSITY` | `pop_density.tif` | Population density raster (long/lat) |
| `-C` | `--COUNTRIES_LIST` | `country_include.txt` | Countries to map with `-H FALSE` |
| `-k` | `--KNOWN_LOCATIONS` | `known_locations.txt` | Known ancestor locations to plot; `NULL` to skip |
| `-g` | `--GGT_DIR` | script's directory | Where `gis/` and `map/` are |
| `-b` | `--DEBUG` | `TRUE` | Write diagnostic plots |
| `-q` | `--QUIET` | | Only print warnings and errors |

To analyze another region, supply a mask covering it (`-m`), and set suitable standard
parallels (`-L`) and map extents (`-l`, `-A`).

## Output files

All files are written to `OUT_DIR` and start with `PREFIX`.

| File | Contents |
|------|----------|
| `_output_hex.pdf` / `.geojson` | Triangulated locations per hexagon (`HexName`, `Freq`) |
| `_output_points.pdf` / `.geojson` | Triangulated locations per birthplace |
| `_triangulated_locations.txt` | One row per segment × triangulated location, with the ancestor used |
| `_matches_geocoded_binned.txt` | The ancestors after cleaning, with their hexagon (`NAME_1`) and population density |
| `_run_parameters.txt` | Settings, input files and software versions, for reproducibility |
| `DEBUG_hex_extent.pdf`, `DEBUG_years_generations.pdf` | The hexagonal grid, and the birth-year model |

With `-H FALSE`, the `_output_hex` files are called `_output_provinces`.

## Bundled data

| File | Contents | Source | License |
|------|----------|--------|---------|
| `gis/provinces.gpkg` | World first-level administrative divisions (`GID_1`, `NAME_0`, `NAME_1`) | [geoBoundaries](https://www.geoboundaries.org) CGAZ | CC BY 4.0 |
| `gis/central_europe_mask.gpkg` | Outline of the 12 countries analyzed in the paper | geoBoundaries CGAZ | CC BY 4.0 |
| `gis/pop_density.tif` | Population density in 2020, 2.5 arc-minute grid | [GPWv4](https://doi.org/10.7927/H49C6VHW), CIESIN, Columbia University | CC BY 4.0 |
| `map/genetic_map_GRCh37.txt.gz` | Genetic map (chromosome, GRCh37 position, cM), autosomes | HapMap Phase II, lifted to GRCh37 ([Beagle](https://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/)) | Public |

The boundary layers and the genetic map are built by
[`data-raw/build_gis.R`](data-raw/build_gis.R) and [`data-raw/build_map.R`](data-raw/build_map.R).
To use other countries, edit the country list in `build_gis.R` and rerun it. For the mask, each country is
limited to the land within 500 km of its largest part. This keeps islands such as
Corsica and Bornholm, and drops overseas territories.

Boundaries: Runfola D. et al. (2020). geoBoundaries: A global database of political
administrative boundaries. *PLoS ONE* 15(4): e0231866.
https://doi.org/10.1371/journal.pone.0231866

Genetic map: The International HapMap Consortium (2007). A second generation human
haplotype map of over 3.1 million SNPs. *Nature* 449, 851–861.
https://doi.org/10.1038/nature06258

Population density: Center for International Earth Science Information Network (CIESIN),
Columbia University (2018). Gridded Population of the World, Version 4 (GPWv4): Population
Density, Revision 11. NASA Socioeconomic Data and Applications Center (SEDAC).
https://doi.org/10.7927/H49C6VHW

## Changes from version 1

Version 2.0 runs on current R, and fixes several errors in version 1.

- **Packages.** It replaces `rgdal`/`rgeos`/`raster`/`sp` with `sf`. It no longer needs
  `dplyr`, `broom` or `viridis`, or the unused `ggmap`, `geosphere` and `reshape2`.
- **Reproducing 1.x.** With version 1's boundary files, the hexagon grid, hexagon counts,
  triangulated locations and birth-year imputation are identical to GGT 1.x.
- **Bundled data.** The bundled boundaries now come from geoBoundaries instead of GADM,
  and the genetic map is the public HapMap GRCh37 map. On the example this gives 87
  triangulated locations instead of 97, with the same peak hexagon. Shifting the hexagon
  lattice by a fraction of a cell changes the count about as much: from 83 to 97.
- **Faster.** It takes about 10 seconds instead of 60.
- **Weights.** The optional weights (`-u TRUE`) had no effect in 1.x, and now work.

See [CHANGELOG.md](CHANGELOG.md) for the details. Version 1, with its original data, is
archived on [Dryad](https://doi.org/10.5061/dryad.k0p2ngfc4).

## Testing

```sh
Rscript tests/testthat.R
```

This needs the `testthat` package. The tests include a full run of the example.

## Citation

If you use GGT, please cite:

> Begg T.J.A., Schmidt A., Kocher A., Larmuseau M.H.D., Runfeldt G., Maier P.A., Wilson
> J.D., Barquera R., Maj C., Szolek A., Sager M., Clayton S., Peltzer A., Hui R., Ronge
> J., Reiter E., Freund C., Burri M., Aron F., Tiliakou A., Osborn J., Behar D.M., Boecker
> M., Brandt G., Cleynen I., Strassburg C., Prüfer K., Kühnert D., Meredith W.R., Nöthen
> M.M., Attenborough R.D., Kivisild T., Krause J. (2023). Genomic analyses of hair from
> Ludwig van Beethoven. *Current Biology*, 33(8), 1431–1447.e22.
> https://doi.org/10.1016/j.cub.2023.02.041

Code and data: https://doi.org/10.5061/dryad.k0p2ngfc4

## License

GGT 2.0 is released under the [PolyForm Noncommercial License 1.0.0](LICENSE). It is free
for research, teaching, personal and other noncommercial use. For commercial use, please
contact the author.

Version 1, archived on [Dryad](https://doi.org/10.5061/dryad.k0p2ngfc4), remains
available under CC0 1.0. The bundled data have their own terms; see
[Bundled data](#bundled-data).
