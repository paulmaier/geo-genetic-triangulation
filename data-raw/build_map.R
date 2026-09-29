#!/usr/bin/env Rscript
#
# Build GGT's bundled genetic map, map/genetic_map_GRCh37.txt.gz, from the
# HapMap Phase II genetic map lifted to GRCh37 (Adam Auton), in the PLINK
# format distributed with Beagle (Brian Browning).
#
# Usage (from the repository root; downloads about 50 MB):
#   Rscript data-raw/build_map.R

url <- "https://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/plink.GRCh37.map.zip"
dl <- file.path(tempdir(), "plink.GRCh37.map.zip")
options(timeout = max(3600, getOption("timeout")))
if (!file.exists(dl)) utils::download.file(url, dl, mode = "wb", quiet = TRUE)
dir <- file.path(tempdir(), "plink_maps")
utils::unzip(dl, exdir = dir)

out <- gzfile("map/genetic_map_GRCh37.txt.gz", "w", compression = 9)
for (chrom in 1:22) {
  # PLINK map columns: chromosome, marker, cM, bp -> GGT: chromosome, bp, cM
  m <- utils::read.table(file.path(dir, sprintf("plink.chr%d.GRCh37.map", chrom)),
                         colClasses = c("character", "NULL", "character", "character"))
  writeLines(paste(m[[1]], m[[3]], m[[2]]), out)
}
close(out)
message("Wrote map/genetic_map_GRCh37.txt.gz")
