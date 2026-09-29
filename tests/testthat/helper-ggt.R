ggt_root <- normalizePath(file.path(testthat::test_path(), "..", ".."))
source(file.path(ggt_root, "R", "ggt.R"))
suppressMessages(sf::sf_use_s2(FALSE))
options(ggt.quiet = TRUE)

# A square polygon (lon/lat) as an sf object
square <- function(xmin, ymin, size, ...) {
  ring <- cbind(c(xmin, xmin + size, xmin + size, xmin, xmin),
                c(ymin, ymin, ymin + size, ymin + size, ymin))
  geometry <- sf::st_sfc(sf::st_polygon(list(ring)), crs = 4326)
  if (...length()) sf::st_sf(data.frame(..., stringsAsFactors = FALSE), geometry = geometry)
  else sf::st_sf(geometry = geometry)
}

# drop attributes added by GGT functions before comparing data frames
plain <- function(x) {
  attr(x, "n_general") <- NULL
  x
}

# Run GGT.R in a separate R process; returns list(status, output).
# `env` is a named character vector of environment variables for the run. They
# are set in this process and inherited, because system2(env = ) prefixes
# "VAR=value" to the command, which Windows' cmd does not understand.
run_ggt <- function(args, env = character()) {
  if (length(env)) {
    old <- Sys.getenv(names(env), unset = NA, names = TRUE)
    do.call(Sys.setenv, as.list(env))
    on.exit({
      restore <- old[!is.na(old)]
      if (length(restore)) do.call(Sys.setenv, as.list(restore))
      if (any(is.na(old))) Sys.unsetenv(names(old)[is.na(old)])
    })
  }
  rscript <- file.path(R.home("bin"), "Rscript")
  out <- suppressWarnings(system2(rscript, c(shQuote(file.path(ggt_root, "GGT.R")), args),
                                  stdout = TRUE, stderr = TRUE))
  list(status = if (is.null(attr(out, "status"))) 0L else attr(out, "status"), output = out)
}
