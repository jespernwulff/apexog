# Convert the Stata reference results written by
#   data-raw/stata_fixtures.do (run in Stata from the package root)
# into tests/testthat/fixtures/stata_fixtures.rds, so the R test suite needs
# neither Stata nor haven.  Run from the package root:
#   Rscript data-raw/make_fixtures.R
# Structure of the .rds:
#   list(data = <data.frame>,
#        specs = list(A1 = list(m = list(<matrix>...), s = c(<scalars>),
#                               mac = c(<macros>)), ...))

stopifnot(requireNamespace("haven", quietly = TRUE))
fxdir <- file.path("data-raw", "stata_fixtures")
files <- list.files(fxdir, pattern = "\\.dta$")

strip <- function(df) {
  df <- as.data.frame(haven::zap_labels(haven::zap_formats(df)))
  for (j in seq_along(df)) attr(df[[j]], "label") <- NULL
  df
}

dat <- strip(haven::read_dta(file.path(fxdir, "fx_data.dta")))

tags <- unique(sub("^fx_([^_]+)_.*$", "\\1", setdiff(files, "fx_data.dta")))
specs <- list()
for (tg in tags) {
  sp <- list(m = list(), s = numeric(0), mac = character(0))
  mf <- files[grepl(paste0("^fx_", tg, "_m_"), files)]
  for (f in mf) {
    nm <- sub(paste0("^fx_", tg, "_m_(.*)\\.dta$"), "\\1", f)
    df <- strip(haven::read_dta(file.path(fxdir, f)))
    sp$m[[nm]] <- unname(as.matrix(df))
  }
  sc <- strip(haven::read_dta(file.path(fxdir, paste0("fx_", tg, "_scalars.dta"))))
  sp$s <- stats::setNames(sc$value, sc$name)
  mc <- strip(haven::read_dta(file.path(fxdir, paste0("fx_", tg, "_macros.dta"))))
  sp$mac <- stats::setNames(as.character(mc$value), mc$name)
  specs[[tg]] <- sp
}

out <- list(data = dat, specs = specs)
saveRDS(out, file.path("tests", "testthat", "fixtures", "stata_fixtures.rds"),
        version = 2)
cat("wrote", length(specs), "specs:", paste(names(specs), collapse = " "), "\n")
