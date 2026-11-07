# Convert the Stata example data (data-raw/example_data.do ->
# inst/extdata/*.dta) to the package data sets data/promotion.rda and
# data/certification.rda.  Run from the package root:
#   Rscript data-raw/example_data.R
stopifnot(requireNamespace("haven", quietly = TRUE))

read_clean <- function(path) {
  d <- haven::read_dta(path)
  d <- as.data.frame(haven::zap_labels(haven::zap_formats(d)))
  for (j in seq_along(d)) attr(d[[j]], "label") <- NULL
  attr(d, "label") <- NULL
  d
}

promotion <- read_clean(file.path("inst", "extdata", "promotion.dta"))
promotion$promoted <- as.integer(promotion$promoted)
promotion$manager <- as.integer(promotion$manager)
promotion$team <- as.integer(promotion$team)

certification <- read_clean(file.path("inst", "extdata", "certification.dta"))
for (v in c("promoted", "certified", "manager"))
  certification[[v]] <- as.integer(certification[[v]])

save(promotion, file = file.path("data", "promotion.rda"), version = 2,
     compress = "xz")
save(certification, file = file.path("data", "certification.rda"), version = 2,
     compress = "xz")
str(promotion)
str(certification)
