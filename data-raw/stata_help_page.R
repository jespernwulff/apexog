# Build vignettes/articles/stata-help.Rmd from the Stata help file rendered to
# text with Stata's smcl2txt translator (data-raw/apexog_help.txt):
#   translate apexog.sthlp apexog_help.txt, translator(smcl2txt)
# Run from the package root: Rscript data-raw/stata_help_page.R
txt <- readLines(file.path("data-raw", "apexog_help.txt"), warn = FALSE,
                 encoding = "UTF-8")
# drop the Stata logo block at the top
start <- grep("^\\s*Title\\s*$", txt)[1]
txt <- txt[start:length(txt)]
# remove the 6-space indentation of the translator
txt <- sub("^      ", "", txt)
# four {opt name(} directives that smcl2txt leaves unrendered
txt <- gsub("{opt fir:ststage(}type)", "firststage(type)   ", txt, fixed = TRUE)
txt <- gsub("{opt aie(}X Z)", "aie(X Z)      ", txt, fixed = TRUE)
txt <- gsub("{opt vce(}vcetype)", "vce(vcetype)      ", txt, fixed = TRUE)
txt <- gsub("{opt graphopt:ions(}string)", "graphoptions(string)      ", txt,
            fixed = TRUE)
stopifnot(!any(grepl("{opt", txt, fixed = TRUE)))
stopifnot(!any(grepl("```", txt, fixed = TRUE)))
while (length(txt) && !nzchar(trimws(txt[length(txt)]))) txt <- txt[-length(txt)]

out <- c(
  "---",
  'title: "Stata help file"',
  "---",
  "",
  "The help file of the Stata command, version 0.21.0, as Stata displays it",
  "with `help apexog` (rendered to text with Stata's `smcl2txt` translator).",
  "",
  "```{.stata-output}",
  txt,
  "```")
writeLines(out, file.path("vignettes", "articles", "stata-help.Rmd"),
           useBytes = TRUE)
cat("wrote vignettes/articles/stata-help.Rmd:", length(txt), "lines of help\n")
