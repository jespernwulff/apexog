# apexog 0.21.2: the notes of the print method and the limit of the structural
# scale.  These are display only; the estimates themselves are covered by the
# Stata-equivalence tests.

printed <- function(x) paste(utils::capture.output(print(x)), collapse = "\n")
has <- function(x, s) grepl(s, printed(x), fixed = TRUE)

test_that("the standard-error label reads two-step and vce is unchanged", {
  f <- fx_calls$A1()
  expect_true(has(f, "SE: two-step"))
  expect_false(has(f, "SE: Murphy-Topel"))
  expect_identical(f$vce, "Murphy-Topel")
  b <- apexog(yb ~ x1 + x2, endog = db ~ z1, data = fxdat, firststage = "biprobit",
              gmin = 0, gmax = 0.1, grid = 3, verbose = FALSE)
  expect_true(has(b, "SE: delta method"))
  expect_identical(b$vce, "MLE")
})

test_that("a 0/1 endogenous regressor with a linear first stage is flagged as descriptive", {
  f <- apexog(yb ~ x1 + x2, endog = db ~ z1, data = fxdat, gmin = 0, gmax = 0.1, grid = 3,
              verbose = FALSE)
  expect_identical(f$binlin, "db")
  expect_true(has(f, "0/1 endogenous regressor with a linear first stage (db):"))
  expect_true(has(f, "descriptive"))
  expect_true(has(f, "firststage = \"biprobit\" or \"probit\""))
  l <- apexog(yb ~ x1 + x2, endog = db ~ z1, data = fxdat, method = "ltz", mu = 0,
              omega = 0.01, iterations = 20, seed = 1, verbose = FALSE)
  expect_true(has(l, "0/1 endogenous regressor with a linear first stage (db):"))

  # no note with a probit first stage, the bivariate probit, or a continuous regressor
  g <- apexog(yb ~ x1 + x2, endog = db ~ z1, data = fxdat, firststage = "probit",
              gmin = 0, gmax = 0.1, grid = 3, verbose = FALSE)
  expect_length(g$binlin, 0)
  expect_false(has(g, "0/1 endogenous regressor with a linear first stage"))
  h <- apexog(yb ~ x1 + x2, endog = db ~ z1, data = fxdat, firststage = "biprobit",
              gmin = 0, gmax = 0.1, grid = 3, verbose = FALSE)
  expect_length(h$binlin, 0)
  expect_false(has(h, "0/1 endogenous regressor with a linear first stage"))
  expect_length(fx_calls$A1()$binlin, 0)

  # inside an endogenous interaction: the treatment is named, the product is
  # not, and the advice to switch first stage is dropped
  a <- apexog(ybb ~ x1 + x2, endog = db + db:x2 ~ z1 + z1:x2, data = fxdat,
              aie = c("db", "x2"), gmin = 0, gmax = 0, grid = 2, verbose = FALSE)
  expect_identical(a$binlin, "db")
  expect_true(has(a, "0/1 endogenous regressor with a linear first stage (db):"))
  expect_true(has(a, "without a causal interpretation (see ?apexog)"))
  expect_false(has(a, "firststage = \"biprobit\" or \"probit\""))
})

test_that("gs_max is |pi_z| / sd(first-stage residual) and NA outside its case", {
  f <- fx_calls$A1()                       # y ~ x1 + x2, d ~ z1
  r <- stats::lm(d ~ z1 + x1 + x2, data = fxdat)
  expect_equal(f$gs_max, abs(unname(stats::coef(r)["z1"])) / stats::sd(stats::resid(r)),
               tolerance = 1e-12)
  expect_true(is.na(fx_calls$A3()$gs_max))   # two instruments
  expect_true(is.na(fx_calls$B1()$gs_max))   # probit first stage
  expect_true(is.na(fx_calls$C7()$gs_max))   # two endogenous regressors
  expect_true(is.na(fx_calls$C5()$gs_max))   # endogenous interaction
  b <- apexog(yb ~ x1 + x2, endog = db ~ z1, data = fxdat, firststage = "biprobit",
              gmin = 0, gmax = 0.1, grid = 3, verbose = FALSE)
  expect_true(is.na(b$gs_max))
  # a 0/1 regressor with a linear first stage: not defined either, and
  # gunits = "structural" gets its own warning instead of the limit
  l <- apexog(yb ~ x1 + x2, endog = db ~ z1, data = fxdat, gmin = 0, gmax = 0.2, grid = 3,
              gunits = "structural", verbose = FALSE)
  expect_true(is.na(l$gs_max))
  expect_true(has(l, "sigma_tilde conversion not derived for a 0/1 regressor"))
  expect_false(has(l, "The conversion uses the baseline sigma_tilde"))
  expect_false(has(l, "|pi_z|/sigma_v"))
})

test_that("gunits = structural reports the baseline conversion, the limit and the warnings", {
  f <- apexog(y ~ x1 + x2, endog = d ~ z1, data = fxdat, gmin = -0.2, gmax = 0.2, grid = 3,
              gunits = "structural", verbose = FALSE)
  lim <- paste0("|pi_z|/sigma_v = ", formatC(f$gs_max, format = "f", digits = 4))
  expect_true(has(f, "The conversion uses the baseline sigma_tilde"))
  expect_true(has(f, lim))
  expect_false(has(f, "the box reaches"))

  big <- 1.5 * f$gs_max
  g <- apexog(y ~ x1 + x2, endog = d ~ z1, data = fxdat, gmin = -big, gmax = big, grid = 3,
              gunits = "structural", verbose = FALSE)
  expect_true(has(g, paste0("Warning:   the box reaches ",
                            formatC(big, format = "f", digits = 4))))
  # the warning changes nothing: same grid as the index-scale call
  gi <- apexog(y ~ x1 + x2, endog = d ~ z1, data = fxdat, gmin = -big / g$stilde,
               gmax = big / g$stilde, grid = 3, verbose = FALSE)
  expect_equal(g$grid, gi$grid, tolerance = 1e-12)
  expect_equal(g$bounds, gi$bounds, tolerance = 1e-12)

  # other units and other cases
  expect_false(has(fx_calls$A1(), "baseline sigma_tilde"))
  p <- apexog(y ~ x1 + x2, endog = d ~ z1, data = fxdat, gmin = -0.02, gmax = 0.02, grid = 3,
              gunits = "pr", verbose = FALSE)
  expect_false(has(p, "|pi_z|/sigma_v"))
  two <- apexog(y ~ x1 + x2, endog = d ~ z1 + z2, data = fxdat, gmin = -0.1, gmax = 0.1,
                grid = 2, gunits = "structural", verbose = FALSE)
  expect_true(has(two, "The conversion uses the baseline sigma_tilde"))
  expect_false(has(two, "|pi_z|/sigma_v"))
  pr <- apexog(yb ~ x1 + x2, endog = db ~ z1, data = fxdat, firststage = "probit",
               gmin = 0, gmax = 0.1, grid = 3, gunits = "structural", verbose = FALSE)
  expect_true(has(pr, "sigma_tilde conversion derived for a LINEAR first stage"))
  expect_false(has(pr, "The conversion uses the baseline sigma_tilde"))

  # LTZ: the printed share is the share of the stored draws at or beyond the limit
  narrow <- apexog(y ~ x1 + x2, endog = d ~ z1, data = fxdat, method = "ltz", mu = 0,
                   omega = (0.05 * f$gs_max)^2, gunits = "structural", iterations = 40,
                   seed = 1, verbose = FALSE)
  expect_true(has(narrow, lim))
  expect_false(has(narrow, "of the prior draws are at or beyond that limit"))
  wide <- apexog(y ~ x1 + x2, endog = d ~ z1, data = fxdat, method = "ltz", mu = 0,
                 omega = f$gs_max^2, gunits = "structural", iterations = 40, seed = 1,
                 verbose = FALSE)
  sh <- mean(abs(wide$sims$gamma[, 1]) * wide$stilde >= wide$gs_max)
  expect_gt(sh, 0)
  expect_lt(sh, 1)
  expect_true(has(wide, paste0("Warning:   ", formatC(100 * sh, format = "f", digits = 1),
                               "% of the prior draws are at or beyond that limit")))
})
