q <- function(...) apexog(..., data = fxdat, verbose = FALSE)

test_that("input errors are caught early", {
  expect_error(q(d ~ x1, endog = df ~ z1, gmin = 0, gmax = 0), "0/1")
  expect_error(q(y ~ x1, endog = d ~ z1, firststage = "probit", gmin = 0, gmax = 0),
               "0/1 endogenous")
  expect_error(q(y ~ x1, endog = dc ~ z1, firststage = "fprobit", gmin = 0, gmax = 0),
               "\\[0, 1\\]")
  expect_error(q(y ~ x1, endog = d ~ z1, firststage = "nonsense", gmin = 0, gmax = 0),
               "Unknown firststage")
  expect_error(q(y ~ x1, endog = d + d2 ~ z1 + z2, firststage = "biprobit",
                 gmin = 0, gmax = 0), "0/1|exactly one")
  expect_error(q(y ~ x1 + d, endog = d ~ z1, gmin = 0, gmax = 0), "both endogenous and exogenous")
  expect_error(q(y ~ x1, endog = d + d2 ~ z1, gmin = 0, gmax = 0), "at least as many")
  expect_error(q(y ~ x1, endog = d ~ z1, gmin = 0.1, gmax = -0.1), "gmin must not exceed")
  expect_error(q(y ~ x1, endog = d ~ z1, gmin = 0, gmax = 0, grid = 1), "grid")
  expect_error(q(y ~ x1, endog = d ~ z1, target = "x1", gmin = 0, gmax = 0),
               "not an endogenous")
  expect_error(q(y ~ x1, endog = d ~ z1, target = "d", aie = c("d", "x1"),
                 gmin = 0, gmax = 0), "mutually exclusive")
  expect_error(q(y ~ x1, endog = d ~ z1, aie = c("d", "x9"), gmin = 0, gmax = 0),
               "not found|interactor")
  expect_error(q(y ~ x1, endog = d ~ z1, method = "ltz", iterations = 20), "mu and omega")
  expect_error(q(y ~ x1, endog = d ~ z1, method = "ltz", mu = 0, omega = -1), "omega")
  expect_error(q(y ~ x1, endog = d ~ z1, method = "ltz", draws = matrix(0, 10, 2)),
               "one column per instrument")
  expect_error(q(y ~ x1, endog = d ~ z1, gmin = 0, gmax = 0, level = 100), "level")
  expect_error(q(y ~ x1, endog = d ~ z1, gmin = 0, gmax = 0, gunits = "foo"), "gunits")
  expect_error(q(y ~ x1, endog = d ~ z1, gmin = c(0, 0, 0), gmax = 0), "gmin")
  expect_error(q(y ~ x1, endog = d ~ z1, method = "ltz", prior = "uniform",
                 iterations = 20), "prior_params")
  expect_error(q(y ~ x1, endog = d ~ z1, gmin = 0, gmax = 0, sweep_omega = 0.01),
               "require method")
  expect_error(q(yx ~ x1 + x3 + I(x1 * x3), endog = d ~ z1, aie = c("x1", "x3"),
                 gmin = 0, gmax = 0), "already among the controls")
})

test_that("level given as a percent equals the fraction", {
  a <- q(y ~ x1, endog = d ~ z1, gmin = 0, gmax = 0, grid = 2, level = 90)
  b <- q(y ~ x1, endog = d ~ z1, gmin = 0, gmax = 0, grid = 2, level = 0.9)
  expect_equal(a$bounds, b$bounds)
})
