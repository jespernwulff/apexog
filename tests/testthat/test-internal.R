test_that("Stata rounding is half away from zero", {
  expect_identical(apexog:::.stata_round(c(2.5, 3.5, -2.5, 2.4999, 0.5)),
                   c(3, 4, -3, 2, 1))
})

test_that("percentile rule follows _apexogsim_pct", {
  x <- c(5, 1, 4, 2, 3, NA, 10, 7, 9, 6, 8)
  # n = 10 valid: lo_pos = max(1, round(0.25)) = 1, hi_pos = min(10, round(9.75)) = 10
  expect_equal(apexog:::.pct_interval(x, 0.95), c(1, 10))
  y <- seq_len(100)
  # n = 100: lo = round(2.5) = 3 (half away from zero), hi = round(97.5) = 98
  expect_equal(apexog:::.pct_interval(y, 0.95), c(3, 98))
  expect_identical(apexog:::.pct_interval(c(NA, NA), 0.95), c(NA_real_, NA_real_))
})

test_that("analytic biprobit Hessian equals the numerical one", {
  S <- list(y = fxdat$yb, D = cbind(db = fxdat$db),
            bp_x1 = cbind(fxdat$x1, fxdat$x2, fxdat$db, 1),
            bp_x2 = cbind(fxdat$z1, fxdat$zb, fxdat$x1, fxdat$x2, 1))
  th <- c(0.3, 0.2, 0.5, -0.1, 0.8, 0.2, 0.3, 0.05, -0.2, 0.4)
  off <- 0.05 * fxdat$z1
  Ha <- apexog:::.bp_derivs(th, S, off, hess = TRUE)$negH
  Hn <- apexog:::.bp_neghess_num(th, S, off)
  expect_lt(max(abs(Ha - Hn)) / max(abs(Hn)), 1e-8)
})

test_that("gunits pr and structural rescale the gamma range exactly", {
  base <- apexog(y ~ x1 + x2, endog = d ~ z1, data = fxdat, gmin = -0.01, gmax = 0.01,
                 grid = 3, verbose = FALSE)
  pr <- apexog(y ~ x1 + x2, endog = d ~ z1, data = fxdat, gunits = "pr",
               gmin = -0.01 * base$mphi, gmax = 0.01 * base$mphi, grid = 3, verbose = FALSE)
  st <- apexog(y ~ x1 + x2, endog = d ~ z1, data = fxdat, gunits = "structural",
               gmin = -0.01 * base$stilde, gmax = 0.01 * base$stilde, grid = 3,
               verbose = FALSE)
  expect_equal(pr$grid[[1]], base$grid[[1]], tolerance = 1e-10)
  expect_equal(st$grid[[1]], base$grid[[1]], tolerance = 1e-10)
  expect_equal(pr$bounds, base$bounds, tolerance = 1e-8)
  idx <- apexog(y ~ x1 + x2, endog = d ~ z1, data = fxdat, gunits = "index",
                gmin = -0.01, gmax = 0.01, grid = 3, verbose = FALSE)
  expect_equal(idx$bounds, base$bounds)
})

test_that("a degenerate LTZ prior collapses on the UCI fit at the same gamma", {
  u <- apexog(y ~ x1 + x2, endog = d ~ z1, data = fxdat, gmin = 0.05, gmax = 0.05,
              grid = 2, verbose = FALSE)
  l <- apexog(y ~ x1 + x2, endog = d ~ z1, data = fxdat, method = "ltz", mu = 0.05,
              omega = 1e-16, iterations = 20, seed = 3, verbose = FALSE)
  expect_equal(l$prior_band$lower, u$grid$est_d[1], tolerance = 1e-6)
  expect_equal(l$prior_band$upper, u$grid$est_d[1], tolerance = 1e-6)
})

test_that("LTZ is reproducible under a seed and every prior type runs", {
  a <- apexog(y ~ x1 + x2, endog = d ~ z1, data = fxdat, method = "ltz", mu = 0,
              omega = 0.0025, iterations = 20, seed = 11, verbose = FALSE)
  b <- apexog(y ~ x1 + x2, endog = d ~ z1, data = fxdat, method = "ltz", mu = 0,
              omega = 0.0025, iterations = 20, seed = 11, verbose = FALSE)
  expect_identical(a$bounds, b$bounds)
  pri <- list(normal = c(0, 0.05), uniform = c(-0.05, 0.05), chi2 = 1, t = 3,
              gamma = c(1, 0.02), poisson = 1, special = c(-0.02, 0, 0.03))
  for (p in names(pri)) {
    f <- apexog(y ~ x1 + x2, endog = d ~ z1, data = fxdat, method = "ltz", prior = p,
                prior_params = list(pri[[p]]), iterations = 10, seed = 1, verbose = FALSE)
    expect_s3_class(f, "apexog")
    # chi2 / t / poisson draws are several latent units: some refits hit
    # fitted probabilities of exactly 0/1 and are excluded (as in Stata)
    expect_lt(f$n_failed, 10)
  }
  f <- apexog(y ~ x1 + x2, endog = d ~ z1 + z2, data = fxdat, method = "ltz",
              prior = function(n) cbind(stats::runif(n, -0.02, 0.02), 0),
              iterations = 10, seed = 1, verbose = FALSE)
  expect_true(all(f$sims$gamma[, 2] == 0))
})

test_that("LTZ interval contains the prior band's centre and the sweep runs", {
  f <- apexog(y ~ x1 + x2, endog = d ~ z1, data = fxdat, method = "ltz", mu = 0,
              omega = 0.0025, iterations = 50, seed = 5, sweep_omega = c(1e-6, 0.01),
              verbose = FALSE)
  expect_lt(f$bounds$lower, f$baseline$estimate)
  expect_gt(f$bounds$upper, f$baseline$estimate)
  expect_equal(dim(f$sweep), c(2L, 6L))
  expect_lt(f$sweep$upper[1] - f$sweep$lower[1], f$sweep$upper[2] - f$sweep$lower[2])
})

test_that("no controls, factors in instruments, and methods work", {
  f <- apexog(y ~ 1, endog = d ~ z1, data = fxdat, gmin = -0.1, gmax = 0.1, grid = 3,
              verbose = FALSE)
  expect_true(is.finite(f$baseline$se))
  g <- apexog(y ~ x1, endog = d ~ z1 + factor(cat), data = fxdat, gmin = 0, gmax = 0,
              grid = 2, verbose = FALSE)
  expect_equal(g$instruments, c("z1", "factor(cat)2", "factor(cat)3"))
  expect_equal(coef(f), c(d = f$baseline$estimate))
  expect_equal(dim(vcov(f)), c(1L, 1L))
  expect_equal(colnames(confint(f)), c("lower", "upper"))
  expect_equal(nobs(f), 3000)
  expect_output(print(f), "UCI union")
  expect_output(summary(f), "gamma grid")
  pdf(NULL); on.exit(grDevices::dev.off(), add = TRUE)
  expect_silent(plot(f))
})

test_that("APE SEs with cfinteract equal the delta method with a numerical gradient", {
  # Rebuild the main-equation design by hand, differentiate the APE
  # numerically in the coefficients, and sandwich with the stored
  # Murphy-Topel covariance: this checks every Jacobian entry, including the
  # control-function-product ones (zero in Stata apexog <= 0.20.0).
  numgrad <- function(fn, b) {
    vapply(seq_along(b), function(j) {
      h <- 1e-6 * max(1, abs(b[j]))
      bp <- b; bp[j] <- bp[j] + h
      bm <- b; bm[j] <- bm[j] - h
      (fn(bp) - fn(bm)) / (2 * h)
    }, numeric(1))
  }
  # (i) two continuous endogenous regressors (fixture B6)
  f <- fx_calls$B6()
  W <- cbind(fxdat$z1, fxdat$z2, fxdat$zb, fxdat$x1, 1)
  c1 <- stats::lm.fit(W, fxdat$d)$residuals
  c2 <- stats::lm.fit(W, fxdat$d2)$residuals
  xm <- cbind(fxdat$x1, fxdat$d, fxdat$d2, c1, c2, c1 * c2, 1)
  b <- f$model$coefficients
  expect_equal(names(b), c("x1", "d", "d2", "cf_d", "cf_d2", "cf_d*cf_d2", "(Intercept)"))
  for (t in 1:2) {
    ape <- function(bb) unname(mean(stats::dnorm(drop(xm %*% bb))) * bb[1 + t])
    expect_equal(ape(b), f$baseline$estimate[t], tolerance = 1e-10)
    J <- numgrad(ape, b)
    expect_equal(sqrt(drop(J %*% f$model$vcov %*% J)), f$baseline$se[t], tolerance = 1e-6)
  }
  # (ii) a 0/1 regressor with a probit first stage (discrete-change APE)
  g <- fx_calls$B8()
  W1 <- cbind(fxdat$z1, fxdat$zb, fxdat$x1, 1)
  fs <- stats::glm.fit(W1, fxdat$db, family = stats::binomial(link = "probit"))
  wpi <- fs$linear.predictors
  g1 <- (fxdat$db - stats::pnorm(wpi)) * stats::dnorm(wpi) /
    (stats::pnorm(wpi) * stats::pnorm(-wpi))
  g2 <- stats::lm.fit(cbind(fxdat$z2, fxdat$zb, fxdat$x1, 1), fxdat$d2)$residuals
  xm <- cbind(fxdat$x1, fxdat$db, fxdat$d2, g1, g2, g1 * g2, 1)
  b <- g$model$coefficients
  dc <- function(bb) {
    eta <- drop(xm %*% bb)
    unname(mean(stats::pnorm(eta + bb[2] * (1 - fxdat$db)) - stats::pnorm(eta - bb[2] * fxdat$db)))
  }
  expect_equal(dc(b), g$baseline$estimate[1], tolerance = 1e-8)
  J <- numgrad(dc, b)
  expect_equal(sqrt(drop(J %*% g$model$vcov %*% J)), g$baseline$se[1], tolerance = 1e-6)
})

test_that("an AIE product written x1*x3 in the formula is reused, not duplicated", {
  a <- apexog(yx ~ x1 + x3, endog = d ~ z1, data = fxdat, aie = c("x1", "x3"),
              gmin = 0, gmax = 0, grid = 2, verbose = FALSE)
  b <- apexog(yx ~ x1 * x3, endog = d ~ z1, data = fxdat, aie = c("x1", "x3"),
              gmin = 0, gmax = 0, grid = 2, verbose = FALSE)
  expect_equal(a$baseline$estimate, b$baseline$estimate, tolerance = 1e-10)
  expect_equal(a$baseline$se, b$baseline$se, tolerance = 1e-10)
})
