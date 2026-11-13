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

test_that("effect SEs equal the full two-step influence-function variance, all derivatives numerical", {
  # Stack the first-stage and probit scores into one M-estimator in
  # theta = (b, pi_1, pi_2), take its influence functions psi_i = H^{-1} s_i
  # with H = -dS/dtheta' by central differences, differentiate the effect
  # q(theta) numerically in ALL of theta (probit coefficients AND first-stage
  # parameters), and form sqrt(sum_i (psi_i' J)^2).  This checks the
  # Murphy-Topel channel and the direct first-stage channel together
  # (apexog <= 0.20.1 omitted the latter) without using any analytic
  # derivative of the package.
  numjac <- function(fn, th) {
    f0 <- fn(th)
    out <- matrix(0, length(f0), length(th))
    for (j in seq_along(th)) {
      h <- 1e-6 * max(1, abs(th[j]))
      tp <- th; tp[j] <- tp[j] + h
      tm <- th; tm[j] <- tm[j] - h
      out[, j] <- (fn(tp) - fn(tm)) / (2 * h)
    }
    out
  }
  full_if_se <- function(score_fn, q_fn, theta, cluster = NULL) {
    S <- function(th) colSums(score_fn(th))
    H <- -numjac(S, theta)
    psi <- score_fn(theta) %*% t(solve(H))
    J <- drop(numjac(q_fn, theta))
    l <- drop(psi %*% J)
    if (!is.null(cluster)) l <- rowsum(l, cluster)
    sqrt(sum(l^2))
  }
  probit_psi <- function(eta, y)
    (y - stats::pnorm(eta)) * stats::dnorm(eta) / (stats::pnorm(eta) * stats::pnorm(-eta))

  # (i) two continuous endogenous regressors, linear first stages, cfinteract (B6)
  f <- fx_calls$B6()
  W <- cbind(fxdat$z1, fxdat$z2, fxdat$zb, fxdat$x1, 1)
  pi1 <- unname(stats::lm.fit(W, fxdat$d)$coefficients)
  pi2 <- unname(stats::lm.fit(W, fxdat$d2)$coefficients)
  b <- unname(f$model$coefficients)
  expect_equal(names(f$model$coefficients),
               c("x1", "d", "d2", "cf_d", "cf_d2", "cf_d*cf_d2", "(Intercept)"))
  nb <- length(b); np <- ncol(W)
  design <- function(th) {
    c1 <- fxdat$d - drop(W %*% th[nb + 1:np])
    c2 <- fxdat$d2 - drop(W %*% th[nb + np + 1:np])
    cbind(fxdat$x1, fxdat$d, fxdat$d2, c1, c2, c1 * c2, 1)
  }
  score <- function(th) {
    xm <- design(th)
    cbind(xm * probit_psi(drop(xm %*% th[1:nb]), fxdat$ym2),
          W * (fxdat$d - drop(W %*% th[nb + 1:np])),
          W * (fxdat$d2 - drop(W %*% th[nb + np + 1:np])))
  }
  theta <- c(b, pi1, pi2)
  expect_lt(max(abs(colSums(score(theta)))), 1e-4)
  for (t in 1:2) {
    q <- function(th) mean(stats::dnorm(drop(design(th) %*% th[1:nb]))) * th[1 + t]
    expect_equal(q(theta), f$baseline$estimate[t], tolerance = 1e-10)
    expect_equal(full_if_se(score, q, theta), f$baseline$se[t], tolerance = 1e-5)
  }

  # (ii) a 0/1 regressor with a probit first stage + a continuous one (B8):
  #      generalized-residual control, cfinteract, discrete-change APE
  g <- fx_calls$B8()
  W1 <- cbind(fxdat$z1, fxdat$zb, fxdat$x1, 1)
  W2 <- cbind(fxdat$z2, fxdat$zb, fxdat$x1, 1)
  fs <- stats::glm.fit(W1, fxdat$db, family = stats::binomial(link = "probit"))
  pi1 <- unname(fs$coefficients)
  pi2 <- unname(stats::lm.fit(W2, fxdat$d2)$coefficients)
  b <- unname(g$model$coefficients); nb <- length(b); n1 <- ncol(W1); n2 <- ncol(W2)
  gres <- function(wpi)
    (fxdat$db - stats::pnorm(wpi)) * stats::dnorm(wpi) / (stats::pnorm(wpi) * stats::pnorm(-wpi))
  design <- function(th) {
    g1 <- gres(drop(W1 %*% th[nb + 1:n1]))
    g2 <- fxdat$d2 - drop(W2 %*% th[nb + n1 + 1:n2])
    cbind(fxdat$x1, fxdat$db, fxdat$d2, g1, g2, g1 * g2, 1)
  }
  score <- function(th) {
    xm <- design(th)
    cbind(xm * probit_psi(drop(xm %*% th[1:nb]), fxdat$ym),
          W1 * gres(drop(W1 %*% th[nb + 1:n1])),
          W2 * (fxdat$d2 - drop(W2 %*% th[nb + n1 + 1:n2])))
  }
  theta <- c(b, pi1, pi2)
  expect_lt(max(abs(colSums(score(theta)))), 1e-4)
  q <- function(th) {
    xm <- design(th); eta <- drop(xm %*% th[1:nb])
    mean(stats::pnorm(eta + th[2] * (1 - fxdat$db)) - stats::pnorm(eta - th[2] * fxdat$db))
  }
  expect_equal(q(theta), g$baseline$estimate[1], tolerance = 1e-8)
  expect_equal(full_if_se(score, q, theta), g$baseline$se[1], tolerance = 1e-5)

  # (iii) an exogenous continuous-by-continuous AIE (C1): the first stage
  #       includes the product control, as apexog's does
  a <- fx_calls$C1()
  expect_equal(names(a$model$coefficients),
               c("x1", "x3", "x1:x3", "d", "cf_d", "(Intercept)"))
  W <- cbind(fxdat$z1, fxdat$x1, fxdat$x3, fxdat$x1 * fxdat$x3, 1)
  pi1 <- unname(stats::lm.fit(W, fxdat$d)$coefficients)
  b <- unname(a$model$coefficients); nb <- length(b); np <- ncol(W)
  design <- function(th)
    cbind(fxdat$x1, fxdat$x3, fxdat$x1 * fxdat$x3, fxdat$d,
          fxdat$d - drop(W %*% th[nb + 1:np]), 1)
  score <- function(th) {
    xm <- design(th)
    cbind(xm * probit_psi(drop(xm %*% th[1:nb]), fxdat$yx),
          W * (fxdat$d - drop(W %*% th[nb + 1:np])))
  }
  theta <- c(b, pi1)
  expect_lt(max(abs(colSums(score(theta)))), 1e-4)
  q <- function(th) {
    eta <- drop(design(th) %*% th[1:nb])
    av <- th[1] + th[3] * fxdat$x3
    bv <- th[2] + th[3] * fxdat$x1
    mean(stats::dnorm(eta) * (th[3] - eta * av * bv))
  }
  expect_equal(q(theta), a$baseline$estimate[1], tolerance = 1e-10)
  expect_equal(full_if_se(score, q, theta), a$baseline$se[1], tolerance = 1e-5)

  # (iv) a binary-by-continuous AIE with the continuous interactor endogenous
  #      (C3b): the product column is excluded from the first stage
  cb <- fx_calls$C3b()
  expect_equal(names(cb$model$coefficients),
               c("x1", "x2", "x2:d", "d", "cf_d", "(Intercept)"))
  W <- cbind(fxdat$z1, fxdat$x1, fxdat$x2, 1)
  pi1 <- unname(stats::lm.fit(W, fxdat$d)$coefficients)
  b <- unname(cb$model$coefficients); nb <- length(b); np <- ncol(W)
  design <- function(th)
    cbind(fxdat$x1, fxdat$x2, fxdat$x2 * fxdat$d, fxdat$d,
          fxdat$d - drop(W %*% th[nb + 1:np]), 1)
  score <- function(th) {
    xm <- design(th)
    cbind(xm * probit_psi(drop(xm %*% th[1:nb]), fxdat$yi),
          W * (fxdat$d - drop(W %*% th[nb + 1:np])))
  }
  theta <- c(b, pi1)
  expect_lt(max(abs(colSums(score(theta)))), 1e-4)
  q <- function(th) {           # derivative in d, differenced across the 0/1 x2
    eta <- drop(design(th) %*% th[1:nb])
    pen <- th[2] + th[3] * fxdat$d
    e1 <- eta + (1 - fxdat$x2) * pen
    e0 <- eta - fxdat$x2 * pen
    mean(stats::dnorm(e1) * (th[4] + th[3]) - stats::dnorm(e0) * th[4])
  }
  expect_equal(q(theta), cb$baseline$estimate[1], tolerance = 1e-10)
  expect_equal(full_if_se(score, q, theta), cb$baseline$se[1], tolerance = 1e-5)

  # (v) cluster-robust: the influence functions are summed within clusters (B5)
  b5 <- fx_calls$B5()
  expect_equal(names(b5$model$coefficients), c("x1", "x2", "d", "cf_d", "(Intercept)"))
  W <- cbind(fxdat$z1, fxdat$x1, fxdat$x2, 1)
  pi1 <- unname(stats::lm.fit(W, fxdat$d)$coefficients)
  b <- unname(b5$model$coefficients); nb <- length(b); np <- ncol(W)
  design <- function(th)
    cbind(fxdat$x1, fxdat$x2, fxdat$d, fxdat$d - drop(W %*% th[nb + 1:np]), 1)
  score <- function(th) {
    xm <- design(th)
    cbind(xm * probit_psi(drop(xm %*% th[1:nb]), fxdat$y),
          W * (fxdat$d - drop(W %*% th[nb + 1:np])))
  }
  theta <- c(b, pi1)
  q <- function(th) mean(stats::dnorm(drop(design(th) %*% th[1:nb]))) * th[3]
  expect_equal(q(theta), b5$baseline$estimate[1], tolerance = 1e-10)
  expect_equal(full_if_se(score, q, theta, cluster = fxdat$team), b5$baseline$se[1],
               tolerance = 1e-5)
  expect_gt(abs(full_if_se(score, q, theta) / b5$baseline$se[1] - 1), 1e-3)
})

test_that("fractional-probit and Poisson first stages: effect SEs equal the numerical two-step variance", {
  # Same construction as above (stacked first-stage and probit scores, every
  # derivative by central differences), for the two generalized-residual first
  # stages it did not cover.  For a fractional regressor the derivative of the
  # generalized residual is NOT the binary-probit identity chat * (wpi + chat)
  # (apexog <= 0.21.0 used it for "fprobit"); this test fails with that weight.
  numjac <- function(fn, th) {
    f0 <- fn(th)
    out <- matrix(0, length(f0), length(th))
    for (j in seq_along(th)) {
      h <- 1e-6 * max(1, abs(th[j]))
      tp <- th; tp[j] <- tp[j] + h
      tm <- th; tm[j] <- tm[j] - h
      out[, j] <- (fn(tp) - fn(tm)) / (2 * h)
    }
    out
  }
  full_if_se <- function(score_fn, q_fn, theta) {
    H <- -numjac(function(th) colSums(score_fn(th)), theta)
    psi <- score_fn(theta) %*% t(solve(H))
    sqrt(sum(drop(psi %*% drop(numjac(q_fn, theta)))^2))
  }
  probit_psi <- function(eta, y)
    (y - stats::pnorm(eta)) * stats::dnorm(eta) / (stats::pnorm(eta) * stats::pnorm(-eta))
  W <- cbind(fxdat$z1, fxdat$x1, 1)
  np <- ncol(W)

  # (i) fractional probit (B2)
  f <- fx_calls$B2()
  expect_equal(names(f$model$coefficients), c("x1", "df", "cf_df", "(Intercept)"))
  pi1 <- unname(suppressWarnings(stats::glm.fit(W, fxdat$df,
    family = stats::quasibinomial(link = "probit"),
    control = list(epsilon = 1e-12, maxit = 100))$coefficients))
  b <- unname(f$model$coefficients); nb <- length(b)
  gres <- function(wpi)
    (fxdat$df - stats::pnorm(wpi)) * stats::dnorm(wpi) / (stats::pnorm(wpi) * stats::pnorm(-wpi))
  design <- function(th) cbind(fxdat$x1, fxdat$df, gres(drop(W %*% th[nb + 1:np])), 1)
  score <- function(th) {
    xm <- design(th)
    cbind(xm * probit_psi(drop(xm %*% th[1:nb]), fxdat$yf), W * gres(drop(W %*% th[nb + 1:np])))
  }
  theta <- c(b, pi1)
  expect_lt(max(abs(colSums(score(theta)))), 1e-4)
  q <- function(th) mean(stats::dnorm(drop(design(th) %*% th[1:nb]))) * th[2]
  expect_equal(q(theta), f$baseline$estimate[1], tolerance = 1e-10)
  expect_equal(full_if_se(score, q, theta), f$baseline$se[1], tolerance = 1e-5)
  # the binary-probit weight gives a visibly different first-stage Hessian here
  wpi <- drop(W %*% pi1)
  lam <- stats::dnorm(wpi) / (stats::pnorm(wpi) * stats::pnorm(-wpi))
  A_bin <- gres(wpi) * (wpi + gres(wpi))
  A_frac <- A_bin + fxdat$df * (1 - fxdat$df) * lam^2
  H_num <- -numjac(function(p) colSums(W * gres(drop(W %*% p))), pi1)
  expect_equal(crossprod(W * A_frac, W), H_num, tolerance = 1e-6, ignore_attr = TRUE)
  expect_gt(max(abs(crossprod(W * A_bin, W) / H_num - 1)), 0.1)

  # (ii) Poisson (B3)
  f <- fx_calls$B3()
  expect_equal(names(f$model$coefficients), c("x1", "dc", "cf_dc", "(Intercept)"))
  pi1 <- unname(stats::glm.fit(W, fxdat$dc, family = stats::poisson(),
    control = list(epsilon = 1e-12, maxit = 100))$coefficients)
  b <- unname(f$model$coefficients); nb <- length(b)
  design <- function(th) cbind(fxdat$x1, fxdat$dc, fxdat$dc - exp(drop(W %*% th[nb + 1:np])), 1)
  score <- function(th) {
    xm <- design(th)
    cbind(xm * probit_psi(drop(xm %*% th[1:nb]), fxdat$yc),
          W * (fxdat$dc - exp(drop(W %*% th[nb + 1:np]))))
  }
  theta <- c(b, pi1)
  expect_lt(max(abs(colSums(score(theta)))), 1e-4)
  q <- function(th) mean(stats::dnorm(drop(design(th) %*% th[1:nb]))) * th[2]
  expect_equal(q(theta), f$baseline$estimate[1], tolerance = 1e-10)
  expect_equal(full_if_se(score, q, theta), f$baseline$se[1], tolerance = 1e-5)
})

test_that("a fractional-probit first stage on a 0/1 regressor reproduces the probit first stage", {
  a <- fx_calls$B1()
  b <- apexog(yb ~ x1 + x2, endog = db ~ z1 + zb, data = fxdat, firststage = "fprobit",
              gmin = -0.1, gmax = 0.1, grid = 3, verbose = FALSE)
  expect_equal(b$baseline$estimate, a$baseline$estimate, tolerance = 1e-8)
  expect_equal(b$baseline$se, a$baseline$se, tolerance = 1e-8)
})

test_that("an AIE product written x1*x3 in the formula is reused, not duplicated", {
  a <- apexog(yx ~ x1 + x3, endog = d ~ z1, data = fxdat, aie = c("x1", "x3"),
              gmin = 0, gmax = 0, grid = 2, verbose = FALSE)
  b <- apexog(yx ~ x1 * x3, endog = d ~ z1, data = fxdat, aie = c("x1", "x3"),
              gmin = 0, gmax = 0, grid = 2, verbose = FALSE)
  expect_equal(a$baseline$estimate, b$baseline$estimate, tolerance = 1e-10)
  expect_equal(a$baseline$se, b$baseline$se, tolerance = 1e-10)
})
