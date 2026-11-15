# apexog 0.22.0: the AIE target with firststage = "biprobit".  The estimates
# and SEs are compared with Stata in test-stata-uci.R / test-stata-ltz.R (specs
# D3-D5, E7); here the SE is checked against a delta method whose gradient is
# taken numerically from an AIE written out from scratch (no package formula),
# and the model and the error paths are checked.

numgrad <- function(fn, th) {
  vapply(seq_along(th), function(j) {
    h <- 1e-6 * max(1, abs(th[j]))
    tp <- th; tp[j] <- tp[j] + h
    tm <- th; tm[j] <- tm[j] - h
    (fn(tp) - fn(tm)) / (2 * h)
  }, numeric(1))
}

test_that("binary moderator: the SE is the delta method with a numerical gradient", {
  f <- apexog(ybb ~ x1 + x2, endog = db ~ z1 + zb, data = fxdat, firststage = "biprobit",
              aie = c("db", "x2"), gmin = 0, gmax = 0, grid = 2, verbose = FALSE)
  nm <- names(f$model$coefficients)
  expect_equal(nm[1:5], c("ybb:x1", "ybb:x2", "ybb:db:x2", "ybb:db", "ybb:(Intercept)"))
  # the treatment equation excludes the product
  expect_equal(nm[6:10], c("db:z1", "db:zb", "db:x1", "db:x2", "db:(Intercept)"))
  expect_equal(nm[11], "athrho")
  th <- unname(f$model$coefficients)
  x1 <- fxdat$x1
  q <- function(t) {
    eta <- function(dv, mv) t[1] * x1 + t[2] * mv + t[3] * dv * mv + t[4] * dv + t[5]
    mean(stats::pnorm(eta(1, 1)) - stats::pnorm(eta(1, 0)) -
           stats::pnorm(eta(0, 1)) + stats::pnorm(eta(0, 0)))
  }
  expect_equal(q(th), f$baseline$estimate, tolerance = 1e-12)
  J <- numgrad(q, th)
  expect_equal(sqrt(drop(J %*% f$model$vcov %*% J)), f$baseline$se, tolerance = 1e-6)
})

test_that("continuous moderator: the SE is the delta method with a numerical gradient", {
  f <- apexog(ybc ~ x1 + x3, endog = db ~ z1 + zb, data = fxdat, firststage = "biprobit",
              aie = c("db", "x3"), gmin = 0, gmax = 0, grid = 2, verbose = FALSE)
  expect_equal(names(f$model$coefficients)[1:5],
               c("ybc:x1", "ybc:x3", "ybc:db:x3", "ybc:db", "ybc:(Intercept)"))
  th <- unname(f$model$coefficients)
  x1 <- fxdat$x1; x3 <- fxdat$x3
  q <- function(t) {
    # the x3-slope of the ASF, averaged, at db = 1 minus at db = 0
    sl <- function(dv) {
      eta <- t[1] * x1 + t[2] * x3 + t[3] * dv * x3 + t[4] * dv + t[5]
      mean(stats::dnorm(eta) * (t[2] + t[3] * dv))
    }
    sl(1) - sl(0)
  }
  expect_equal(q(th), f$baseline$estimate, tolerance = 1e-12)
  J <- numgrad(q, th)
  expect_equal(sqrt(drop(J %*% f$model$vcov %*% J)), f$baseline$se, tolerance = 1e-6)
  # both target orders give the same number
  g <- apexog(ybc ~ x1 + x3, endog = db ~ z1 + zb, data = fxdat, firststage = "biprobit",
              aie = c("x3", "db"), gmin = 0, gmax = 0, grid = 2, verbose = FALSE)
  expect_equal(g$baseline$estimate, f$baseline$estimate, tolerance = 1e-10)
  expect_equal(g$baseline$se, f$baseline$se, tolerance = 1e-10)
})

test_that("the bivariate-probit AIE rejects what it does not support", {
  q <- function(...) apexog(ybb ~ x1 + x2, data = fxdat, firststage = "biprobit",
                            gmin = 0, gmax = 0, grid = 2, verbose = FALSE, ...)
  expect_error(q(endog = db + db:x2 ~ z1 + z1:x2, aie = c("db", "x2")),
               "needs no first stage of its own")
  expect_error(q(endog = db ~ z1 + zb, aie = c("x1", "x2")),
               "must interact the 0/1 treatment")
  expect_error(q(endog = db ~ z1 + zb, aie = c("db", "x2"), aie_binary = c(FALSE, TRUE)),
               "enters the interaction as a 0/1 variable")
})

test_that("control-function binary-by-binary AIE: SE = two-step influence function, all derivatives numerical", {
  # 0.22.0 fixed the Jacobian entries of the two interactors and their product
  # in the binary-by-binary AIE (they used j - X_i etc. instead of the
  # counterfactual values j, k, jk).  This stacks the linear first stage and
  # the probit into one M-estimator and differentiates everything numerically.
  f <- apexog(ybb ~ x1 + x2, endog = db ~ z1 + zb, data = fxdat, aie = c("db", "x2"),
              gmin = 0, gmax = 0, grid = 2, verbose = FALSE)
  expect_equal(names(f$model$coefficients),
               c("x1", "x2", "db:x2", "db", "cf_db", "(Intercept)"))
  W <- cbind(fxdat$z1, fxdat$zb, fxdat$x1, fxdat$x2, 1)
  pi1 <- unname(stats::lm.fit(W, fxdat$db)$coefficients)
  b <- unname(f$model$coefficients)
  nb <- length(b); np <- ncol(W)
  x1 <- fxdat$x1; x2 <- fxdat$x2; db <- fxdat$db; y <- fxdat$ybb
  cf <- function(th) db - drop(W %*% th[nb + 1:np])
  score <- function(th) {
    xm <- cbind(x1, x2, db * x2, db, cf(th), 1)
    eta <- drop(xm %*% th[1:nb])
    psi <- (y - stats::pnorm(eta)) * stats::dnorm(eta) /
      (stats::pnorm(eta) * stats::pnorm(-eta))
    cbind(xm * psi, W * cf(th))
  }
  theta <- c(b, pi1)
  expect_lt(max(abs(colSums(score(theta)))), 1e-4)
  q <- function(th) {
    t <- th[1:nb]; c0 <- cf(th)
    eta <- function(dv, mv) t[1] * x1 + t[2] * mv + t[3] * dv * mv + t[4] * dv +
      t[5] * c0 + t[6]
    mean(stats::pnorm(eta(1, 1)) - stats::pnorm(eta(1, 0)) -
           stats::pnorm(eta(0, 1)) + stats::pnorm(eta(0, 0)))
  }
  expect_equal(q(theta), f$baseline$estimate, tolerance = 1e-10)
  H <- -sapply(seq_along(theta), function(j) {
    h <- 1e-6 * max(1, abs(theta[j]))
    tp <- theta; tp[j] <- tp[j] + h
    tm <- theta; tm[j] <- tm[j] - h
    (colSums(score(tp)) - colSums(score(tm))) / (2 * h)
  })
  psi <- score(theta) %*% t(solve(H))
  J <- numgrad(q, theta)
  expect_equal(sqrt(sum(drop(psi %*% J)^2)), f$baseline$se, tolerance = 1e-5)
})
