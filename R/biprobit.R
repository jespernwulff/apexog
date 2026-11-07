# Recursive bivariate probit engine (firststage = "biprobit"): the port of
# apexog's v0.19 engine.  At each gamma it fits
#     y = 1{x1'b1 + gamma'z + e1 > 0},  x1 = [X, d, 1]
#     d = 1{x2'b2 + e2 > 0},            x2 = [Z, X, 1]
#     corr(e1, e2) = rho = tanh(athrho)
# jointly by maximum likelihood (Stata's -biprobit- with offset() in the outcome
# equation), then computes the ASF discrete-change APE of d with a delta-method
# SE from the robust (or cluster-robust) MLE covariance.  Parameter order
# matches Stata's e(b): [eq1: X, d, _cons | eq2: Z, X, _cons | athrho].

# Log likelihood, per-observation scores and (optionally) the analytic
# negative Hessian.  Notation: w1 = q1*xb1, w2 = q2*xb2, r = q1*q2*rho,
# P = Phi2(w1, w2, r); derivatives of Phi2: P_w1 = phi(w1)Phi(v1),
# P_w2 = phi(w2)Phi(v2), P_r = phi2, P_w1w1 = -w1 P_w1 - r phi2,
# P_w2w2 = -w2 P_w2 - r phi2, P_w1w2 = phi2, P_w1r = phi2 (r w2 - w1)/(1-r^2),
# P_w2r = phi2 (r w1 - w2)/(1-r^2),
# P_rr = phi2 [r/(1-r^2) + (w1 w2 (1-r^2) - r Q)/(1-r^2)^2],
# Q = w1^2 - 2 r w1 w2 + w2^2.
.bp_derivs <- function(theta, S, off, hess = FALSE) {
  p1 <- ncol(S$bp_x1); p2 <- ncol(S$bp_x2)
  b1 <- theta[seq_len(p1)]
  b2 <- theta[p1 + seq_len(p2)]
  ath <- theta[p1 + p2 + 1L]
  rho <- tanh(ath)
  xb1 <- drop(S$bp_x1 %*% b1) + off
  xb2 <- drop(S$bp_x2 %*% b2)
  q1 <- 2 * S$y - 1
  q2 <- 2 * S$D[, 1] - 1
  w1 <- q1 * xb1
  w2 <- q2 * xb2
  r <- q1 * q2 * rho
  P <- pbivnorm::pbivnorm(w1, w2, r)
  P <- pmax(P, 1e-300)
  s2 <- 1 - r^2
  s <- sqrt(s2)
  Pw1 <- stats::dnorm(w1) * stats::pnorm((w2 - r * w1) / s)
  Pw2 <- stats::dnorm(w2) * stats::pnorm((w1 - r * w2) / s)
  Q <- w1^2 - 2 * r * w1 * w2 + w2^2
  phi2 <- exp(-Q / (2 * s2)) / (2 * pi * s)
  drho <- 1 - rho^2                      # d rho / d athrho
  l1 <- q1 * Pw1 / P                     # dL/dxb1
  l2 <- q2 * Pw2 / P                     # dL/dxb2
  lr <- q1 * q2 * phi2 / P               # dL/drho
  out <- list(ll = sum(log(P)),
              sc = cbind(S$bp_x1 * l1, S$bp_x2 * l2, lr * drho))
  if (hess) {
    Pw1w1 <- -w1 * Pw1 - r * phi2
    Pw2w2 <- -w2 * Pw2 - r * phi2
    Pw1r <- phi2 * (r * w2 - w1) / s2
    Pw2r <- phi2 * (r * w1 - w2) / s2
    Prr <- phi2 * (r / s2 + (w1 * w2 * s2 - r * Q) / s2^2)
    a1 <- Pw1 / P; a2 <- Pw2 / P; ar <- phi2 / P
    h11 <- Pw1w1 / P - a1^2
    h22 <- Pw2w2 / P - a2^2
    h12 <- q1 * q2 * (phi2 / P - a1 * a2)
    h1r <- q2 * (Pw1r / P - a1 * ar)
    h2r <- q1 * (Pw2r / P - a2 * ar)
    hrr <- Prr / P - ar^2
    H <- matrix(0, p1 + p2 + 1L, p1 + p2 + 1L)
    i1 <- seq_len(p1); i2 <- p1 + seq_len(p2); ir <- p1 + p2 + 1L
    H[i1, i1] <- crossprod(S$bp_x1 * h11, S$bp_x1)
    H[i2, i2] <- crossprod(S$bp_x2 * h22, S$bp_x2)
    H[i1, i2] <- crossprod(S$bp_x1 * h12, S$bp_x2)
    H[i2, i1] <- t(H[i1, i2])
    H[i1, ir] <- H[ir, i1] <- colSums(S$bp_x1 * h1r) * drho
    H[i2, ir] <- H[ir, i2] <- colSums(S$bp_x2 * h2r) * drho
    H[ir, ir] <- sum(hrr) * drho^2 + sum(lr) * (-2 * rho * drho)
    out$negH <- -H
  }
  out
}

.bp_grad <- function(theta, S, off) colSums(.bp_derivs(theta, S, off)$sc)

# Negative Hessian by central differences of the analytic gradient (used only
# to check the analytic Hessian in the test suite).
.bp_neghess_num <- function(theta, S, off) {
  k <- length(theta)
  H <- matrix(0, k, k)
  for (j in seq_len(k)) {
    h <- 1e-5 * max(1, abs(theta[j]))
    tp <- theta; tp[j] <- tp[j] + h
    tm <- theta; tm[j] <- tm[j] - h
    H[, j] <- (.bp_grad(tp, S, off) - .bp_grad(tm, S, off)) / (2 * h)
  }
  -(H + t(H)) / 2
}

.bp_start <- function(S, off) {
  f1 <- suppressWarnings(stats::glm.fit(S$bp_x1, S$y, offset = off,
                                        family = stats::binomial(link = "probit")))
  f2 <- suppressWarnings(stats::glm.fit(S$bp_x2, S$D[, 1],
                                        family = stats::binomial(link = "probit")))
  b1 <- f1$coefficients; b1[is.na(b1)] <- 0
  b2 <- f2$coefficients; b2[is.na(b2)] <- 0
  c(b1, b2, 0)
}

# Damped Newton-Raphson with the analytic score and Hessian.
.bp_mle <- function(S, off, start = NULL, maxit = 100) {
  theta <- if (is.null(start)) .bp_start(S, off) else start
  cur <- .bp_derivs(theta, S, off, hess = TRUE)
  if (!is.finite(cur$ll)) {
    theta <- .bp_start(S, off)
    cur <- .bp_derivs(theta, S, off, hess = TRUE)
  }
  conv <- FALSE
  for (it in seq_len(maxit)) {
    g <- colSums(cur$sc)
    step <- tryCatch(solve(cur$negH, g), error = function(e) NULL)
    if (is.null(step) || any(!is.finite(step))) break
    crit <- sum(g * step)
    if (is.finite(crit) && abs(crit) < 1e-16) { conv <- TRUE; break }
    lam <- 1
    repeat {
      cand <- theta + lam * step
      nxt <- .bp_derivs(cand, S, off, hess = TRUE)
      if (is.finite(nxt$ll) && nxt$ll >= cur$ll - 1e-10) break
      lam <- lam / 2
      if (lam < 1e-8) break
    }
    if (lam < 1e-8) break
    theta <- cand
    cur <- nxt
    if (max(abs(lam * step)) < 1e-12) { conv <- TRUE; break }
  }
  if (!conv) {
    g <- colSums(cur$sc)
    crit <- tryCatch(sum(g * solve(cur$negH, g)), error = function(e) Inf)
    conv <- is.finite(crit) && abs(crit) < 1e-10
  }
  list(theta = theta, converged = conv, sc = cur$sc, ll = cur$ll, negH = cur$negH)
}

.fit_bp <- function(S, gam, start = NULL) {
  K <- ncol(S$Z)
  off <- drop(S$Z %*% gam)
  m <- tryCatch(.bp_mle(S, off, start), error = function(e) NULL)
  if (is.null(m) || !m$converged) return(.na_fit(1L, K, "biprobit did not converge"))
  theta <- m$theta
  Hi <- tryCatch(solve(m$negH), error = function(e) NULL)
  if (is.null(Hi)) return(.na_fit(1L, K, "biprobit Hessian singular"))
  n <- length(S$y)
  if (is.null(S$cluster)) {
    V <- Hi %*% crossprod(m$sc) %*% Hi * (n / (n - 1))
  } else {
    Sc <- rowsum(m$sc, S$cluster, reorder = FALSE)
    G <- nrow(Sc)
    V <- Hi %*% crossprod(Sc) %*% Hi * (G / (G - 1))
  }
  p1 <- ncol(S$bp_x1)
  b1 <- theta[seq_len(p1)]
  di <- p1 - 1L                       # d is the column before _cons
  d <- S$D[, 1]
  eta <- drop(S$bp_x1 %*% b1) + off
  e1 <- eta + b1[di] * (1 - d)
  e0 <- eta - b1[di] * d
  f1 <- stats::dnorm(e1)
  f0 <- stats::dnorm(e0)
  dphi <- f1 - f0
  est <- mean(stats::pnorm(e1) - stats::pnorm(e0))
  Jb <- numeric(length(theta))
  Jb[seq_len(p1)] <- colMeans(S$bp_x1 * dphi)
  Jb[di] <- mean(f1)
  se2 <- drop(Jb %*% V %*% Jb)
  se <- if (se2 < 0) NA_real_ else sqrt(se2)
  z <- stats::qnorm(1 - (1 - S$level) / 2)
  list(ok = TRUE, msg = "", est = est, se = se, lo = est - z * se, hi = est + z * se,
       vtgt = matrix(se^2, 1, 1), mphi = mean(stats::dnorm(eta)), stilde = 1,
       gpr = .gpr(eta, S$Z, gam, S$zbin), coef = stats::setNames(theta, S$bp_names), V = V,
       rho = tanh(theta[length(theta)]))
}
