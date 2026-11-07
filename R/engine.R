# Control-function engine: one fit of the main probit at a given gamma-offset,
# the Murphy-Topel two-step covariance, and the APE / AIE delta method.
# A line-by-line port of _apexogfit + Mata _apexogfit_apese_mt (Stata apexog
# v0.20.1).  `S` is the model setup built once by apexog().

# --- first stages ------------------------------------------------------------
# Returns list(chat, A): the control function (residual / generalized residual)
# and the per-observation Murphy-Topel weight A_i (Mata comment block).
.first_stage <- function(d, W, type, name) {
  if (type == "linear") {
    f <- stats::lm.fit(W, d)
    wpi <- d - f$residuals
    return(list(chat = f$residuals, A = rep(1, length(d)), wpi = wpi))
  }
  fam <- switch(type,
    probit  = stats::binomial(link = "probit"),
    fprobit = stats::quasibinomial(link = "probit"),
    poisson = stats::poisson()
  )
  f <- withCallingHandlers(
    stats::glm.fit(W, d, family = fam,
                   control = list(epsilon = 1e-12, maxit = 100, trace = FALSE)),
    warning = function(w) {
      if (grepl("fitted probabilities numerically 0 or 1", conditionMessage(w)))
        stop("The ", type, " first stage for ", name,
             " shows (quasi-)perfect prediction; apexog cannot proceed.", call. = FALSE)
      if (grepl("non-integer", conditionMessage(w))) invokeRestart("muffleWarning")
    })
  if (!isTRUE(f$converged))
    stop("The ", type, " first stage for ", name, " did not converge.", call. = FALSE)
  wpi <- f$linear.predictors
  if (type %in% c("probit", "fprobit")) {
    Phi <- stats::pnorm(wpi)
    Phim <- stats::pnorm(-wpi)
    chat <- (d - Phi) * stats::dnorm(wpi) / (Phi * Phim)
    A <- chat * (wpi + chat)
  } else {
    mu <- exp(wpi)
    chat <- d - mu
    A <- mu
  }
  if (any(!is.finite(chat)))
    stop("The control function for ", name, " is undefined in ",
         sum(!is.finite(chat)), " observation(s) (extreme first-stage fitted values).",
         call. = FALSE)
  list(chat = chat, A = A, wpi = wpi)
}

# --- main probit at an offset --------------------------------------------------
.probit_offset <- function(x, y, off, start) {
  fit <- tryCatch(
    suppressWarnings(stats::glm.fit(x, y, offset = off, start = start,
                                    family = stats::binomial(link = "probit"),
                                    control = list(epsilon = 1e-12, maxit = 100,
                                                   trace = FALSE))),
    error = function(e) NULL)
  if (is.null(fit) || !isTRUE(fit$converged)) {
    # retry from zero if a warm start failed
    fit <- tryCatch(
      suppressWarnings(stats::glm.fit(x, y, offset = off,
                                      family = stats::binomial(link = "probit"),
                                      control = list(epsilon = 1e-12, maxit = 200,
                                                     trace = FALSE))),
      error = function(e) NULL)
  }
  fit
}

.na_fit <- function(J, K, msg) {
  list(ok = FALSE, msg = msg, est = rep(NA_real_, J), se = rep(NA_real_, J),
       lo = rep(NA_real_, J), hi = rep(NA_real_, J),
       vtgt = matrix(NA_real_, J, J), mphi = NA_real_, stilde = NA_real_,
       gpr = rep(NA_real_, K), coef = NULL)
}

# Exact probability-scale value of each instrument's violation at a refit
# (_apexoggpr): 0/1 instrument -> average discrete change, else average
# derivative.  `eta` is the fitted index INCLUDING the offset.
.gpr <- function(eta, Z, gam, zbin) {
  K <- ncol(Z)
  out <- numeric(K)
  for (k in seq_len(K)) {
    g <- gam[k]
    if (zbin[k]) {
      base <- eta - g * Z[, k]
      out[k] <- mean(stats::pnorm(base + g) - stats::pnorm(base))
    } else {
      out[k] <- mean(stats::dnorm(eta)) * g
    }
  }
  out
}

# --- one control-function fit --------------------------------------------------
.fit_cf <- function(S, gam, start = NULL) {
  J <- if (is.null(S$aie)) ncol(S$D) else 1L
  K <- ncol(S$Z)
  off <- drop(S$Z %*% gam)
  xm <- S$xmain
  nx <- ncol(xm)
  fit <- .probit_offset(xm, S$y, off, start)
  if (is.null(fit)) return(.na_fit(J, K, "probit failed"))
  if (!isTRUE(fit$converged)) return(.na_fit(J, K, "probit did not converge"))
  b <- fit$coefficients
  aliased <- is.na(b)
  if (any(aliased[S$load_cols]))
    return(.na_fit(J, K, paste0(colnames(xm)[S$load_cols][aliased[S$load_cols]][1],
                                " omitted from the probit (collinear)")))
  b[aliased] <- 0

  n <- length(S$y)
  eta <- drop(xm %*% b) + off
  phi <- stats::dnorm(eta)
  Phi <- stats::pnorm(eta)
  Phim <- stats::pnorm(-eta)
  mphi <- mean(phi)
  Jn <- ncol(S$D)
  bcf <- b[S$chat_idx]
  stilde <- 1 / sqrt(1 + drop(t(bcf) %*% stats::var(S$CHAT) %*% bcf))
  psi <- (S$y - Phi) * phi / (Phi * Phim)
  if (any(!is.finite(psi))) return(.na_fit(J, K, "probit score not finite"))
  pep <- psi * (eta + psi)

  # probit Hessian (observed information) and its g-inverse
  Hb <- crossprod(xm * pep, xm)
  U2 <- .ginv_sym(Hb, !aliased)

  # first-stage contributions (Murphy-Topel)
  stot <- xm * psi
  for (k in seq_len(Jn)) {
    Wk <- S$W[, S$wcols[[k]], drop = FALSE]
    chat_k <- S$CHAT[, k]
    A_k <- S$A[, k]
    wd_k <- -A_k
    M_k <- .ginv_sym(crossprod(Wk * A_k, Wk))
    psi_pi <- (Wk * chat_k) %*% t(M_k)
    g_k <- rep(b[S$chat_idx[k]], n)
    partners <- which(S$cf_a == k | S$cf_b == k)
    for (cc in partners) {
      m <- if (S$cf_a[cc] == k) S$cf_b[cc] else S$cf_a[cc]
      g_k <- g_k + b[S$cfp_idx[cc]] * S$CHAT[, m]
    }
    A1 <- crossprod(xm * (pep * g_k * wd_k), Wk)
    A2 <- crossprod(Wk, psi * wd_k)
    J2 <- -A1
    J2[S$chat_idx[k], ] <- J2[S$chat_idx[k], ] + drop(A2)
    for (cc in partners) {
      m <- if (S$cf_a[cc] == k) S$cf_b[cc] else S$cf_a[cc]
      J2[S$cfp_idx[cc], ] <- J2[S$cfp_idx[cc], ] +
        drop(crossprod(Wk, psi * S$CHAT[, m] * wd_k))
    }
    stot <- stot + psi_pi %*% t(J2)
  }
  meat <- if (is.null(S$cluster)) crossprod(stot) else {
    Sc <- rowsum(stot, S$cluster, reorder = FALSE)
    crossprod(Sc)
  }
  Vb <- U2 %*% meat %*% U2
  z <- stats::qnorm(1 - (1 - S$level) / 2)

  if (!is.null(S$aie)) {
    r <- .aie_target(S, xm, b, eta, phi)
    se2 <- drop(r$J %*% Vb %*% r$J)
    se <- if (se2 < 0) NA_real_ else sqrt(se2)
    out <- list(ok = TRUE, msg = "", est = r$est, se = se,
                lo = r$est - z * se, hi = r$est + z * se,
                vtgt = matrix(se^2, 1, 1), mphi = mphi, stilde = stilde,
                coef = stats::setNames(b, colnames(xm)), Vbeta = Vb)
  } else {
    est <- se <- numeric(Jn)
    Jall <- matrix(0, nx, Jn)
    p <- S$p
    for (t in seq_len(Jn)) {
      di <- S$d_idx[t]
      bd <- b[di]
      Jb <- numeric(nx)
      Dt <- S$D[, t]
      if (S$dbin[t]) {
        e1 <- eta + bd * (1 - Dt)
        e0 <- eta - bd * Dt
        f1 <- stats::dnorm(e1)
        f0 <- stats::dnorm(e0)
        dphi <- f1 - f0
        est[t] <- mean(stats::pnorm(e1) - stats::pnorm(e0))
        Jb <- colMeans(xm * dphi)
        Jb[di] <- mean(f1)
      } else {
        est[t] <- mphi * bd
        Jb <- bd * colMeans(xm * (-eta * phi))
        Jb[di] <- Jb[di] + mphi
      }
      # Jb covers every column of xmain, including the cfinteract products
      # (Stata apexog <= 0.20.0 left those entries at 0; fixed in 0.20.1).
      Jall[, t] <- Jb
      se2 <- drop(Jb %*% Vb %*% Jb)
      se[t] <- if (se2 < 0) NA_real_ else sqrt(se2)
    }
    Vt <- t(Jall) %*% Vb %*% Jall
    Vt <- (Vt + t(Vt)) / 2
    out <- list(ok = TRUE, msg = "", est = est, se = se,
                lo = est - z * se, hi = est + z * se, vtgt = Vt,
                mphi = mphi, stilde = stilde,
                coef = stats::setNames(b, colnames(xm)), Vbeta = Vb)
  }
  out$gpr <- .gpr(eta, S$Z, gam, S$zbin)
  out
}

# --- AIE target (Mata AIE branch; four cases) ---------------------------------
.aie_target <- function(S, xm, b, eta, phi) {
  a <- S$aie
  ix <- a$x_idx; iz <- a$z_idx; ixz <- a$xz_idx
  X <- xm[, ix]; Zv <- xm[, iz]; XZ <- xm[, ixz]
  bX <- b[ix]; bZ <- b[iz]; bXZ <- b[ixz]
  nx <- ncol(xm)
  if (a$case == 0L) {
    av <- bX + bXZ * Zv
    bv <- bZ + bXZ * X
    w <- bXZ - eta * av * bv
    est <- mean(phi * w)
    Tm <- -eta * w - av * bv
    Jb <- colMeans(xm * (phi * Tm))
    Jb[ix] <- Jb[ix] - mean(phi * eta * bv)
    Jb[iz] <- Jb[iz] - mean(phi * eta * av)
    Jb[ixz] <- Jb[ixz] + mean(phi) - mean(phi * eta * Zv * bv) -
      mean(phi * eta * av * X)
  } else if (a$case %in% c(1L, 2L)) {
    if (a$case == 1L) {           # X cts, Z binary
      cts_i <- ix; bin_i <- iz; C <- X; B <- Zv; bC <- bX; bB <- bZ
    } else {                      # X binary, Z cts
      cts_i <- iz; bin_i <- ix; C <- Zv; B <- X; bC <- bZ; bB <- bX
    }
    pen <- bB + bXZ * C
    e1 <- eta + (1 - B) * pen
    e0 <- eta - B * pen
    f1 <- stats::dnorm(e1)
    f0 <- stats::dnorm(e0)
    est <- mean(f1 * (bC + bXZ) - f0 * bC)
    G <- -e1 * f1 * (bC + bXZ) + e0 * f0 * bC
    Jb <- colMeans(xm * G)
    Jb[cts_i] <- Jb[cts_i] + mean(f1) - mean(f0)
    Jb[bin_i] <- mean(-e1 * f1 * (bC + bXZ))
    Jb[ixz] <- mean(-e1 * f1 * (bC + bXZ) * C) + mean(f1)
  } else {
    e11 <- eta + (1 - X) * bX + (1 - Zv) * bZ + (1 - X * Zv) * bXZ
    e10 <- eta + (1 - X) * bX + (0 - Zv) * bZ + (0 - X * Zv) * bXZ
    e01 <- eta + (0 - X) * bX + (1 - Zv) * bZ + (0 - X * Zv) * bXZ
    e00 <- eta + (0 - X) * bX + (0 - Zv) * bZ + (0 - X * Zv) * bXZ
    est <- mean(stats::pnorm(e11) - stats::pnorm(e10) - stats::pnorm(e01) +
                  stats::pnorm(e00))
    f11 <- stats::dnorm(e11); f10 <- stats::dnorm(e10)
    f01 <- stats::dnorm(e01); f00 <- stats::dnorm(e00)
    Dphi <- f11 - f10 - f01 + f00
    Jb <- colMeans(xm * Dphi)
    Jb[ix] <- mean((1 - X) * (f11 - f10)) - mean(X * (f01 - f00))
    Jb[iz] <- mean((1 - Zv) * (f11 - f01)) - mean(Zv * (f10 - f00))
    Jb[ixz] <- mean(f11) - mean(XZ * Dphi)
  }
  list(est = est, J = Jb)
}
