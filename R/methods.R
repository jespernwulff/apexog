.fmt <- function(x, d = 6) formatC(x, format = "f", digits = d)

.tlabel <- function(x) {
  if (x$target == "aie") paste0("AIE(", paste(x$aie$vars, collapse = " x "), ")")
  else paste0("APE(", x$baseline$target, ")")
}

#' Print and extract an apexog fit
#'
#' \code{print} shows the results display: the baseline estimate at
#' \eqn{\gamma = 0}, the scale factors, and the UCI union or the LTZ interval
#' with the prior-only band.  \code{coef} and \code{vcov} return the baseline
#' target estimates and their (Murphy-Topel or MLE) covariance; \code{nobs}
#' the number of observations used.
#' @param x,object An \code{"apexog"} object.
#' @param ... Unused.
#' @return \code{print} returns \code{x} invisibly; \code{coef} a named
#'   vector; \code{vcov} a matrix; \code{nobs} an integer.
#' @examples
#' data(promotion)
#' fit <- apexog(promoted ~ tenure + manager, endog = training ~ subsidy,
#'               data = promotion, gmin = -0.1, gmax = 0.1, grid = 3,
#'               verbose = FALSE)
#' fit
#' coef(fit)
#' sqrt(diag(vcov(fit)))
#' nobs(fit)
#' @name apexog-methods
#' @export
print.apexog <- function(x, ...) {
  hr <- strrep("-", 78)
  what <- if (x$target == "aie") "AIE" else "APE"
  meth <- if (x$method == "uci") "UCI method" else
    paste0("LTZ method, ", x$prior$type, " prior")
  cat(hr, "\n", sep = "")
  cat("Plausibly-exogenous ", what, " bounds (CHR 2012, ", meth, ")\n", sep = "")
  if (x$estimator == "biprobit") {
    cat("Estimator: biprobit (recursive bivariate probit, MLE)  |  SE: delta method  |  N = ",
        x$N, "\n", sep = "")
  } else {
    cat("Estimator: cfprobit  |  First stage: ", paste(x$firststage, collapse = " "),
        "  |  SE: Murphy-Topel  |  N = ", x$N, "\n", sep = "")
    if (any(x$firststage %in% c("probit", "fprobit", "poisson")))
      cat("Note:      the generalized-residual control function (probit/fprobit/poisson\n",
          "           first stage) is consistent under Wooldridge's sufficiency assumption\n",
          "           and only approximate under joint normality of the structural errors;\n",
          "           for a 0/1 treatment firststage = \"biprobit\" is exact.\n", sep = "")
  }
  if (!is.null(x$N_clust)) cat("Clusters:  ", x$N_clust, " (", x$clustvar, ")\n", sep = "")
  if (x$target == "aie") {
    cat("Target:    Average Interaction Effect (", paste(x$aie$vars, collapse = " x "),
        "): ", x$aie$label, "\n", sep = "")
  } else {
    cat("Endogenous: ", paste(x$endog, collapse = " "), "\n", sep = "")
    cat("Target:    APE of ", paste(x$targets, collapse = " "), "\n", sep = "")
  }
  cat(hr, "\n", sep = "")
  cat("Baseline (gamma = 0):\n")
  b <- x$baseline
  lab <- .tlabel(x)
  show <- if (x$target == "aie") seq_len(nrow(b)) else which(b$target %in% x$targets)
  for (i in show)
    cat("    ", format(lab[i], width = 18), " = ", .fmt(b$estimate[i]), "  (SE ",
        .fmt(b$se[i]), ")  [", .fmt(b$lower[i]), ", ", .fmt(b$upper[i]), "]\n", sep = "")
  cat("\nScale:     1 latent-index unit of gamma = ", .fmt(x$mphi, 4),
      " pr. points per unit of Z (avg., baseline)\n", sep = "")
  cat("           1 structural-latent unit of gamma = ", .fmt(1 / x$stilde, 4),
      " latent-index units  (sigma_tilde = ", .fmt(x$stilde, 4), ")\n", sep = "")
  if (x$estimator == "biprobit")
    cat("           biprobit: structural and index scales coincide (sigma_tilde = 1);\n",
        "           rho(eps_y, eps_d) at the baseline fit = ", .fmt(x$rho, 4), "\n", sep = "")
  if (x$gunits == "pr")
    cat("           gunits = \"pr\": gamma values were probability points, divided by m_phi\n",
        "           (first order); values below are latent units, exact Pr scale alongside.\n",
        sep = "")
  if (x$gunits == "structural") {
    cat("           gunits = \"structural\": gamma values were divided by sigma_tilde;\n",
        "           values below are latent-index units.\n", sep = "")
    if (any(x$firststage %in% c("probit", "fprobit", "poisson")))
      cat("Warning:   sigma_tilde conversion derived for a LINEAR first stage.\n")
  }
  lvl <- format(100 * x$level)
  inst <- x$instruments
  if (x$method == "uci") {
    cat("\nUCI union over gamma grid (", x$gridn, " point(s) per instrument):\n", sep = "")
    G <- as.matrix(x$grid[, seq_along(inst), drop = FALSE])
    for (k in seq_along(inst)) {
      pr <- x$grid_pr[, k]
      cat("    gamma_", format(inst[k], width = 8), " range: [", .fmt(min(G[, k]), 4), ", ",
          .fmt(max(G[, k]), 4), "]   exact Pr scale: [",
          .fmt(suppressWarnings(min(pr, na.rm = TRUE)), 4), ", ",
          .fmt(suppressWarnings(max(pr, na.rm = TRUE)), 4), "]\n", sep = "")
    }
    cat("\n")
    bl <- x$bounds
    for (i in seq_len(nrow(bl))) {
      li <- if (x$target == "aie") lab else paste0("APE(", bl$target[i], ")")
      cat("    ", format(li, width = 18), " bounds = [", .fmt(bl$lower[i]), ", ",
          .fmt(bl$upper[i]), "]\n", sep = "")
    }
  } else {
    cat("\nLTZ prior:\n")
    for (k in seq_along(inst)) {
      if (!is.null(x$prior$mu)) {
        cat("    gamma_", inst[k], " ~ N(", .fmt(x$prior$mu[k], 4), ", ",
            .fmt(x$prior$omega[k], 4), ")", sep = "")
        if (x$gunits != "latent") cat("  [", x$gunits, " units]", sep = "")
        cat("\n")
      } else {
        cat("    gamma_", inst[k], " ~ ", x$prior$type, "\n", sep = "")
      }
      v <- x$sims$gamma_pr[, k]
      cat("      exact Pr scale of the draws: mean ", .fmt(mean(v, na.rm = TRUE), 4),
          ", sd ", .fmt(stats::sd(v, na.rm = TRUE), 4), "\n", sep = "")
    }
    cat("    iterations = ", x$iterations, if (x$n_failed > 0)
      paste0("  (", x$n_failed, " failed, excluded)"), "\n", sep = "")
    cat("\nLTZ ", lvl, "% interval (sampling noise convolved with the prior; CHR 2012 eq. 7):\n",
        sep = "")
    bl <- x$bounds
    for (i in seq_len(nrow(bl))) {
      li <- if (x$target == "aie") lab else paste0("APE(", bl$target[i], ")")
      cat("    ", format(li, width = 18), " ", lvl, "% CI = [", .fmt(bl$lower[i]), ", ",
          .fmt(bl$upper[i]), "]\n", sep = "")
    }
    cat("\nPrior-only spread of the point estimates (sensitivity band, no sampling noise):\n")
    pb <- x$prior_band
    for (i in seq_len(nrow(pb))) {
      li <- if (x$target == "aie") lab else paste0("APE(", pb$target[i], ")")
      cat("    ", format(li, width = 18), " ", lvl, "% range = [", .fmt(pb$lower[i]), ", ",
          .fmt(pb$upper[i]), "]\n", sep = "")
    }
    if (!is.null(x$sweep)) {
      cat("\nDelta sweep (first target):\n")
      print(format(x$sweep, digits = 6), row.names = FALSE)
    }
  }
  cat(hr, "\n", sep = "")
  invisible(x)
}

#' Summarize an apexog fit
#'
#' Prints the fit and, for \code{method = "uci"}, the estimates over the
#' \eqn{\gamma} grid; for \code{method = "ltz"}, quantiles of the refit
#' estimates and of the convolved draws.
#' @param object An \code{"apexog"} object.
#' @param ... Unused.
#' @return \code{object}, invisibly.
#' @export
summary.apexog <- function(object, ...) {
  print(object)
  if (object$method == "uci") {
    cat("\nEstimates over the gamma grid (latent units):\n")
    print(format(object$grid, digits = 6), row.names = FALSE)
  } else {
    cat("\nRefit estimates across the draws:\n")
    qs <- c(0, 0.025, 0.25, 0.5, 0.75, 0.975, 1)
    tab <- rbind(
      t(apply(object$sims$est, 2, stats::quantile, probs = qs, na.rm = TRUE)),
      t(apply(object$sims$est_cnv, 2, stats::quantile, probs = qs, na.rm = TRUE)))
    rownames(tab) <- c(paste(colnames(object$sims$est), "(point)"),
                       paste(colnames(object$sims$est_cnv), "(convolved)"))
    print(round(tab, 6))
  }
  invisible(object)
}

#' @rdname apexog-methods
#' @export
coef.apexog <- function(object, ...) object$coefficients

#' @rdname apexog-methods
#' @export
vcov.apexog <- function(object, ...) object$vcov

#' @rdname apexog-methods
#' @export
nobs.apexog <- function(object, ...) object$N

#' Plausibly-exogenous bounds
#'
#' Returns the UCI union (or the LTZ interval) for each target.
#' @param object An \code{"apexog"} object.
#' @param parm Optional subset of targets.
#' @param level Ignored except for a check: the level is fixed when fitting.
#' @param ... Unused.
#' @return A matrix with columns \code{lower} and \code{upper}.
#' @export
confint.apexog <- function(object, parm, level = object$level, ...) {
  if (!isTRUE(all.equal(level, object$level)))
    warning("the bounds were computed at level ", object$level,
            "; refit apexog() with level = ", level, ".", call. = FALSE)
  b <- object$bounds
  m <- cbind(lower = b$lower, upper = b$upper)
  rownames(m) <- b$target
  if (!missing(parm)) m <- m[parm, , drop = FALSE]
  m
}

#' Plot an apexog fit
#'
#' UCI: the point estimate and confidence band against the assumed
#' \eqn{\gamma} (single instrument; on the exact probability scale when
#' \code{gunits = "pr"}).  LTZ: a histogram of the convolved draws with the
#' interval and the baseline marked, or, after a delta sweep, the bounds
#' against the prior width.
#' @param x An \code{"apexog"} object.
#' @param which \code{"auto"} (default), \code{"grid"}, \code{"hist"}, or
#'   \code{"sweep"}.
#' @param ... Further arguments passed to the base graphics call.
#' @return \code{x}, invisibly.
#' @export
plot.apexog <- function(x, which = c("auto", "grid", "hist", "sweep"), ...) {
  which <- match.arg(which)
  if (which == "auto")
    which <- if (x$method == "uci") "grid" else if (!is.null(x$sweep)) "sweep" else "hist"
  lab <- if (x$target == "aie") "AIE" else "APE"
  if (which == "grid") {
    if (x$method != "uci") stop("which = \"grid\" needs method = \"uci\".", call. = FALSE)
    if (length(x$instruments) > 1L) {
      message("plot() of the UCI grid needs a single instrument; nothing drawn.")
      return(invisible(x))
    }
    t1 <- x$baseline$target[if (x$target == "aie") 1L else
      match(x$targets[1], x$baseline$target)]
    gx <- if (x$gunits == "pr") x$grid_pr[, 1] else x$grid[[1]]
    xl <- paste0("gamma_", x$instruments, if (x$gunits == "pr") " (exact Pr points)" else "")
    o <- order(gx)
    est <- x$grid[[paste0("est_", t1)]][o]
    lo <- x$grid[[paste0("lo_", t1)]][o]
    hi <- x$grid[[paste0("hi_", t1)]][o]
    gx <- gx[o]
    graphics::plot(gx, est, type = "l", col = "blue", ylim = range(c(lo, hi, 0), na.rm = TRUE),
                   xlab = xl, ylab = paste(lab, "(probability scale)"),
                   main = paste0("UCI ", lab, " bounds vs assumed gamma"), ...)
    graphics::lines(gx, lo, lty = 2)
    graphics::lines(gx, hi, lty = 2)
    graphics::abline(h = 0, lty = 3, col = "gray")
    graphics::legend("topright", c("Point estimate", paste0(100 * x$level, "% CI")),
                     lty = c(1, 2), col = c("blue", "black"), bty = "n")
  } else if (which == "hist") {
    if (x$method != "ltz") stop("which = \"hist\" needs method = \"ltz\".", call. = FALSE)
    v <- x$sims$est_cnv[, 1]
    graphics::hist(v, breaks = "FD", main = paste("LTZ", lab, "distribution across draws"),
                   xlab = paste(lab, "(probability scale)"), ...)
    graphics::abline(v = c(x$bounds$lower[1], x$bounds$upper[1]), lty = 2, col = "red")
    graphics::abline(v = x$baseline$estimate[1], col = "blue")
  } else {
    if (is.null(x$sweep)) stop("no delta sweep stored; use sweep_omega.", call. = FALSE)
    s <- x$sweep[order(x$sweep$delta), ]
    graphics::plot(s$delta, s$upper, type = "l", col = "blue",
                   ylim = range(c(s$lower, s$upper), na.rm = TRUE),
                   xlab = "delta = prior variance", ylab = paste(lab, "(probability scale)"),
                   main = "LTZ bounds vs prior width", ...)
    graphics::lines(s$delta, s$lower, col = "blue")
    graphics::lines(s$delta, s$lower_prior, lty = 2, col = "gray40")
    graphics::lines(s$delta, s$upper_prior, lty = 2, col = "gray40")
    graphics::abline(h = x$baseline$estimate[1], lty = 3)
    graphics::legend("topleft", c("LTZ CI", "prior-only band"), lty = c(1, 2),
                     col = c("blue", "gray40"), bty = "n")
  }
  invisible(x)
}
