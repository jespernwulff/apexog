# Internal helpers --------------------------------------------------------------

# Stata's round(): half away from zero (R's round() is half to even).
.stata_round <- function(x) sign(x) * floor(abs(x) + 0.5)

# Percentile positions used by apexog's LTZ intervals (Mata _apexogsim_pct):
#   lo_pos = max(1, round(a/2 * n)),  hi_pos = min(n, round((1 - a/2) * n))
# on the sorted non-missing draws.
.pct_interval <- function(x, level) {
  x <- sort(x[!is.na(x)])
  n <- length(x)
  if (n == 0L) return(c(NA_real_, NA_real_))
  a <- 1 - level
  lo <- max(1, .stata_round(a / 2 * n))
  hi <- min(n, .stata_round((1 - a / 2) * n))
  c(x[lo], x[hi])
}

# Generalized inverse in the spirit of Mata's invsym(): invert the block of
# the columns flagged in `keep`, zeros elsewhere.  With `keep` all TRUE this
# is solve(), falling back to a pseudo-inverse if the matrix is singular.
.ginv_sym <- function(H, keep = rep(TRUE, ncol(H))) {
  out <- matrix(0, nrow(H), ncol(H))
  if (!any(keep)) return(out)
  Hk <- H[keep, keep, drop = FALSE]
  inv <- tryCatch(solve(Hk), error = function(e) NULL)
  if (is.null(inv)) {
    s <- svd(Hk)
    tol <- max(dim(Hk)) * max(s$d) * .Machine$double.eps
    dinv <- ifelse(s$d > tol, 1 / s$d, 0)
    inv <- s$v %*% (dinv * t(s$u))
  }
  out[keep, keep] <- inv
  out
}

# Columns of `cand` that are linearly independent of `base` (and of earlier
# columns of `cand`), mirroring Stata's _rmcoll on [base, cand]: later
# collinear columns are dropped.
.independent_cols <- function(base, cand, tol = 1e-7) {
  if (ncol(cand) == 0L) return(integer(0))
  M <- cbind(base, cand)
  q <- qr(M, tol = tol)
  kept <- sort(q$pivot[seq_len(q$rank)])
  kept <- kept[kept > ncol(base)] - ncol(base)
  kept
}

.is01 <- function(x) all(x == 0 | x == 1)

.msg <- function(verbose, ...) if (isTRUE(verbose)) message(...)

# Prior draws for method = "ltz".  Returns an n x K matrix in the user's gamma
# units.  `prior` is "normal" (Gaussian via mu/omega), one of the named
# distribution types (with prior_params), or a function(n) returning n x K.
.draw_prior <- function(n, K, prior, mu, omega, prior_params) {
  if (is.function(prior)) {
    G <- as.matrix(prior(n))
    if (nrow(G) != n || ncol(G) != K)
      stop("prior() must return an ", n, " x ", K, " matrix (draws x instruments).",
           call. = FALSE)
    return(G)
  }
  G <- matrix(NA_real_, n, K)
  if (prior == "normal" && is.null(prior_params)) {
    for (i in seq_len(n)) G[i, ] <- mu + sqrt(omega) * stats::rnorm(K)
    return(G)
  }
  pp <- prior_params
  for (i in seq_len(n)) {
    for (k in seq_len(K)) {
      p <- pp[[k]]
      G[i, k] <- switch(prior,
        normal  = stats::rnorm(1, p[1], p[2]),
        uniform = p[1] + (p[2] - p[1]) * stats::runif(1),
        chi2    = stats::rchisq(1, p[1]),
        t       = stats::rt(1, p[1]),
        gamma   = stats::rgamma(1, shape = p[1], scale = p[2]),
        poisson = stats::rpois(1, p[1]),
        special = p[ceiling(stats::runif(1) * length(p))]
      )
    }
  }
  G
}

# Validate / broadcast prior_params to one parameter vector per instrument.
.check_prior <- function(prior, prior_params, K) {
  types <- c("normal", "uniform", "chi2", "t", "gamma", "poisson", "special")
  if (!prior %in% types)
    stop("prior must be a function or one of: ", paste(types, collapse = ", "),
         call. = FALSE)
  if (prior == "normal" && is.null(prior_params)) return(NULL)
  if (is.null(prior_params))
    stop("prior = \"", prior, "\" needs prior_params.", call. = FALSE)
  if (!is.list(prior_params)) prior_params <- list(prior_params)
  if (length(prior_params) == 1L && K > 1L) prior_params <- rep(prior_params, K)
  if (length(prior_params) != K)
    stop("prior_params must have one element per instrument (", K,
         ") or a single element (broadcast).", call. = FALSE)
  npar <- c(normal = 2, uniform = 2, chi2 = 1, t = 1, gamma = 2, poisson = 1)
  for (k in seq_len(K)) {
    p <- prior_params[[k]]
    if (!is.numeric(p)) stop("prior_params must be numeric.", call. = FALSE)
    if (prior == "special") {
      p <- p[!is.na(p)]
      if (length(p) < 1L) stop("prior = \"special\" needs a non-empty donor vector.",
                               call. = FALSE)
      prior_params[[k]] <- p
      next
    }
    if (length(p) != npar[[prior]])
      stop("prior = \"", prior, "\" takes ", npar[[prior]], " parameter(s) per instrument.",
           call. = FALSE)
    if (prior == "normal" && p[2] <= 0) stop("normal sd must be > 0.", call. = FALSE)
    if (prior == "uniform" && p[2] < p[1]) stop("uniform needs min <= max.", call. = FALSE)
    if (prior %in% c("chi2", "t") && p[1] < 1) stop(prior, " df must be >= 1.", call. = FALSE)
    if (prior == "poisson" && p[1] < 1) stop("poisson mean must be >= 1.", call. = FALSE)
    if (prior == "gamma" && any(p <= 0)) stop("gamma shape and scale must be > 0.", call. = FALSE)
  }
  prior_params
}

.broadcast <- function(x, K, what) {
  if (is.null(x)) return(NULL)
  if (length(x) == 1L) return(rep(x, K))
  if (length(x) != K)
    stop(what, " needs ", K, " value(s) (one per instrument) or a single value.",
         call. = FALSE)
  x
}
