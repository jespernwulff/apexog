#' Plausibly-exogenous bounds for average partial and interaction effects
#'
#' Sensitivity analysis for a binary-outcome instrumental-variables model in the
#' spirit of Conley, Hansen and Rossi (2012).  The exclusion restriction is
#' relaxed to \eqn{y^* = x'\beta + d\delta + z'\gamma + \varepsilon}, and the
#' target is the average partial effect (APE) of each endogenous regressor (or
#' an average interaction effect, AIE) on the probability scale, estimated by a
#' control-function probit (two-step standard errors: the Murphy-Topel
#' parameter covariance plus the direct first-stage term) or, for one
#' 0/1 treatment, by a recursive bivariate probit.  At every assumed violation
#' \eqn{\gamma} the model is refit with \eqn{z'\gamma} as an offset in the
#' outcome index.
#'
#' \describe{
#'   \item{\code{method = "uci"}}{Union of confidence intervals over a grid of
#'     assumed \eqn{\gamma} values (\code{gmin}, \code{gmax}, \code{grid}).}
#'   \item{\code{method = "ltz"}}{Local-to-zero: \eqn{\gamma} is drawn from a
#'     prior, the model is refit per draw, and the reported interval is the
#'     percentile interval of the refit estimates convolved with their sampling
#'     noise (CHR 2012, eq. 7).  The prior-only spread is reported as well.}
#' }
#'
#' This is an R port of the Stata command \code{apexog} (v0.21.2); results agree
#' with Stata to optimizer tolerance (about 1e-7 relative on the validation
#' specifications).  Random draws differ between R and Stata, so LTZ intervals
#' agree only in distribution unless the same \code{draws} are supplied.
#' Stata apexog up to v0.20.0 omitted the control-function-product
#' coefficients from the APE standard error when \code{cfinteract} was used
#' with an APE target (0.02\% of the SE in the validation design); v0.20.1 and
#' this package include them.  Since v0.21.0 (Stata and R) the standard error
#' of every two-step effect also includes the effect's direct dependence on the
#' estimated first-stage parameters through the control function (influence
#' function \code{J_beta' psi_beta + J_pi' psi_pi}); earlier versions, and
#' Stata's \code{cfprobit} followed by \code{margins}, use the first term only.
#' Point estimates are unchanged; the difference is about 1\% of the SE in the
#' validation designs.
#'
#' @param formula Two-sided formula \code{y ~ x1 + x2} with the 0/1 outcome and
#'   the exogenous controls (factors and interactions allowed; use
#'   \code{y ~ 1} for no controls).  A constant is always included.
#' @param endog Two-sided formula \code{d ~ z1 + z2} (endogenous regressor(s) on
#'   the left, excluded instruments on the right), or a list of such formulas,
#'   one per group of endogenous regressors sharing an instrument set (Stata's
#'   multi-clause syntax).  An endogenous interaction of a treatment with a
#'   moderator is written on the left as \code{d + d:x}, with its instruments
#'   on the right (e.g. \code{d + d:x ~ z + z:x}); combine it with
#'   \code{aie = c("d", "x")}.
#' @param data A data frame.
#' @param method \code{"uci"} or \code{"ltz"}.
#' @param gmin,gmax UCI: lower and upper end of the assumed \eqn{\gamma} range,
#'   one value per instrument or a single value (broadcast).
#' @param grid UCI: number of grid points per instrument (default 11).
#' @param mu,omega LTZ, Gaussian prior: mean and variance of each instrument's
#'   \eqn{\gamma} (single values are broadcast).
#' @param prior LTZ prior.  \code{"normal"} (default; uses \code{mu} and
#'   \code{omega} unless \code{prior_params} is given), one of \code{"uniform"},
#'   \code{"chi2"}, \code{"t"}, \code{"gamma"}, \code{"poisson"},
#'   \code{"special"} (parameters in \code{prior_params}), or a function
#'   \code{function(n)} returning an \code{n x K} matrix of draws.
#' @param prior_params List with one parameter vector per instrument (or one,
#'   broadcast): normal \code{c(mean, sd)}, uniform \code{c(min, max)}, chi2
#'   \code{df}, t \code{df}, gamma \code{c(shape, scale)}, poisson \code{mean},
#'   special: a numeric donor vector sampled with replacement.
#' @param draws Optional \code{iterations x K} matrix of \eqn{\gamma} draws
#'   (in \code{gunits}) used instead of a prior, e.g. a joint prior of your own.
#'   Passing Stata's \code{e(gamma_sim)} (latent units, so with the default
#'   \code{gunits = "latent"}) reproduces a Stata LTZ run draw by draw; only
#'   the convolution noise then differs.
#' @param iterations LTZ: number of prior draws (default 1000).
#' @param seed Optional seed passed to \code{set.seed()} before drawing.
#' @param firststage First-stage model: \code{"linear"} (default),
#'   \code{"probit"} (0/1 regressor, generalized residual), \code{"fprobit"}
#'   (fractional regressor), \code{"poisson"} (count), or \code{"biprobit"}
#'   (recursive bivariate probit for a single 0/1 regressor; exact under joint
#'   normality).  A single value applies to every endogenous regressor,
#'   including a declared interaction such as \code{d:x} (which needs
#'   \code{"linear"}); otherwise give one value per endogenous regressor, in
#'   order of appearance, or a named vector (unnamed regressors default to
#'   \code{"linear"}), e.g. \code{c(d = "probit")}.
#' @param target Character vector of endogenous regressors whose APE bounds are
#'   reported (default: all).
#' @param aie Character vector \code{c(X, Z)}: switch the target to the average
#'   interaction effect of \code{X} and \code{Z}.  Each must be an exogenous
#'   control column or an endogenous regressor (use numeric 0/1 variables, not
#'   \code{factor()}).  If the product \code{X:Z} is a declared endogenous
#'   regressor (\code{endog = d + d:x ~ ...}) that regressor is the target;
#'   otherwise the product is taken from \code{formula} when it is written as
#'   \code{x:z} or \code{x*z}, or added as an exogenous control.  Do not write
#'   it as \code{I(x*z)}: that is detected as a duplicate and refused.
#' @param aie_binary Logical length-2 vector forcing the type of \code{X} and
#'   \code{Z} (\code{TRUE} = binary, \code{FALSE} = continuous, \code{NA} =
#'   detect from the data).  Continuous x continuous gives the cross-partial,
#'   continuous x binary the first difference of the slope, binary x binary the
#'   second discrete difference.
#' @param cfinteract Add the products of the control functions to the main
#'   probit (as \code{cfprobit}'s \code{cfinteract}).  Automatic for a declared
#'   endogenous interaction.
#' @param gunits Units of all user-supplied \eqn{\gamma} values:
#'   \code{"latent"} (probit-index units; default), \code{"pr"} (probability
#'   points, divided by \eqn{m_\phi}, the average density at the baseline fit),
#'   or \code{"structural"} (structural-latent units, divided by
#'   \eqn{\tilde\sigma}; derived for a linear first stage).  The structural
#'   conversion uses the baseline \eqn{\tilde\sigma}, which is exact only at
#'   \eqn{\gamma = 0}.  With one continuous endogenous regressor, one
#'   instrument and a linear first stage, a structural value has a unique
#'   index-scale counterpart only below \eqn{|\pi_z|/\sigma_v} in absolute
#'   value (returned as \code{gs_max}); \code{print} reports that limit and
#'   warns when the UCI box or LTZ prior draws reach it.  For a 0/1 regressor
#'   with a linear first stage the conversion is not derived and \code{print}
#'   warns of that instead.  Results are not changed.
#' @param level Confidence level, as a fraction or a percent (default 0.95).
#' @param cluster Optional cluster identifier: a one-sided formula
#'   (\code{~ id}), a column name, or a vector of length \code{nrow(data)}.
#' @param subset Optional logical expression evaluated in \code{data}.
#'   Observations where it is \code{NA} are dropped (Stata's \code{if} treats
#'   a missing value as \code{+Inf} and may keep them).
#' @param sweep_omega,sweep_mu,sweep_delta LTZ, Gaussian prior: rerun the
#'   simulation for each prior variance in \code{sweep_omega} (means
#'   \code{sweep_mu}, default 0; x-axis labels \code{sweep_delta}, default the
#'   variances) and store the bounds-versus-width table in \code{$sweep}.
#' @param verbose Print notes while fitting (default \code{TRUE}).
#'
#' @return An object of class \code{"apexog"}: a list with, among others,
#'   \code{baseline} (estimates at \eqn{\gamma = 0}), \code{bounds} (UCI union or
#'   LTZ interval per target), \code{prior_band} (LTZ), \code{grid} and
#'   \code{grid_pr} (UCI), \code{sims} (LTZ draws and refit estimates),
#'   \code{mphi}, \code{stilde}, \code{gs_max} (the limit of the structural
#'   scale, \code{NA} when not defined), \code{binlin} (0/1 endogenous
#'   regressors that got a linear first stage, for which the effect is a
#'   descriptive fitted contrast rather than a causal estimate; \code{print}
#'   says so), \code{coefficients} and \code{vcov} (baseline
#'   targets), and \code{model} (the baseline main-equation coefficients with
#'   their Murphy-Topel covariance, or the bivariate-probit coefficients with
#'   their robust MLE covariance).  Methods: \code{print}, \code{summary}, \code{coef}, \code{vcov},
#'   \code{confint}, \code{plot}.
#'
#' @references Conley, T. G., Hansen, C. B. and Rossi, P. E. (2012).
#'   Plausibly exogenous.  \emph{Review of Economics and Statistics} 94(1),
#'   260-272.
#'
#'   Murphy, K. M. and Topel, R. H. (1985).  Estimation and inference in
#'   two-step econometric models.  \emph{Journal of Business and Economic
#'   Statistics} 3(4), 370-379.
#'
#'   Wooldridge, J. M. (2014).  Quasi-maximum likelihood estimation and testing
#'   for nonlinear models with endogenous explanatory variables.
#'   \emph{Journal of Econometrics} 182(1), 226-234.
#'
#' @examples
#' set.seed(1)
#' n <- 2000
#' dat <- data.frame(x = rnorm(n), z = rnorm(n), v = rnorm(n))
#' dat$d <- 0.7 * dat$z + 0.3 * dat$x + dat$v
#' dat$y <- as.integer(-0.3 + 0.5 * dat$d + 0.4 * dat$x +
#'                     0.5 * dat$v + sqrt(0.75) * rnorm(n) > 0)
#' fit <- apexog(y ~ x, endog = d ~ z, data = dat, method = "uci",
#'               gmin = -0.1, gmax = 0.1, grid = 5)
#' fit
#' confint(fit)
#' fit2 <- apexog(y ~ x, endog = d ~ z, data = dat, method = "ltz",
#'                mu = 0, omega = 0.0025, iterations = 100, seed = 1)
#' summary(fit2)
#' @export
apexog <- function(formula, endog, data, method = c("uci", "ltz"),
                   gmin = NULL, gmax = NULL, grid = 11,
                   mu = NULL, omega = NULL, prior = "normal", prior_params = NULL,
                   draws = NULL, iterations = 1000, seed = NULL,
                   firststage = "linear", target = NULL,
                   aie = NULL, aie_binary = NULL, cfinteract = FALSE,
                   gunits = c("latent", "pr", "structural"),
                   level = 0.95, cluster = NULL, subset = NULL,
                   sweep_omega = NULL, sweep_mu = NULL, sweep_delta = NULL,
                   verbose = TRUE) {
  cl <- match.call()
  method <- match.arg(method)
  if (length(gunits) > 1L) gunits <- gunits[1L]
  gunits <- switch(gunits, index = "latent", probability = "pr", gunits)
  if (!gunits %in% c("latent", "pr", "structural"))
    stop("gunits must be \"latent\" (synonym \"index\"), \"pr\" or \"structural\".",
         call. = FALSE)
  if (level >= 1 && level < 100) level <- level / 100
  if (level <= 0 || level >= 1)
    stop("level must be in (0, 1) (or 1-99.99 as a percent).", call. = FALSE)
  if (!is.data.frame(data)) stop("data must be a data frame.", call. = FALSE)
  if (!inherits(formula, "formula") || length(formula) != 3L)
    stop("formula must be two-sided: y ~ controls (y ~ 1 for none).", call. = FALSE)
  if (inherits(endog, "formula")) endog <- list(endog)
  if (!is.list(endog) || length(endog) < 1L ||
      !all(vapply(endog, function(f) inherits(f, "formula") && length(f) == 3L,
                  logical(1))))
    stop("endog must be a two-sided formula (d ~ z) or a list of them.", call. = FALSE)
  if (!is.null(sweep_omega) && method != "ltz")
    stop("sweep_omega/sweep_mu/sweep_delta (the delta sweep) require method = \"ltz\".",
         call. = FALSE)
  if (method == "ltz" && (!is.null(gmin) || !is.null(gmax)))
    .msg(verbose, "  [note] gmin/gmax are UCI options; ignored with method = \"ltz\".")
  if (method == "uci" && (!is.null(mu) || !is.null(omega) || !is.null(draws) ||
                          !is.null(prior_params) || !is.null(seed) ||
                          !identical(prior, "normal")))
    .msg(verbose, "  [note] mu/omega/prior/draws/seed are LTZ options; ignored with method = \"uci\".")

  # ---- estimation sample ----------------------------------------------------
  clvar <- NULL
  if (!is.null(cluster)) {
    if (inherits(cluster, "formula")) clvar <- all.vars(cluster)
    else if (is.character(cluster) && length(cluster) == 1L) clvar <- cluster
    else {
      if (length(cluster) != nrow(data))
        stop("cluster must be a formula, a column name, or a vector of length nrow(data).",
             call. = FALSE)
      data$.apexog_cluster <- cluster
      clvar <- ".apexog_cluster"
    }
    if (length(clvar) != 1L) stop("cluster must name one variable.", call. = FALSE)
  }
  vars <- unique(c(all.vars(formula), unlist(lapply(endog, all.vars)), clvar))
  miss <- setdiff(vars, names(data))
  if (length(miss)) stop("variable(s) not found in data: ", paste(miss, collapse = ", "),
                         call. = FALSE)
  sub <- substitute(subset)
  dat <- data
  if (!is.null(sub)) {
    keep <- eval(sub, data, parent.frame())
    if (!is.logical(keep) || length(keep) != nrow(data))
      stop("subset must be a logical expression evaluated in data.", call. = FALSE)
    dat <- dat[!is.na(keep) & keep, , drop = FALSE]
  }
  dat <- dat[stats::complete.cases(dat[, vars, drop = FALSE]), , drop = FALSE]
  dat <- droplevels(dat)
  N <- nrow(dat)
  if (N == 0L) stop("No observations remain after subset and missing-data filtering.",
                    call. = FALSE)

  # ---- outcome and exogenous controls --------------------------------------
  mf <- stats::model.frame(formula, dat, na.action = stats::na.pass)
  y <- stats::model.response(mf)
  if (is.logical(y)) y <- as.numeric(y)
  if (is.factor(y) || !is.numeric(y)) stop("The outcome must be numeric 0/1.", call. = FALSE)
  if (!.is01(y)) stop("The outcome ", deparse(formula[[2]]), " takes values other than 0/1 in ",
                      sum(!(y == 0 | y == 1)), " observation(s).", call. = FALSE)
  y <- as.numeric(y)
  ylab <- deparse(formula[[2]])
  X <- stats::model.matrix(attr(mf, "terms"), mf)
  X <- X[, colnames(X) != "(Intercept)", drop = FALSE]
  storage.mode(X) <- "double"
  if (ncol(X) > 0L) {
    kx <- .independent_cols(matrix(1, N, 1), X)
    if (length(kx) < ncol(X)) {
      .msg(verbose, "  [note] control(s) dropped as constant or collinear: ",
           paste(colnames(X)[-kx], collapse = ", "))
      X <- X[, kx, drop = FALSE]
    }
  }
  X0 <- X

  # ---- endogenous regressors and instruments --------------------------------
  endog_names <- character(0)
  dcols <- list()
  clause_of <- integer(0)
  zmats <- list()
  aieendog <- NULL
  for (c in seq_along(endog)) {
    f <- endog[[c]]
    lhs <- attr(stats::terms(stats::as.formula(call("~", f[[2]]))), "term.labels")
    if (!length(lhs)) stop("Endogenous clause ", c, " has no regressor.", call. = FALSE)
    for (lab in lhs) {
      if (lab %in% endog_names) next
      if (grepl(":", lab, fixed = TRUE)) {
        ops <- strsplit(lab, ":", fixed = TRUE)[[1]]
        if (length(ops) != 2L)
          stop("Only two-way interaction endogenous regressors are supported; got ", lab,
               call. = FALSE)
        if (!is.null(aieendog))
          stop("Only one interaction endogenous regressor is supported.", call. = FALSE)
        v1 <- eval(str2lang(ops[1]), dat, environment(f))
        v2 <- eval(str2lang(ops[2]), dat, environment(f))
        if (!is.numeric(v1) || !is.numeric(v2))
          stop("Interaction operands must be numeric: ", lab, call. = FALSE)
        val <- as.numeric(v1) * as.numeric(v2)
        aieendog <- list(name = lab, ops = ops)
      } else {
        val <- eval(str2lang(lab), dat, environment(f))
        if (is.logical(val)) val <- as.numeric(val)
        if (!is.numeric(val))
          stop("Endogenous regressor ", lab, " must be numeric.", call. = FALSE)
        val <- as.numeric(val)
      }
      endog_names <- c(endog_names, lab)
      dcols[[lab]] <- val
      clause_of <- c(clause_of, c)
    }
    rhs <- f[-2L]
    Zc <- stats::model.matrix(rhs, stats::model.frame(rhs, dat, na.action = stats::na.pass))
    Zc <- Zc[, colnames(Zc) != "(Intercept)", drop = FALSE]
    storage.mode(Zc) <- "double"
    if (ncol(Zc) == 0L) stop("Endogenous clause ", c, " has no instrument.", call. = FALSE)
    zmats[[c]] <- Zc
  }
  J <- length(endog_names)
  D <- do.call(cbind, dcols)
  colnames(D) <- endog_names
  if (ylab %in% c(colnames(X), endog_names))
    stop("The outcome also appears among the regressors.", call. = FALSE)
  ov <- intersect(endog_names, colnames(X))
  if (length(ov)) stop("Variable(s) listed as both endogenous and exogenous: ",
                       paste(ov, collapse = ", "), call. = FALSE)

  # ---- instruments: full rank against the controls, pooled -------------------
  Zpool <- NULL
  raw_n <- length(unique(unlist(lapply(zmats, colnames))))
  iv_of <- vector("list", J)
  for (j in seq_len(J)) {
    Zc <- zmats[[clause_of[j]]]
    ov <- intersect(endog_names, colnames(Zc))
    if (length(ov)) stop("Variable(s) listed as both endogenous regressor and instrument: ",
                         paste(ov, collapse = ", "), call. = FALSE)
    if (ylab %in% colnames(Zc)) stop("The outcome also appears among the instruments.",
                                     call. = FALSE)
    keep <- .independent_cols(cbind(1, X0), Zc)
    nm <- colnames(Zc)[keep]
    for (k in nm) if (is.null(Zpool) || !k %in% colnames(Zpool))
      Zpool <- cbind(Zpool, Zc[, k, drop = FALSE])
    iv_of[[j]] <- nm
  }
  if (is.null(Zpool)) stop("No usable instruments remain after dropping collinear terms.",
                           call. = FALSE)
  Z <- Zpool
  K <- ncol(Z)
  if (raw_n > K) {
    .msg(verbose, "  [note] ", raw_n - K, " of ", raw_n,
         " instrument term(s) dropped (constant or collinear with the controls / other instruments).")
    .msg(verbose, "         ", K, " full-rank instrument(s) retained: ",
         paste(colnames(Z), collapse = " "))
  }
  if (K < J) stop("Need at least as many instruments (", K, ") as endogenous regressors (",
                  J, ").", call. = FALSE)

  # ---- first-stage types and engine -------------------------------------------
  types <- c("linear", "probit", "fprobit", "poisson", "biprobit")
  fs <- firststage
  if (!is.null(names(fs)) && all(names(fs) != "")) {
    if (!all(names(fs) %in% endog_names))
      stop("Names of firststage must be endogenous regressors.", call. = FALSE)
    fs2 <- rep("linear", J); names(fs2) <- endog_names
    fs2[names(fs)] <- fs
    fs <- unname(fs2)
  } else if (length(fs) == 1L) {
    fs <- rep(fs, J)
  } else if (length(fs) != J) {
    stop("firststage must have 1 value or one per endogenous regressor (", J, ").",
         call. = FALSE)
  }
  if (!all(fs %in% types))
    stop("Unknown firststage type; valid: ", paste(types, collapse = ", "), call. = FALSE)
  engine <- if (any(fs == "biprobit")) "biprobit" else "cf"
  if (engine == "biprobit" && J > 1L)
    stop("firststage = \"biprobit\" supports exactly one endogenous regressor.", call. = FALSE)
  for (j in seq_len(J)) {
    dj <- D[, j]
    if (fs[j] %in% c("probit", "biprobit") && !.is01(dj))
      stop("firststage = \"", fs[j], "\" requires a 0/1 endogenous regressor; ",
           endog_names[j], " has other values.", call. = FALSE)
    if (fs[j] == "fprobit" && any(dj < 0 | dj > 1))
      stop("firststage = \"fprobit\" requires a regressor in [0, 1]; ", endog_names[j],
           " has values outside that range.", call. = FALSE)
    if (fs[j] == "poisson" && any(dj < 0))
      stop("firststage = \"poisson\" requires a nonnegative regressor; ", endog_names[j],
           " has negative values.", call. = FALSE)
  }

  # ---- targets ---------------------------------------------------------------
  targets <- endog_names
  if (!is.null(target)) {
    if (!is.null(aie)) stop("target and aie are mutually exclusive.", call. = FALSE)
    bad <- setdiff(target, endog_names)
    if (length(bad)) stop("target ", paste(bad, collapse = ", "),
                          " is not an endogenous regressor.", call. = FALSE)
    targets <- endog_names[endog_names %in% target]
  }
  aie_info <- NULL
  aie_xz_is_endog <- FALSE
  if (!is.null(aie)) {
    if (length(aie) != 2L) stop("aie must name exactly two variables.", call. = FALSE)
    if (aie[1] == aie[2]) stop("The two AIE interactors must be distinct.", call. = FALSE)
    if (engine == "biprobit")
      stop("The AIE target is not available with firststage = \"biprobit\".", call. = FALSE)
    getv <- function(v) {
      if (v %in% colnames(X)) return(X[, v])
      if (v %in% endog_names && !grepl(":", v, fixed = TRUE)) return(D[, v])
      stop("AIE interactor ", v, " must be an exogenous control column or an endogenous regressor.",
           call. = FALSE)
    }
    xv <- getv(aie[1]); zv <- getv(aie[2])
    if (is.null(aie_binary)) aie_binary <- c(NA, NA)
    if (length(aie_binary) != 2L) stop("aie_binary must have length 2.", call. = FALSE)
    xb <- if (is.na(aie_binary[1])) .is01(xv) else isTRUE(aie_binary[1])
    zb <- if (is.na(aie_binary[2])) .is01(zv) else isTRUE(aie_binary[2])
    case <- 2L * as.integer(xb) + as.integer(zb)
    lbl <- c("continuous-by-continuous (cross-partial)",
             paste0("continuous-by-binary (first difference in ", aie[2], ")"),
             paste0("binary-by-continuous (first difference in ", aie[1], ")"),
             "binary-by-binary (second discrete difference)")[case + 1L]
    if (!is.null(aieendog) && setequal(aieendog$ops, aie)) {
      aie_xz_is_endog <- TRUE
      xzname <- aieendog$name
    } else {
      cand <- c(paste0(aie[1], ":", aie[2]), paste0(aie[2], ":", aie[1]))
      hit <- intersect(cand, colnames(X))
      if (length(hit)) {
        xzname <- hit[1]
      } else {
        xzname <- cand[1]
        if (!length(.independent_cols(cbind(1, X), cbind(xv * zv))))
          stop("The product of ", aie[1], " and ", aie[2], " is already among the ",
               "controls under another name (e.g. I(x*z)); write it as ", xzname,
               " or leave it out.", call. = FALSE)
        X <- cbind(X, xv * zv)
        colnames(X)[ncol(X)] <- xzname
      }
    }
    aie_info <- list(vars = aie, case = case, label = lbl, xz = xzname,
                     binary = c(xb, zb))
    .msg(verbose, "  [aie] target = ", lbl)
    if (any(aie %in% endog_names))
      .msg(verbose, "  [note] AIE: one or both interactors are endogenous; the MT correction propagates first-stage uncertainty into the AIE SE.")
    targets <- "AIE"
  }
  p <- ncol(X)

  do_cf <- isTRUE(cfinteract) || aie_xz_is_endog
  if (J < 2L) {
    if (isTRUE(cfinteract)) .msg(verbose, "  [note] cfinteract needs at least two endogenous regressors; ignored.")
    do_cf <- FALSE
  }

  # ---- masks and first stages -------------------------------------------------
  IVMASK <- matrix(FALSE, J, K)
  for (j in seq_len(J)) IVMASK[j, ] <- colnames(Z) %in% iv_of[[j]]
  EXOGMASK <- matrix(TRUE, J, p)
  if (!is.null(aie_info) && !aie_xz_is_endog) {
    xzpos <- match(aie_info$xz, colnames(X))
    for (j in seq_len(J)) if (endog_names[j] %in% aie) EXOGMASK[j, xzpos] <- FALSE
  }
  W <- cbind(Z, X, 1)
  wcols <- lapply(seq_len(J), function(k) c(which(IVMASK[k, ]), K + which(EXOGMASK[k, ]),
                                           K + p + 1L))
  CHAT <- AA <- matrix(NA_real_, N, J)
  pi_first <- NULL
  if (engine == "cf") {
    for (j in seq_len(J)) {
      fsj <- .first_stage(D[, j], W[, wcols[[j]], drop = FALSE], fs[j], endog_names[j])
      CHAT[, j] <- fsj$chat
      AA[, j] <- fsj$A
      if (j == 1L) pi_first <- fsj$coef
    }
  }
  # v0.21.2 (display only; no computation uses either).  binlin: 0/1 endogenous
  # regressors that got a linear first stage (the product regressor of an
  # endogenous interaction is skipped; its treatment carries the flag).
  # gs_max: |pi_z| / sigma_v, the limit of the structural scale, defined for one
  # continuous (not 0/1) endogenous regressor with one instrument and a linear
  # first stage.
  binlin <- character(0)
  gs_max <- NA_real_
  if (engine == "cf") {
    binlin <- endog_names[fs == "linear" & apply(D, 2, .is01) &
                            !grepl(":", endog_names, fixed = TRUE)]
    if (J == 1L && K == 1L && fs[1] == "linear" && !is.null(pi_first) && !.is01(D[, 1])) {
      sdv <- stats::sd(CHAT[, 1])
      if (is.finite(sdv) && sdv > 0) gs_max <- abs(unname(pi_first[1])) / sdv
    }
  }
  colnames(CHAT) <- paste0("cf_", endog_names)
  cf_a <- cf_b <- integer(0)
  CFPROD <- NULL
  if (do_cf) {
    for (a in seq_len(J - 1L)) for (b in (a + 1L):J) {
      cf_a <- c(cf_a, a); cf_b <- c(cf_b, b)
      CFPROD <- cbind(CFPROD, CHAT[, a] * CHAT[, b])
    }
    colnames(CFPROD) <- paste0("cf_", endog_names[cf_a], "*cf_", endog_names[cf_b])
    .msg(verbose, "  [cfinteract] added ", ncol(CFPROD),
         " control-function product term(s) to the main probit.")
  }
  if (J > 1L)
    .msg(verbose, "  [note] ", J, " endogenous regressors: all get a control function; ",
         "the MT correction accumulates uncertainty from all first stages.")

  xmain <- cbind(X, D, if (engine == "cf") CHAT, CFPROD, `(Intercept)` = 1)
  colnames(xmain)[ncol(xmain)] <- "(Intercept)"
  d_idx <- p + seq_len(J)
  chat_idx <- p + J + seq_len(J)
  cfp_idx <- if (do_cf) p + 2L * J + seq_along(cf_a) else integer(0)

  clid <- NULL
  if (!is.null(clvar)) clid <- as.integer(factor(dat[[clvar]]))

  S <- list(y = y, X = X, D = D, CHAT = CHAT, A = AA, Z = Z, W = W, wcols = wcols,
            xmain = xmain, p = p, d_idx = d_idx, chat_idx = chat_idx,
            cfp_idx = cfp_idx, cf_a = cf_a, cf_b = cf_b,
            load_cols = c(d_idx, if (engine == "cf") chat_idx, cfp_idx),
            cluster = clid, level = level,
            dbin = apply(D, 2, .is01), zbin = apply(Z, 2, .is01))
  if (!is.null(aie_info)) {
    S$aie <- list(case = aie_info$case,
                  x_idx = match(aie[1], colnames(xmain)),
                  z_idx = match(aie[2], colnames(xmain)),
                  xz_idx = match(aie_info$xz, colnames(xmain)))
  }
  if (engine == "biprobit") {
    S$bp_x1 <- cbind(X, D, 1)
    S$bp_x2 <- cbind(Z, X, 1)
    S$bp_names <- c(paste0(ylab, ":", c(colnames(X), endog_names, "(Intercept)")),
                    paste0(endog_names, ":", c(colnames(Z), colnames(X), "(Intercept)")),
                    "athrho")
  }
  fitter <- if (engine == "cf") .fit_cf else .fit_bp
  ntg <- if (is.null(aie_info)) J else 1L
  tnames <- if (is.null(aie_info)) endog_names else "AIE"

  # ---- baseline (gamma = 0) ----------------------------------------------------
  gam0 <- rep(0, K)
  f0 <- fitter(S, gam0)
  if (!f0$ok) .msg(TRUE, "  [warn] baseline fit failed: ", f0$msg)
  warm <- f0$coef
  rho0 <- if (engine == "biprobit") f0$rho else NULL
  baseline <- data.frame(target = tnames, estimate = f0$est, se = f0$se,
                         lower = f0$lo, upper = f0$hi, stringsAsFactors = FALSE)
  vt0 <- f0$vtgt
  dimnames(vt0) <- list(tnames, tnames)

  ginv <- 1
  if (gunits == "pr") {
    if (!is.finite(f0$mphi)) stop("gunits = \"pr\" needs a successful baseline fit.", call. = FALSE)
    ginv <- 1 / f0$mphi
  } else if (gunits == "structural") {
    if (!is.finite(f0$stilde))
      stop("gunits = \"structural\" needs a successful baseline fit.", call. = FALSE)
    ginv <- 1 / f0$stilde
  }

  safe_fit <- function(g) {
    r <- tryCatch(fitter(S, g, warm), error = function(e) .na_fit(ntg, K, conditionMessage(e)))
    if (!r$ok) .msg(verbose, "  [warn] ", r$msg, " at gamma = (",
                    paste(signif(g, 4), collapse = ", "), "); returning NA for this fit.")
    r
  }

  res <- list(call = cl, method = method,
              estimator = if (engine == "cf") "cfprobit" else "biprobit",
              firststage = fs, vce = if (engine == "cf") "Murphy-Topel" else "MLE",
              N = N, N_clust = if (is.null(clid)) NULL else length(unique(clid)),
              clustvar = clvar, level = level, gunits = gunits,
              mphi = f0$mphi, stilde = f0$stilde, gs_max = gs_max, binlin = binlin,
              rho = rho0,
              depvar = ylab, endog = endog_names, exog = colnames(X0),
              instruments = colnames(Z), targets = targets,
              target = if (is.null(aie_info)) "ape" else "aie", aie = aie_info,
              baseline = baseline, coefficients = stats::setNames(f0$est, tnames),
              vcov = vt0, ginv = ginv,
              model = list(coefficients = f0$coef,
                           vcov = if (engine == "cf") f0$Vbeta else f0$V))
  if (!is.null(res$model$vcov) && !is.null(names(f0$coef)))
    dimnames(res$model$vcov) <- list(names(f0$coef), names(f0$coef))

  # ---- UCI ---------------------------------------------------------------------
  if (method == "uci") {
    if (is.null(gmin) || is.null(gmax)) stop("method = \"uci\" needs gmin and gmax.",
                                            call. = FALSE)
    gmin <- .broadcast(gmin, K, "gmin")
    gmax <- .broadcast(gmax, K, "gmax")
    if (any(gmin > gmax)) stop("gmin must not exceed gmax.", call. = FALSE)
    if (grid < 2) stop("grid must be at least 2.", call. = FALSE)
    total <- grid^K
    if (total > 11000)
      .msg(verbose, "  [note] the UCI grid has ", total, " points (grid^#instruments); ",
           "consider a coarser grid or method = \"ltz\".")
    .msg(verbose, "Running UCI grid: ", grid, " point(s) per instrument, ", K,
         " instrument(s), ", total, " total grid point(s).")
    lo_l <- gmin * ginv
    hi_l <- gmax * ginv
    idx <- as.matrix(expand.grid(rep(list(0:(grid - 1)), K)))
    G <- matrix(NA_real_, total, K)
    est <- lo <- hi <- matrix(NA_real_, total, ntg)
    gpr <- matrix(NA_real_, total, K)
    for (r in seq_len(total)) {
      g <- lo_l + ((hi_l - lo_l) / (grid - 1)) * idx[r, ]
      G[r, ] <- g
      f <- safe_fit(g)
      est[r, ] <- f$est; lo[r, ] <- f$lo; hi[r, ] <- f$hi
      if (f$ok) gpr[r, ] <- f$gpr
    }
    colnames(G) <- paste0("g_", colnames(Z))
    grid_df <- data.frame(G, check.names = FALSE)
    for (t in seq_len(ntg)) {
      grid_df[[paste0("est_", tnames[t])]] <- est[, t]
      grid_df[[paste0("lo_", tnames[t])]] <- lo[, t]
      grid_df[[paste0("hi_", tnames[t])]] <- hi[, t]
    }
    colnames(gpr) <- paste0("gpr_", colnames(Z))
    lb <- apply(lo, 2, function(x) if (all(is.na(x))) NA_real_ else min(x, na.rm = TRUE))
    ub <- apply(hi, 2, function(x) if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE))
    res$bounds <- data.frame(target = tnames, lower = lb, upper = ub,
                             stringsAsFactors = FALSE)[tnames %in% targets, , drop = FALSE]
    res$grid <- grid_df
    res$grid_pr <- gpr
    res$gmin <- gmin; res$gmax <- gmax; res$gridn <- grid
    res$n_failed <- sum(rowSums(is.na(est)) > 0)
  }

  # ---- LTZ ---------------------------------------------------------------------
  if (method == "ltz") {
    gauss <- FALSE
    if (!is.null(draws)) {
      Gu <- as.matrix(draws)
      if (ncol(Gu) != K) stop("draws must have one column per instrument (", K, ").",
                              call. = FALSE)
      iterations <- nrow(Gu)
      prior_lab <- "user-supplied draws"
    } else {
      if (is.function(prior)) {
        prior_lab <- "user-supplied prior function"
      } else {
        prior_params <- .check_prior(prior, prior_params, K)
        if (prior == "normal" && is.null(prior_params)) {
          if (is.null(mu) || is.null(omega))
            stop("method = \"ltz\" needs mu and omega (or prior/prior_params/draws).",
                 call. = FALSE)
          mu <- .broadcast(mu, K, "mu")
          omega <- .broadcast(omega, K, "omega")
          if (any(omega <= 0)) stop("All omega values must be > 0.", call. = FALSE)
          gauss <- TRUE
          prior_lab <- "gaussian"
        } else prior_lab <- prior
      }
      if (iterations < 10) stop("iterations must be >= 10.", call. = FALSE)
      if (!is.null(seed)) set.seed(seed)
      Gu <- .draw_prior(iterations, K, prior, mu, omega, prior_params)
    }
    .msg(verbose, "Running LTZ: ", iterations, " draw(s) from the ", prior_lab,
         " prior over ", K, " instrument violation(s).")
    G <- Gu * ginv
    est <- se <- matrix(NA_real_, iterations, ntg)
    gpr <- matrix(NA_real_, iterations, K)
    for (i in seq_len(iterations)) {
      f <- safe_fit(G[i, ])
      est[i, ] <- f$est; se[i, ] <- f$se
      if (f$ok) gpr[i, ] <- f$gpr
      if (verbose && i %% 100 == 0) message("  ... ", i, " / ", iterations, " draws done")
    }
    nfail <- sum(rowSums(is.na(cbind(est, se))) > 0)
    if (nfail > 0)
      .msg(TRUE, "  [note] ", nfail, " of ", iterations,
           " LTZ draws failed; excluded from the percentile interval.")
    eps <- matrix(stats::rnorm(iterations * ntg), iterations, ntg)
    cnv <- est + se * eps
    pct <- t(apply(cnv, 2, .pct_interval, level = level))
    prb <- t(apply(est, 2, .pct_interval, level = level))
    colnames(G) <- colnames(gpr) <- colnames(Z)
    colnames(est) <- colnames(se) <- colnames(cnv) <- tnames
    res$bounds <- data.frame(target = tnames, lower = pct[, 1], upper = pct[, 2],
                             stringsAsFactors = FALSE)[tnames %in% targets, , drop = FALSE]
    res$prior_band <- data.frame(target = tnames, lower = prb[, 1], upper = prb[, 2],
                                 stringsAsFactors = FALSE)[tnames %in% targets, , drop = FALSE]
    res$sims <- list(gamma = G, gamma_pr = gpr, est = est, se = se, est_cnv = cnv)
    res$prior <- list(type = prior_lab, mu = if (gauss) mu, omega = if (gauss) omega,
                      params = if (!gauss && !is.function(prior)) prior_params)
    res$iterations <- iterations
    res$n_failed <- nfail

    # delta sweep
    if (!is.null(sweep_omega)) {
      if (!gauss) stop("sweep_omega requires the Gaussian mu/omega prior.", call. = FALSE)
      ns <- length(sweep_omega)
      if (any(sweep_omega <= 0)) stop("sweep_omega values must be > 0.", call. = FALSE)
      smu <- if (is.null(sweep_mu)) rep(0, ns) else sweep_mu
      if (length(smu) != ns) stop("sweep_mu must have as many values as sweep_omega.",
                                  call. = FALSE)
      sde <- if (is.null(sweep_delta)) sweep_omega else sweep_delta
      if (length(sde) != ns) stop("sweep_delta must have as many values as sweep_omega.",
                                  call. = FALSE)
      sw <- matrix(NA_real_, ns, 6,
                   dimnames = list(NULL, c("delta", "mu", "lower", "upper",
                                           "lower_prior", "upper_prior")))
      .msg(verbose, "Delta sweep: ", ns, " prior width(s) x ", iterations,
           " draws each (first target).")
      for (s in seq_len(ns)) {
        a1 <- b1 <- numeric(iterations)
        for (i in seq_len(iterations)) {
          g <- (smu[s] + sqrt(sweep_omega[s]) * stats::rnorm(K)) * ginv
          f <- safe_fit(g)
          a1[i] <- f$est[1]; b1[i] <- f$se[1]
        }
        ok <- !is.na(a1) & !is.na(b1)
        tt <- a1[ok] + b1[ok] * stats::rnorm(sum(ok))
        sw[s, ] <- c(sde[s], smu[s], .pct_interval(tt, level), .pct_interval(a1[ok], level))
        .msg(verbose, "  sweep ", s, "/", ns, ": omega = ", format(sweep_omega[s]),
             "  CI [", format(round(sw[s, 3], 6)), ", ", format(round(sw[s, 4], 6)), "]")
      }
      res$sweep <- as.data.frame(sw)
    }
  }
  class(res) <- "apexog"
  res
}
