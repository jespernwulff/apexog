#' apexog: plausibly exogenous bounds for average partial effects
#'
#' R port of the Stata command \code{apexog}: Conley-Hansen-Rossi (2012)
#' union-of-confidence-intervals and local-to-zero sensitivity bounds for the
#' average partial effect (or average interaction effect) of an endogenous
#' regressor in a binary-outcome model, estimated by a control-function probit
#' with Murphy-Topel standard errors or by a recursive bivariate probit.  The
#' main function is \code{\link{apexog}}.
#' @keywords internal
#' @importFrom stats coef vcov confint nobs
"_PACKAGE"
