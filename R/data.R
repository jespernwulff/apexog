#' Simulated promotions: a continuous endogenous treatment
#'
#' Employees choose hours of training partly on an unobserved ability that also
#' drives promotion, so training is endogenous.  A randomly assigned training
#' subsidy shifts training; in this simulation it has no direct effect on
#' promotion, so the true exclusion-restriction violation is zero.  Training
#' pays off more under a supportive manager.  The same file is shipped as
#' \code{inst/extdata/promotion.dta} for the Stata command.
#'
#' @format A data frame with 4,000 rows and 6 variables:
#' \describe{
#'   \item{promoted}{Promoted (0/1); the outcome.}
#'   \item{training}{Hours of training (standardized); the endogenous treatment.}
#'   \item{subsidy}{Training subsidy (standardized); the instrument.}
#'   \item{tenure}{Tenure (standardized); a control.}
#'   \item{manager}{Supportive manager (0/1); a control and moderator.}
#'   \item{team}{Team identifier, 200 teams of 20; for clustered standard errors.}
#' }
#' @details Data-generating process (\code{data-raw/example_data.do}, seed
#'   20260613): \code{training = 0.7 subsidy + 0.3 tenure + 0.2 manager + a},
#'   \code{promoted = 1[-0.3 + 0.5 training + 0.4 tenure + 0.3 manager +
#'   0.4 training x manager + e > 0]}, with \code{a ~ N(0, 1)} the unobserved
#'   ability and \code{e = 0.5 a + sqrt(0.75) u}, \code{u ~ N(0, 1)}.
#' @source Simulated.
#' @examples
#' data(promotion)
#' summary(promotion)
"promotion"

#' Simulated promotions: a binary endogenous treatment
#'
#' Employees obtain a certification partly on an unobservable that is
#' correlated (0.5) with the promotion error, so certification is endogenous.
#' A randomly assigned certification subsidy shifts certification and has no
#' direct effect on promotion.  The errors are jointly normal, the case in
#' which the recursive bivariate probit (\code{firststage = "biprobit"}) is
#' exact and the generalized-residual control function only approximate.  The
#' same file is shipped as \code{inst/extdata/certification.dta}.
#'
#' @format A data frame with 4,000 rows and 5 variables:
#' \describe{
#'   \item{promoted}{Promoted (0/1); the outcome.}
#'   \item{certified}{Certified (0/1); the endogenous treatment.}
#'   \item{subsidy}{Certification subsidy (standardized); the instrument.}
#'   \item{tenure}{Tenure (standardized); a control.}
#'   \item{manager}{Supportive manager (0/1); a control.}
#' }
#' @details Data-generating process (\code{data-raw/example_data.do}, seed
#'   90210): \code{certified = 1[-0.2 + 0.9 subsidy + 0.3 tenure + 0.2 manager
#'   + v > 0]}, \code{promoted = 1[-0.2 + 0.6 certified + 0.4 tenure +
#'   0.3 manager + e > 0]}, with \code{(e, v)} standard bivariate normal with
#'   correlation 0.5.  The true average effect of certification on the
#'   probability of promotion is about 0.21.
#' @source Simulated.
#' @examples
#' data(certification)
#' table(certification$certified, certification$promoted)
"certification"
