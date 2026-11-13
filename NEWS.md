# apexog 0.21.1

* Fractional-probit first stage (`firststage = "fprobit"`): the first-stage
  weight, minus the derivative of the generalized residual in the first-stage
  index, now uses the formula for a fractional regressor,
  `chat * (wpi + chat) + d * (1 - d) * (phi / (Phi * (1 - Phi)))^2`.  Versions
  through 0.21.0 used the binary-probit identity `chat * (wpi + chat)`, which
  holds only for a 0/1 regressor.  Point estimates are unchanged; the mean
  standard error with a fractional-probit first stage moves by at most about
  3%, in either direction, in the validation design, because the two places
  the weight enters largely offset each other.  All other first stages, and a fractional-probit first
  stage on a 0/1 regressor, are unchanged.  Mirrors Stata apexog 0.21.1.
* New tests check the fractional-probit and Poisson first stages against a
  stacked two-step estimator differentiated numerically in every parameter.

# apexog 0.21.0

* The standard error of every two-step effect (APE and AIE, all first stages,
  `cfinteract`, cluster-robust) now adds the direct dependence of the effect
  on the estimated first-stage parameters through the control function: the
  influence function is `J_beta' psi_beta + J_pi' psi_pi`, where versions
  through 0.20.1 used only the first term (as does Stata's `cfprobit` followed
  by `margins`).  Point estimates are unchanged; standard errors move by about
  1% in the validation designs, in either direction.  The bivariate-probit
  engine is unaffected.  Mirrors Stata apexog 0.21.0.

# apexog 0.20.1

First public release of the R package, a port of the Stata command `apexog`
0.20.1 that gives the same results.

* `apexog()` with `method = "uci"` (union of confidence intervals over a range
  for the violation) and `method = "ltz"` (local-to-zero interval under a
  prior: normal, uniform, chi-squared, t, gamma, Poisson, an empirical donor,
  a function, or user-supplied draws), with a sweep over prior widths.
* Control-function probit with a linear, probit, fractional-probit or Poisson
  first stage and Murphy-Topel two-step standard errors; recursive bivariate
  probit for one binary treatment.
* Average partial effects of every endogenous regressor, or an average
  interaction effect (continuous by continuous, continuous by binary, binary
  by binary), including an endogenous treatment interacted with a moderator.
* The violation on the latent index, in probability points or on the
  structural latent scale, with the exact probability-scale value of every
  refit.
* Robust or cluster-robust standard errors; several instruments and several
  endogenous regressors.
* Example data `promotion` and `certification`, shipped also as Stata files in
  `inst/extdata`.
* Since Stata apexog 0.20.1 the standard error of an average partial effect
  with `cfinteract` includes the sampling variance of the control-function
  product terms; the R package does the same.
