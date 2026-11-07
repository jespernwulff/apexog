# apexog

<div class="apx-hero">
<p class="apx-kicker">R package &middot; Stata command</p>
<p class="apx-question">How large would a violation of the exclusion restriction have to be before the average partial effect from your probit changes sign, or stops being significant?</p>
<p>apexog answers that question for instrumental-variables models with a binary outcome. It carries the plausibly exogenous framework of Conley, Hansen and Rossi (2012) from linear 2SLS coefficients to average partial effects and average interaction effects on the probability scale, estimated by a control-function probit or by a recursive bivariate probit.</p>
<div class="apx-buttons">
<a class="btn btn-light" href="articles/apexog.html">Get started in R</a>
<a class="btn btn-light" href="articles/stata.html">Get started in Stata</a>
<a class="btn btn-outline-light" href="articles/gamma-scale.html">Which scale for the violation?</a>
</div>
</div>

## What it does

An instrument *z* is valid if it moves the endogenous regressor *d* and has no direct effect on the outcome. That second condition, the exclusion restriction, cannot be tested, and in applied work it is rarely beyond doubt. Conley, Hansen and Rossi (2012) replace it by an assumption about how large the direct effect, *γ*, can be:

$$y^* = x'\beta + d\,\delta + z'\gamma + \varepsilon, \qquad y = 1[y^* > 0].$$

With *γ* = 0 this is the usual model. apexog asks what the average partial effect (APE) of *d* on the probability that *y* = 1 would be for other values of *γ*. At each assumed violation it refits the model with *z'γ* as an offset in the latent index and recomputes the APE and its standard error: two-step (Murphy and Topel 1985) standard errors for the control-function probit, maximum-likelihood standard errors for the bivariate probit. The results are summarized in one of two ways:

- **Union of confidence intervals** (`uci`). You give a range for *γ*. The union of the confidence intervals over that range covers the APE with at least the nominal probability as long as the true *γ* lies in the range.
- **Local to zero** (`ltz`). You give a prior for *γ*. Each draw from the prior is refitted, and the reported interval combines the spread that the prior induces with the sampling uncertainty of each refit.

## Why a probit version?

The existing implementation of Conley, Hansen and Rossi's methods, `plausexog` (Clarke and Matta 2018), works with 2SLS. In a linear model a direct effect *z'γ* can simply be subtracted from the outcome before re-estimating. With a binary outcome that is not possible, because the violation acts on the latent index, not on the observed 0/1 outcome. Running 2SLS on the binary outcome turns the problem into a linear probability model, with a different target and a violation stated on a different scale. apexog keeps the nonlinear model: the violation shifts the latent index, and the target remains the average partial effect on the probability scale. You can still state the violation in probability points, the scale `plausexog` uses for a binary outcome, so the two can be compared for the same belief.

## Install {.tabset .tabset-pills}

### R

```r
# install.packages("remotes")
remotes::install_github("jespernwulff/apexog")
```

### Stata

```stata
net install apexog, from("https://raw.githubusercontent.com/jespernwulff/stata-apexog/main/")
```

The Stata command requires Stata 18 or later. It implements the control-function probit itself and does not need StataNow's `cfprobit`.

## One call {.tabset .tabset-pills}

The example data are simulated. Employees choose hours of training partly on an unobserved ability that also drives promotion, so training is endogenous. A training subsidy shifts training. Could the subsidy also affect promotion directly, for instance because it signals that the firm values the employee? Suppose that direct effect is at most 0.10 in absolute value on the latent index.

### R

```r
library(apexog)
data(promotion)

apexog(promoted ~ tenure + manager, endog = training ~ subsidy,
       data = promotion, method = "uci", gmin = -0.10, gmax = 0.10, grid = 11)
#> ------------------------------------------------------------------------------
#> Plausibly-exogenous APE bounds (CHR 2012, UCI method)
#> Estimator: cfprobit  |  First stage: linear  |  SE: Murphy-Topel  |  N = 4000
#> Endogenous: training
#> Target:    APE of training
#> ------------------------------------------------------------------------------
#> Baseline (gamma = 0):
#>     APE(training)      = 0.144457  (SE 0.008323)  [0.128145, 0.160769]
#>
#> Scale:     1 latent-index unit of gamma = 0.2091 pr. points per unit of Z (avg., baseline)
#>            1 structural-latent unit of gamma = 1.1986 latent-index units  (sigma_tilde = 0.8343)
#>
#> UCI union over gamma grid (11 point(s) per instrument):
#>     gamma_subsidy  range: [-0.1000, 0.1000]   exact Pr scale: [-0.0209, 0.0209]
#>
#>     APE(training)      bounds = [0.097219, 0.190375]
#> ------------------------------------------------------------------------------
```

### Stata

```stata
use "https://raw.githubusercontent.com/jespernwulff/apexog/main/inst/extdata/promotion.dta", clear

apexog uci promoted tenure manager (training = subsidy), gmin(-0.10) gmax(0.10) grid(11)
```

```{.stata-output}
------------------------------------------------------------------------------
Plausibly-exogenous APE bounds (CHR 2012, UCI method)
Estimator: cfprobit  |  First stage: linear  |  SE: Murphy-Topel  |  N = 4000
Endogenous: training
Target:    APE of every endogenous regressor (default)
------------------------------------------------------------------------------
Baseline (gamma = 0):
    APE(training)              =  0.144457  (SE 0.008323)  [ 0.128145,  0.160769]

Scale:     1 latent-index unit of gamma = 0.2091 pr. points per unit of Z (avg., baseline)
           1 structural-latent unit of gamma = 1.1986 latent-index units  (sigma_tilde = 0.8343)

UCI union over gamma grid (11 point(s) per IV):
    gamma_subsidy range  : [-0.1000,  0.1000]   exact Pr scale: [-0.0209,  0.0209]

    APE(training) bounds       = [ 0.097219,  0.190375]
------------------------------------------------------------------------------
```

## Reading the result

Without a violation, an extra unit of training raises the probability of promotion by 0.144 on average (standard error 0.008). If the subsidy's direct effect could be anywhere between −0.10 and 0.10 on the latent index (about ±2 percentage points on the probability scale, as the "exact Pr scale" line shows), the effect lies between 0.097 and 0.190 with at least 95% confidence. It stays clearly positive. R and Stata analyse the same file and print the same numbers.

## What you get

<div class="apx-cards">
<div class="apx-card"><h4>Two estimators</h4><p>A control-function probit with a linear, probit, fractional-probit or Poisson first stage and two-step standard errors, or, for one binary treatment, the recursive bivariate probit.</p></div>
<div class="apx-card rust"><h4>Partial and interaction effects</h4><p>The APE of every endogenous regressor, or an average interaction effect: continuous by continuous, continuous by binary, binary by binary, including an endogenous treatment interacted with a moderator.</p></div>
<div class="apx-card"><h4>A range or a prior</h4><p>The union of confidence intervals over a range for the violation, or the local-to-zero interval under a Gaussian, uniform, other or user-supplied prior, with a sweep over prior widths.</p></div>
<div class="apx-card rust"><h4>Three scales for the violation</h4><p>State the violation on the latent index, in probability points, or on the structural latent scale, and see its exact effect on the probability scale at every refit.</p></div>
<div class="apx-card"><h4>The inference you need</h4><p>Robust or cluster-robust standard errors, several instruments and several endogenous regressors.</p></div>
<div class="apx-card rust"><h4>One method, two languages</h4><p>The R package is a port of the Stata command and is tested against it on 29 specifications; on the same data the two agree to about seven significant digits.</p></div>
</div>

## Where to go from here

- [Get started in R](articles/apexog.html) and [in Stata](articles/stata.html): a first analysis of the example data, the output line by line.
- [Which scale for the violation?](articles/gamma-scale.html) How to state *γ*, and how the scales relate.
- [Union of intervals or a prior?](articles/uci-ltz.html) Choosing between the two methods and their settings.
- [Binary treatments and interaction effects](articles/binary-and-interactions.html): the control function or the bivariate probit, and interaction effects.
- [R and Stata side by side](articles/r-and-stata.html): every option and stored result, matched.

## Citation

If you use apexog, please cite the method it implements and the software:

Conley, T. G., Hansen, C. B., & Rossi, P. E. (2012). Plausibly exogenous. *Review of Economics and Statistics*, 94(1), 260–272. <https://doi.org/10.1162/REST_a_00139>

Wulff, J. N. (2026). *apexog: Plausibly exogenous bounds for average partial effects in binary-outcome IV models.* R package and Stata command, version 0.20.1. <https://jespernwulff.github.io/apexog/>

The two-step standard errors follow Murphy, K. M., & Topel, R. H. (1985). Estimation and inference in two-step econometric models. *Journal of Business & Economic Statistics*, 3(4), 370–379. <https://doi.org/10.1080/07350015.1985.10509471>. The linear-model version of the method is `plausexog`: Clarke, D., & Matta, B. (2018). Practical considerations for questionable IVs. *Stata Journal*, 18(3), 663–691. <https://doi.org/10.1177/1536867X1801800308>.
