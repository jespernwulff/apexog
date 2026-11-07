*==============================================================================
* example_data.do -- the two simulated example data sets of apexog.
*
* Generated in Stata (so that Stata users and R users analyse the identical
* file); data-raw/example_data.R converts them to the package's .rda files.  Run from this
* folder:  stata /e do example_data.do
*
*   promotion      continuous endogenous treatment (training), instrument
*                  subsidy, controls tenure and manager, clusters team
*   certification  binary endogenous treatment (certified), instrument subsidy
*==============================================================================
clear all
set more off
version 18
capture mkdir "../inst"
capture mkdir "../inst/extdata"

*--- promotion -----------------------------------------------------------------
* Employees choose training partly on unobserved ability, which also drives
* promotion, so training is endogenous.  A randomly assigned training subsidy
* shifts training; in this simulation it has no direct effect on promotion
* (the true violation is zero).
set seed 20260613
set obs 4000
gen double tenure   = rnormal()
gen byte   manager  = runiform() < 0.5
gen double subsidy  = rnormal()
gen double ability  = rnormal()
gen double eps      = 0.5*ability + sqrt(0.75)*rnormal()
gen double training = 0.7*subsidy + 0.3*tenure + 0.2*manager + ability
gen double pstar    = -0.3 + 0.5*training + 0.4*tenure + 0.3*manager ///
                      + 0.4*training*manager + eps
gen byte   promoted = pstar > 0
gen int    team     = ceil(_n/20)
keep promoted training subsidy tenure manager team
order promoted training subsidy tenure manager team
label variable promoted "Promoted (0/1)"
label variable training "Hours of training (standardized; endogenous)"
label variable subsidy  "Training subsidy (standardized; instrument)"
label variable tenure   "Tenure (standardized)"
label variable manager  "Supportive manager (0/1)"
label variable team     "Team identifier (200 teams of 20)"
label data "apexog example data: simulated promotions (continuous treatment)"
save "../inst/extdata/promotion.dta", replace

*--- certification -------------------------------------------------------------
* A binary treatment: certification, chosen on an unobservable correlated with
* the promotion error (correlation 0.5).  The subsidy is again truly excluded.
clear
set seed 90210
set obs 4000
gen double tenure    = rnormal()
gen byte   manager   = runiform() < 0.5
gen double subsidy   = rnormal()
gen double v         = rnormal()
gen double e         = 0.5*v + sqrt(0.75)*rnormal()
gen byte   certified = (-0.2 + 0.9*subsidy + 0.3*tenure + 0.2*manager + v) > 0
gen double pstar     = -0.2 + 0.6*certified + 0.4*tenure + 0.3*manager + e
gen byte   promoted  = pstar > 0
keep promoted certified subsidy tenure manager
order promoted certified subsidy tenure manager
label variable promoted  "Promoted (0/1)"
label variable certified "Certified (0/1; endogenous)"
label variable subsidy   "Certification subsidy (standardized; instrument)"
label variable tenure    "Tenure (standardized)"
label variable manager   "Supportive manager (0/1)"
label data "apexog example data: simulated promotions (binary treatment)"
save "../inst/extdata/certification.dta", replace
