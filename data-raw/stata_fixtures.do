*==============================================================================
* stata_fixtures.do
*
* Reference results for the R package's test suite.  Runs a set of apexog
* specifications covering every engine / target / option on ONE simulated data
* set and exports, per specification, every e() matrix (svmat double, exact),
* every e() scalar and every e() macro to .dta files in
*   data-raw/stata_fixtures/
* data-raw/make_fixtures.R converts them into the .rds file the R tests read
* (so R CMD check needs neither Stata nor haven).
*
* LTZ specifications store e(gamma_sim): the R tests refit at Stata's own gamma
* draws, which removes the RNG difference between Stata and R.
*
* Requires the Stata command apexog (net install apexog, from(...); see
* https://github.com/jespernwulff/stata-apexog).  Run in Stata from the root of
* the R package:
*   do data-raw/stata_fixtures.do
* The fixtures shipped with the package were written by apexog 0.20.1 in
* StataNow 19.5.
*==============================================================================
clear all
set more off
global FXDIR "data-raw/stata_fixtures"
capture mkdir "$FXDIR"
which apexog

*------------------------------------------------------------------------------
* fxsave TAG: dump e() matrices, scalars and macros of the last apexog call
*------------------------------------------------------------------------------
capture program drop fxsave
program fxsave
    args tag
    local mats : e(matrices)
    local k = 0
    foreach m of local mats {
        local ++k
        tempname M`k'
        matrix `M`k'' = e(`m')
        local mname`k' "`m'"
    }
    forvalues j = 1/`k' {
        preserve
        clear
        qui svmat double `M`j'', names(c)
        qui save "$FXDIR/fx_`tag'_m_`mname`j''.dta", replace
        restore
    }
    local scal : e(scalars)
    tempname h
    postfile `h' str40 name double value using "$FXDIR/fx_`tag'_scalars.dta", replace
    foreach s of local scal {
        post `h' ("`s'") (e(`s'))
    }
    postclose `h'
    local mac : e(macros)
    postfile `h' str40 name str2045 value using "$FXDIR/fx_`tag'_macros.dta", replace
    foreach s of local mac {
        post `h' ("`s'") (`"`e(`s')'"')
    }
    postclose `h'
    di as result "  fixture `tag': `k' matrices saved"
end

*------------------------------------------------------------------------------
* Data
*------------------------------------------------------------------------------
set seed 20261107
set obs 3000
gen double x1   = rnormal()
gen byte   x2   = runiform() < 0.5
gen double x3   = rnormal()
gen byte   cat  = 1 + floor(3*runiform())
gen double z1   = rnormal()
gen double z2   = rnormal()
gen byte   zb   = runiform() < 0.5
gen double v    = rnormal()
gen double v2   = rnormal()
gen double e    = 0.5*v + sqrt(0.75)*rnormal()
gen int    team = ceil(_n/15)
* continuous endog
gen double d    = 0.7*z1 + 0.4*z2 + 0.3*x1 + 0.2*x2 + v
gen byte   y    = (-0.3 + 0.5*d + 0.4*x1 + 0.3*x2 + 0.2*(cat==2) - 0.2*(cat==3) + e) > 0
* binary endog
gen byte   db   = (-0.2 + 0.9*z1 + 0.3*zb + 0.3*x1 + v) > 0
gen byte   yb   = (-0.2 + 0.6*db + 0.4*x1 + 0.3*x2 + e) > 0
* fractional endog
gen double df   = normal(0.6*z1 + 0.3*x1 + 0.5*v)
gen byte   yf   = (-0.3 + 1.0*df + 0.4*x1 + e) > 0
* count endog
gen double dc   = rpoisson(exp(0.2 + 0.4*z1 + 0.2*x1 + 0.3*v))
gen byte   yc   = (-0.5 + 0.3*dc + 0.4*x1 + e) > 0
* second continuous endog
gen double d2   = 0.6*z2 + 0.4*zb + 0.2*x1 + v2 + 0.3*v
gen byte   ym   = (-0.2 + 0.5*db + 0.3*d2 + 0.3*x1 + e) > 0
gen byte   ym2  = (-0.2 + 0.4*d + 0.3*d2 + 0.3*x1 + e) > 0
* interaction designs
gen byte   yx   = (-0.2 + 0.3*x1 + 0.3*x3 + 0.3*x1*x3 + 0.5*d + e) > 0
gen byte   yx2  = (-0.2 + 0.4*x1 + 0.3*x3 + 0.5*d + 0.4*d*x3 + e) > 0
gen byte   yi   = (-0.2 + 0.4*x1 + 0.3*x2 + 0.5*d + 0.4*d*x2 + e) > 0
gen byte   ybb  = (-0.2 + 0.3*x1 + 0.3*x2 + 0.5*db + 0.4*db*x2 + e) > 0
* a control with missing values
gen double x1m  = x1
replace    x1m  = . if runiform() < 0.03
* string cluster id
gen str6   teams = "T" + string(team)
* a constant control: Stata omits it, so the fit equals a no-controls model
* (apexog's parser cannot take an empty control list; the R port can)
gen byte   one  = 1

save "$FXDIR/fx_data.dta", replace

*------------------------------------------------------------------------------
* (A) linear first stage
*------------------------------------------------------------------------------
apexog uci y x1 x2 (d = z1), gmin(-0.1) gmax(0.1) grid(5)
fxsave A1

apexog ltz y x1 x2 (d = z1), mu(0) omega(0.0025) iterations(100) seed(101)
fxsave A2

apexog uci y x1 x2 i.cat (d = z1 z2), gmin(-0.05) gmax(0.05) grid(3)
fxsave A3

apexog uci y x1m x2 (d = z1) if x3 > -1.5, gmin(-0.05) gmax(0.05) grid(3)
fxsave A5

apexog uci y one (d = z1), gmin(-0.1) gmax(0.1) grid(3)
fxsave A6

*------------------------------------------------------------------------------
* (B) nonlinear first stages, multi-endog, cluster, cfinteract
*------------------------------------------------------------------------------
apexog uci yb x1 x2 (db = z1 zb), firststage(probit) gmin(-0.1) gmax(0.1) grid(3)
fxsave B1

apexog uci yf x1 (df = z1), firststage(fprobit) gmin(-0.1) gmax(0.1) grid(5)
fxsave B2

apexog uci yc x1 (dc = z1), firststage(poisson) gmin(-0.1) gmax(0.1) grid(5)
fxsave B3

apexog uci ym x1 (db = z1 zb) (d2 = z2 zb), firststage(probit linear) gmin(-0.05) gmax(0.05) grid(3)
fxsave B4

apexog uci y x1 x2 (d = z1), vce(cluster team) gmin(-0.1) gmax(0.1) grid(3)
fxsave B5

apexog uci ym2 x1 (d d2 = z1 z2 zb), cfinteract gmin(-0.05) gmax(0.05) grid(2)
fxsave B6

apexog uci y x1 x2 (d = z1), vce(cluster teams) gmin(-0.1) gmax(0.1) grid(3)
fxsave B7

* cfinteract with a 0/1 regressor (discrete-change APE; v0.20.1 Jacobian fix)
apexog uci ym x1 (db = z1 zb) (d2 = z2 zb), firststage(probit linear) cfinteract gmin(-0.05) gmax(0.05) grid(2)
fxsave B8

*------------------------------------------------------------------------------
* (C) AIE targets
*------------------------------------------------------------------------------
* exog cts x cts
apexog uci yx x1 x3 (d = z1), aie(x1 x3) gmin(-0.1) gmax(0.1) grid(3)
fxsave C1
* cts endog d x binary exog x2 via aie() (product as exog control)
apexog uci yi x1 x2 (d = z1), aie(d x2) gmin(-0.1) gmax(0.1) grid(3)
fxsave C2
* binary exog x2 x binary endog db (2nd discrete difference)
apexog uci ybb x1 x2 (db = z1 zb), firststage(probit) aie(x2 db) gmin(-0.1) gmax(0.1) grid(3)
fxsave C3
* binary exog x2 x cts endog d (firstdiff in x2)
apexog uci yi x1 x2 (d = z1), aie(x2 d) gmin(-0.1) gmax(0.1) grid(3)
fxsave C3b
* endogenous interaction cts x cts, auto cfinteract
apexog uci yx2 x1 x3 (c.d c.d#c.x3 = z1 c.z1#c.x3), target(c.d#c.x3) gmin(-0.05) gmax(0.05) grid(3)
fxsave C4
* endogenous interaction cts x binary
apexog uci yi x1 x2 (c.d c.d#i.x2 = z1 c.z1#i.x2), target(c.d#i.x2) gmin(-0.05) gmax(0.05) grid(3)
fxsave C5
* forced cts formula on a 0/1 moderator
apexog uci yi x1 x2 (c.d c.d#c.x2 = z1 c.z1#c.x2), target(c.d#c.x2) gmin(-0.05) gmax(0.05) grid(3)
fxsave C6
* APE subset via target()
apexog uci ym x1 (db = z1 zb) (d2 = z2 zb), firststage(probit linear) target(d2) gmin(-0.05) gmax(0.05) grid(3)
fxsave C7

*------------------------------------------------------------------------------
* (D) recursive bivariate probit engine
*------------------------------------------------------------------------------
apexog uci yb x1 x2 (db = z1 zb), firststage(biprobit) gmin(-0.1) gmax(0.1) grid(3)
fxsave D1

apexog uci yb x1 x2 (db = z1 zb), firststage(biprobit) vce(cluster team) gmin(0) gmax(0) grid(2)
fxsave D2

*------------------------------------------------------------------------------
* (E) gunits, priors, LTZ with AIE / biprobit
*------------------------------------------------------------------------------
apexog uci yb x1 x2 (db = z1 zb), firststage(probit) gunits(pr) gmin(-0.02) gmax(0.02) grid(3)
fxsave E1

apexog ltz y x1 x2 (d = z1), gunits(structural) mu(0.01) omega(0.0004) iterations(60) seed(102)
fxsave E2

apexog ltz y x1 x2 (d = z1 z2), distribution(uniform -0.1 0.1 -0.05 0.05) iterations(60) seed(103)
fxsave E3

apexog ltz yx2 x1 x3 (c.d c.d#c.x3 = z1 c.z1#c.x3), target(c.d#c.x3) mu(0) omega(0.0004) iterations(40) seed(104)
fxsave E4

apexog ltz yb x1 x2 (db = z1 zb), firststage(biprobit) gunits(pr) mu(0) omega(0.0001) iterations(30) seed(105)
fxsave E5

apexog uci y x1 x2 (d = z1), gunits(structural) gmin(-0.1) gmax(0.1) grid(3) level(90)
fxsave E6

di as result "ALL FIXTURES WRITTEN"
