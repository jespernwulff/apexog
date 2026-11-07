*==============================================================================
* site_examples.do -- the Stata runs whose output is shown on the package site
* (articles stata, binary-and-interactions, gamma-scale, uci-ltz, r-and-stata).
* Uses the shipped example data in ../inst/extdata; the site loads the same
* files from GitHub.  Run from this folder with apexog >= 0.20.1 on the adopath.
* The Stata graph is exported to ../man/figures (copied into the site).
*==============================================================================
clear all
set more off
set linesize 120
version 18
which apexog

*--- 1. UCI, continuous treatment ------------------------------------------------
use "../inst/extdata/promotion.dta", clear
describe, short
apexog uci promoted tenure manager (training = subsidy), gmin(-0.10) gmax(0.10) grid(11)
lincom training

*--- 2. graph over a wider range -------------------------------------------------
apexog uci promoted tenure manager (training = subsidy), gmin(-0.20) gmax(0.20) grid(21) graph
graph export "../man/figures/stata-uci.png", replace width(1400)

*--- 3. LTZ, Gaussian prior ------------------------------------------------------
apexog ltz promoted tenure manager (training = subsidy), mu(0) omega(0.0025) iterations(1000) seed(101)

*--- 4. probability-point units --------------------------------------------------
apexog uci promoted tenure manager (training = subsidy), gunits(pr) gmin(-0.02) gmax(0.02) grid(11)

*--- 5. clustered standard errors ------------------------------------------------
apexog uci promoted tenure manager (training = subsidy), gmin(-0.10) gmax(0.10) grid(11) vce(cluster team)

*--- 6. interaction effect -------------------------------------------------------
apexog uci promoted tenure manager (c.training c.training#c.manager = subsidy c.subsidy#i.manager), target(c.training#i.manager) gmin(-0.05) gmax(0.05) grid(5)

*--- 7. binary treatment: control function vs bivariate probit -------------------
use "../inst/extdata/certification.dta", clear
apexog uci promoted tenure manager (i.certified = subsidy), gmin(-0.10) gmax(0.10) grid(11)
apexog uci promoted tenure manager (certified = subsidy), firststage(biprobit) gmin(-0.10) gmax(0.10) grid(11)
