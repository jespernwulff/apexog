# R vs Stata, method LTZ.  R refits at Stata's own gamma draws (e(gamma_sim)),
# so the per-draw estimates and SEs are compared one to one; the percentile
# rule is compared on Stata's convolved draws (e(ape_sim_cnv)).
for (tag in names(fx_ltz_conv)) {
  test_that(paste("LTZ spec", tag, "matches Stata at Stata's draws"), {
    sp <- fx$specs[[tag]]
    cv <- fx_ltz_conv[[tag]]
    conv <- if (is.character(cv)) 1 / sp$s[[cv]] else cv
    f <- fx_calls[[tag]](sp$m$gamma_sim / conv)
    expect_equal(unname(f$sims$gamma), sp$m$gamma_sim, tolerance = 1e-6)
    expect_equal(unname(f$sims$est), sp$m$ape_sim, tolerance = 1e-5)
    expect_equal(unname(f$sims$se), sp$m$se_sim, tolerance = 1e-5)
    expect_equal(unname(f$sims$gamma_pr), sp$m$gamma_pr_sim, tolerance = 1e-5)
    expect_equal(f$baseline$estimate, drop(sp$m$baseline_ape_vec), tolerance = 1e-6)
    pc <- apexog:::.pct_interval(sp$m$ape_sim_cnv[, 1], f$level)
    expect_identical(pc, c(sp$m$lo_pct_vec[1, 1], sp$m$hi_pct_vec[1, 1]))
    pr <- apexog:::.pct_interval(sp$m$ape_sim[, 1], f$level)
    expect_identical(pr, c(sp$m$lo_prior_vec[1, 1], sp$m$hi_prior_vec[1, 1]))
    # the R prior band uses R's own refits: equal up to optimizer tolerance
    expect_equal(c(f$prior_band$lower, f$prior_band$upper),
                 c(sp$m$lo_prior_vec[1, 1], sp$m$hi_prior_vec[1, 1]), tolerance = 1e-5)
  })
}
