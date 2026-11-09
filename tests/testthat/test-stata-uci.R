# R vs Stata apexog v0.21.0, method UCI.  Tolerances are relative and sit at
# optimizer tolerance (Stata's ml stopping rules vs glm.fit / the Newton loop
# here); the observed agreement is about 1e-7 or better on estimates and SEs.
uci_tags <- c("A1", "A3", "A5", "A6", "B1", "B2", "B3", "B4", "B5", "B6", "B7", "B8",
              "C1", "C2", "C3", "C3b", "C4", "C5", "C6", "C7",
              "D1", "D2", "E1", "E6")

for (tag in uci_tags) {
  test_that(paste("UCI spec", tag, "matches Stata"), {
    sp <- fx$specs[[tag]]
    f <- fx_calls[[tag]]()
    K <- length(f$instruments)
    expect_equal(f$N, unname(sp$s[["N"]]))
    expect_equal(K, length(strsplit(sp$mac[["instruments"]], " ")[[1]]))
    expect_equal(f$baseline$estimate, drop(sp$m$baseline_ape_vec), tolerance = 1e-6)
    expect_equal(f$baseline$se, drop(sp$m$baseline_se_vec), tolerance = 1e-6)
    expect_equal(unname(f$vcov), sp$m$V, tolerance = 1e-6)
    expect_equal(f$mphi, unname(sp$s[["mphi"]]), tolerance = 1e-6)
    expect_equal(f$stilde, unname(sp$s[["stilde"]]), tolerance = 1e-6)
    g <- unname(as.matrix(f$grid))
    expect_equal(g[, seq_len(K), drop = FALSE], sp$m$grid[, seq_len(K), drop = FALSE],
                 tolerance = 1e-8)
    expect_equal(g[, -seq_len(K), drop = FALSE], sp$m$grid[, -seq_len(K), drop = FALSE],
                 tolerance = 1e-5)
    expect_equal(unname(f$grid_pr), sp$m$grid_pr, tolerance = 1e-6)
    lb <- sp$s[grep("^lb_", names(sp$s))]
    ub <- sp$s[grep("^ub_", names(sp$s))]
    expect_equal(sort(f$bounds$lower), sort(unname(unique(lb))), tolerance = 1e-6)
    expect_equal(sort(f$bounds$upper), sort(unname(unique(ub))), tolerance = 1e-6)
    if ("rho" %in% names(sp$s)) expect_equal(unname(f$rho), unname(sp$s[["rho"]]), tolerance = 1e-6)
  })
}
