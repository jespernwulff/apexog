# Hex logo for apexog: a probability curve inside a shaded band, the picture of
# a partial effect bounded over a range of assumed violations.  Writes
# man/figures/logo.svg.  Run from the package root: Rscript data-raw/logo.R
teal <- "#0e5e5a"; rust <- "#b85c2c"; amber <- "#e2a35a"; paper <- "#fbf8f2"

hex <- "259,8 510,153 510,447 259,592 8,447 8,153"
t <- seq(-3.2, 3.2, length.out = 120)
xs <- 259 + t * 48
ymap <- function(p) 392 - 236 * p
pts <- function(x, y) paste(sprintf("%.1f,%.1f", x, y), collapse = " ")
mid <- ymap(stats::pnorm(t))
hi <- ymap(stats::pnorm(t + 0.62))
lo <- ymap(stats::pnorm(t - 0.62))
band <- paste(pts(xs, hi), pts(rev(xs), rev(lo)))

svg <- c(
  '<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 518 600" width="518" height="600">',
  '  <title>apexog</title>',
  sprintf('  <polygon points="%s" fill="%s" stroke="%s" stroke-width="14" stroke-linejoin="round"/>',
          hex, teal, rust),
  sprintf('  <line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" stroke-opacity="0.35" stroke-width="3"/>',
          min(xs), ymap(0), max(xs), ymap(0), paper),
  sprintf('  <line x1="259" y1="%.1f" x2="259" y2="%.1f" stroke="%s" stroke-opacity="0.55" stroke-width="3" stroke-dasharray="10 9"/>',
          ymap(1.05), ymap(-0.05), paper),
  sprintf('  <polygon points="%s" fill="%s" fill-opacity="0.75"/>', band, amber),
  sprintf('  <polyline points="%s" fill="none" stroke="%s" stroke-width="9" stroke-linecap="round"/>',
          pts(xs, mid), paper),
  sprintf('  <circle cx="259" cy="%.1f" r="13" fill="%s" stroke="%s" stroke-width="5"/>',
          ymap(0.5), rust, paper),
  sprintf('  <text x="259" y="498" text-anchor="middle" font-family="Georgia, \'Times New Roman\', serif" font-size="84" font-weight="700" fill="%s" letter-spacing="1">apexog</text>',
          paper),
  '</svg>')
dir.create(file.path("man", "figures"), showWarnings = FALSE, recursive = TRUE)
writeLines(svg, file.path("man", "figures", "logo.svg"))
cat("wrote man/figures/logo.svg\n")
