#!/usr/bin/env bash
# Build the pkgdown site OUTSIDE Dropbox (Dropbox and the virus scanner lock
# freshly written files, which breaks pkgdown's write-then-rewrite of each page)
# and copy the result into ./docs.  Keep the build root short: Windows' 260-
# character path limit is hit by the bundled Google-font files otherwise.
# The package is installed into a private library inside the build folder,
# never into the user library.
#
# Usage (Git Bash, from anywhere):  bash r/apexog/tools/build_site.sh [build_root]
set -euo pipefail
SRC="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
ROOT="${1:?usage: build_site.sh <short build root outside Dropbox, e.g. a scratch folder>}"
BUILD="$ROOT/$(date +%H%M%S)"
RBIN="/c/Program Files/R/R-4.2.1/bin/x64"
mkdir -p "$BUILD/lib" "$BUILD/apexog"
(cd "$SRC" && tar cf - --exclude=./docs --exclude=./data-raw/stata_fixtures .) |
  (cd "$BUILD/apexog" && tar xf -)
BW="$(cygpath -m "$BUILD")"
"$RBIN/R.exe" CMD INSTALL --no-test-load --library="$BW/lib" "$BW/apexog" \
  > "$BUILD/install.log" 2>&1
# R must NOT run from inside the build copy: pkgdown 2.0.6 then builds an
# over-long path to pkgdown/pre-reference.R and the reference step fails.
"$RBIN/Rscript.exe" -e ".libPaths(c('$BW/lib', .libPaths()));
  pkgdown::build_site('$BW/apexog', new_process = FALSE, install = FALSE,
                      preview = FALSE)" 2>&1 | tee "$BUILD/build.log" | tail -5
grep -q "Error" "$BUILD/build.log" && { echo "BUILD FAILED: $BUILD/build.log"; exit 1; }
# pkgdown 2.0.6 resolves each article's "Source:" path against the working
# directory instead of the package; point the links back at vignettes/.
for f in "$BUILD"/apexog/docs/articles/*.html; do
  sed -i -E 's#blob/HEAD/[^"]*/vignettes/#blob/HEAD/vignettes/#; s#<code>[^<]*/vignettes/#<code>vignettes/#' "$f"
done
touch "$BUILD/apexog/docs/.nojekyll"   # GitHub Pages: serve the files as they are
mkdir -p "$SRC/docs"
(cd "$BUILD/apexog/docs" && tar cf - .) | (cd "$SRC/docs" && tar xf -)
echo "site built in $BUILD and copied to $SRC/docs"
