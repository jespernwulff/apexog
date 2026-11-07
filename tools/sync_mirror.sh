#!/usr/bin/env bash
# Copy the R package source (./r/apexog, the single source of truth) into its
# Git mirror ./github/apexog.  Files are added or overwritten, never deleted;
# the mirror's .git folder is not touched.
#
# Usage (Git Bash):  bash r/apexog/tools/sync_mirror.sh
set -euo pipefail
SRC="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
MIRROR="$(cd "$SRC/../.." && pwd)/github/apexog"
mkdir -p "$MIRROR"
(cd "$SRC" && tar cf - --exclude=./.git .) | (cd "$MIRROR" && tar xf -)
echo "synced $SRC -> $MIRROR"
