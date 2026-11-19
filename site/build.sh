#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SITE="$ROOT/site"
QUARTZ_REF="6c85c0950ac0bd9567d2509dd1d50d38e4e6696c"
QUARTZ_DIR="${QUARTZ_DIR:-$SITE/.quartz}"

if [[ ! -d "$QUARTZ_DIR/.git" ]]; then
  git clone --filter=blob:none --no-checkout https://github.com/MathStruct/quartz.git "$QUARTZ_DIR"
fi

if ! git -C "$QUARTZ_DIR" cat-file -e "$QUARTZ_REF^{commit}" 2>/dev/null; then
  git -C "$QUARTZ_DIR" fetch --filter=blob:none origin "$QUARTZ_REF"
fi
git -C "$QUARTZ_DIR" checkout --detach "$QUARTZ_REF"

exec bash "$QUARTZ_DIR/mathstruct/site.sh" "$SITE" "$ROOT/vault" "$@"
