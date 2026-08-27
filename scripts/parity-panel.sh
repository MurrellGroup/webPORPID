#!/usr/bin/env bash
set -euo pipefail

repo_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
julia_bin="${JULIA_BIN:-julia}"
julia_project="${PORPID_JULIA_PROJECT:-}"
if [[ -z "$julia_project" ]]; then
  echo "Set PORPID_JULIA_PROJECT to the attached PORPIDpipeline source tree." >&2
  exit 2
fi

work_dir="$(mktemp -d)"
trap 'rm -rf "$work_dir"' EXIT
export TMPDIR="$work_dir"
node "$repo_dir/tests/panel-parity-driver.mjs" > "$work_dir/native.tsv"
(cd "$julia_project" && "$julia_bin" --startup-file=no --project=. "$repo_dir/tests/julia/panel_parity.jl") > "$work_dir/julia.tsv"
node "$repo_dir/scripts/compare-panel-parity.mjs" "$work_dir/native.tsv" "$work_dir/julia.tsv"
