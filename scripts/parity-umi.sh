#!/usr/bin/env bash
set -euo pipefail

repo_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
julia_bin="${JULIA_BIN:-julia}"
julia_project="${PORPID_JULIA_PROJECT:-}"
if [[ -z "$julia_project" ]]; then
  echo "Set PORPID_JULIA_PROJECT to the supplied PORPIDpipeline source tree." >&2
  exit 2
fi

work_dir="$(mktemp -d)"
trap 'rm -rf "$work_dir"' EXIT
export TMPDIR="$work_dir"
"${CXX:-c++}" -std=c++20 -O3 -I"$repo_dir/wasm/include" \
  "$repo_dir/tests/cpp/umi_parity_driver.cpp" \
  "$repo_dir/wasm/src/core.cpp" "$repo_dir/wasm/src/preprocess.cpp" "$repo_dir/wasm/src/umi.cpp" \
  "$repo_dir/wasm/src/alignment.cpp" "$repo_dir/wasm/src/consensus.cpp" \
  -o "$work_dir/umi-parity"
"$work_dir/umi-parity" > "$work_dir/native.tsv"
(cd "$julia_project" && "$julia_bin" --startup-file=no --project=. "$repo_dir/tests/julia/umi_parity.jl") > "$work_dir/julia.tsv"
node "$repo_dir/scripts/compare-umi-parity.mjs" "$work_dir/native.tsv" "$work_dir/julia.tsv"
