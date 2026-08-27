#!/usr/bin/env bash
set -euo pipefail
project_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
build_dir="$project_dir/.build/native-tests"
mkdir -p "$build_dir"
export TMPDIR="$build_dir"
"${CXX:-g++}" -std=c++20 -O2 -Wall -Wextra -Wpedantic -I"$project_dir/wasm/include" \
  "$project_dir/wasm/src/core.cpp" "$project_dir/wasm/src/preprocess.cpp" \
  "$project_dir/wasm/src/umi.cpp" "$project_dir/wasm/src/alignment.cpp" \
  "$project_dir/wasm/src/consensus.cpp" "$project_dir/tests/cpp/core_tests.cpp" \
  -o "$build_dir/core-tests"
"$build_dir/core-tests"
