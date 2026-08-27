#!/usr/bin/env bash
set -euo pipefail

project_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
wasi_sdk="${WASI_SDK:-}"
if [[ -z "$wasi_sdk" ]]; then
  echo "Set WASI_SDK to an extracted WASI SDK directory (WASI SDK 25 or newer)." >&2
  exit 2
fi
compiler="$wasi_sdk/bin/clang++"
if [[ ! -x "$compiler" ]]; then
  echo "No executable clang++ was found at $compiler." >&2
  exit 2
fi
wasm_opt="$project_dir/node_modules/.bin/wasm-opt"
if [[ ! -x "$wasm_opt" ]]; then
  echo "Run npm install before building so the pinned Binaryen wasm-opt is available." >&2
  exit 2
fi

build_dir="$project_dir/.build/alivibe-msa"
mkdir -p "$build_dir" "$project_dir/public"
export TMPDIR="$build_dir"
export LD_LIBRARY_PATH="$wasi_sdk/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"

"$compiler" --target=wasm32-wasip1 -O3 -flto -msimd128 -mbulk-memory -DNDEBUG \
  -std=c++20 -fno-exceptions -fno-rtti "$project_dir/wasm/vendor/alivibe_msa.cpp" \
  -mexec-model=reactor -Wl,--lto-O3 -Wl,--gc-sections -Wl,--export-memory -Wl,--strip-all \
  -o "$build_dir/alivibe-msa.raw.wasm"
"$wasm_opt" "$build_dir/alivibe-msa.raw.wasm" -O4 --converge --enable-bulk-memory \
  --enable-simd --strip-debug --strip-producers -o "$build_dir/alivibe-msa.wasm"
install -m 755 "$build_dir/alivibe-msa.wasm" "$project_dir/public/alivibe-msa.wasm"
echo "Built optimized Alivibe-compatible public/alivibe-msa.wasm"
