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

build_dir="$project_dir/.build/wasm"
mkdir -p "$build_dir" "$project_dir/public"
export TMPDIR="$build_dir"
export LD_LIBRARY_PATH="$wasi_sdk/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"

sources=(
  "$project_dir/wasm/src/core.cpp"
  "$project_dir/wasm/src/preprocess.cpp"
  "$project_dir/wasm/src/umi.cpp"
  "$project_dir/wasm/src/alignment.cpp"
  "$project_dir/wasm/src/consensus.cpp"
  "$project_dir/wasm/src/web_adapter.cpp"
)

"$compiler" --target=wasm32-wasip1 -O3 -msimd128 -mbulk-memory -DNDEBUG \
  -std=c++20 -fno-exceptions -fno-rtti -I"$project_dir/wasm/include" "${sources[@]}" \
  -mexec-model=reactor -Wl,--gc-sections -Wl,--export-memory -Wl,--strip-all \
  -Wl,--export=wpp_alloc -Wl,--export=wpp_free -Wl,--export=wpp_version \
  -Wl,--export=wpp_init_config -Wl,--export=wpp_preprocess -Wl,--export=wpp_partition_counts \
  -Wl,--export=wpp_count_families -Wl,--export=wpp_build_family_model -Wl,--export=wpp_init_family_model \
  -Wl,--export=wpp_consensus_partition -Wl,--export=wpp_stats -Wl,--export=wpp_result_ptr \
  -Wl,--export=wpp_result_len -Wl,--export=wpp_error_ptr -Wl,--export=wpp_error_len \
  -o "$build_dir/webporpid.raw.wasm"

"$wasm_opt" "$build_dir/webporpid.raw.wasm" -O4 --converge --enable-bulk-memory \
  --enable-simd --enable-nontrapping-float-to-int --strip-debug --strip-producers -o "$build_dir/webporpid.wasm"
install -m 755 "$build_dir/webporpid.wasm" "$project_dir/public/webporpid.wasm"
echo "Built optimized SIMD public/webporpid.wasm"
