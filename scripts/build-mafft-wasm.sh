#!/usr/bin/env bash
set -euo pipefail

MAFFT_COMMIT=52b59f064c600da59bca8233736418fb8bb35d5e
REPOSITORY_ROOT=$(cd "$(dirname "$0")/.." && pwd)
PATCH_FILE="$REPOSITORY_ROOT/wasm/patches/mafft-7.520-wasm.patch"
OUTPUT_DIRECTORY="$REPOSITORY_ROOT/public/biowasm/mafft"
TEMPORARY_PARENT=${TMPDIR:-/tmp}
if [[ ! -d "$TEMPORARY_PARENT" ]]; then TEMPORARY_PARENT="$REPOSITORY_ROOT/.build"; mkdir -p "$TEMPORARY_PARENT"; fi
BUILD_ROOT=$(mktemp -d "$TEMPORARY_PARENT/webporpid-mafft.XXXXXXXX")
trap 'rm -rf "$BUILD_ROOT"' EXIT

if ! command -v emcc >/dev/null 2>&1; then
  echo "emcc is required (the release asset was built with Emscripten 3.1.64)." >&2
  exit 1
fi

SOURCE_DIRECTORY="$BUILD_ROOT/mafft"
if [[ -n "${MAFFT_SOURCE_DIR:-}" ]]; then
  mkdir -p "$SOURCE_DIRECTORY"
  git -C "$MAFFT_SOURCE_DIR" archive "$MAFFT_COMMIT" | tar -x -C "$SOURCE_DIRECTORY"
else
  git clone --filter=blob:none https://gitlab.com/sysimm/mafft.git "$SOURCE_DIRECTORY"
  git -C "$SOURCE_DIRECTORY" checkout --detach "$MAFFT_COMMIT"
fi
git -C "$SOURCE_DIRECTORY" apply "$PATCH_FILE"
mkdir -p "$SOURCE_DIRECTORY/build" "$OUTPUT_DIRECTORY"

COMMON_FLAGS="-O3 -flto -msimd128 -Denablemultithread"
LINK_FLAGS="-O3 -flto -msimd128 -sUSE_PTHREADS=0 -sALLOW_MEMORY_GROWTH=1 -sINITIAL_MEMORY=67108864 -sMAXIMUM_MEMORY=4294967296 -sSTACK_SIZE=2097152 -sMODULARIZE=1 -sEXPORT_ES6=1 -sEXPORT_NAME=createMafftDisttbfast -sENVIRONMENT=web,worker,node -sEXPORTED_RUNTIME_METHODS=FS,callMain -sINVOKE_RUN=0 -sEXIT_RUNTIME=0 -sFILESYSTEM=1 -sASSERTIONS=0"

make -C "$SOURCE_DIRECTORY/core" disttbfast CC=emcc CFLAGS="$COMMON_FLAGS" LDFLAGS="$LINK_FLAGS" LIBS=-lm
cp "$SOURCE_DIRECTORY/build/disttbfast.mjs" "$OUTPUT_DIRECTORY/disttbfast.mjs"
cp "$SOURCE_DIRECTORY/build/disttbfast.wasm" "$OUTPUT_DIRECTORY/disttbfast.wasm"
cp "$SOURCE_DIRECTORY/license" "$OUTPUT_DIRECTORY/LICENSE"

echo "Built MAFFT 7.520 FFT-NS-2 WebAssembly assets in $OUTPUT_DIRECTORY"
sha256sum "$OUTPUT_DIRECTORY/LICENSE" "$OUTPUT_DIRECTORY/disttbfast.mjs" "$OUTPUT_DIRECTORY/disttbfast.wasm"
