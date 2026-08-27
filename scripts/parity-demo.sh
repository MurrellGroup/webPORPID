#!/usr/bin/env bash
set -euo pipefail

repo_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
julia_bin="${JULIA_BIN:-julia}"
julia_project="${PORPID_JULIA_PROJECT:-}"
if [[ -z "$julia_project" ]]; then
  echo "Set PORPID_JULIA_PROJECT to the attached PORPIDpipeline source tree." >&2; exit 2
fi
work_dir="$(mktemp -d)"; trap 'rm -rf "$work_dir"' EXIT
export TMPDIR="$work_dir"
(cd "$repo_dir" && "$julia_bin" --startup-file=no --project="$julia_project" tests/julia/demo_parity.jl demo/synthetic_reads.fastq) \
  | rg '^(counts|family|consensus)\t' > "$work_dir/julia.tsv"
node "$repo_dir/scripts/build-cli.mjs" >/dev/null
node "$repo_dir/cli/porpid-cli.mjs" run "$repo_dir/demo/synthetic_reads.fastq" --config "$repo_dir/demo/synthetic_config.yaml" \
  --output "$work_dir/result.webporpid" --workers 2 >/dev/null
node "$repo_dir/scripts/compare-demo-parity.mjs" "$work_dir/julia.tsv" "$work_dir/result.webporpid"
