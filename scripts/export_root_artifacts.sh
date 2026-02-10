#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 1 ]]; then
  echo "Usage: $0 <root_file> [output_base_dir]" >&2
  exit 1
fi

ROOT_FILE="$1"
OUT_BASE="${2:-results/root/png}"

root -l -b -q "analysis/root/export_run_artifacts.C(\"${ROOT_FILE}\",\"${OUT_BASE}\")"
