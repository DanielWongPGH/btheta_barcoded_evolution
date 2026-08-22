#!/bin/bash
# Set isolate manifest headers from project_config.yaml without changing samples.
# Usage: bash generate_fastq_stems.sh <project_config.yaml> [manifest ...]

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
config_file="${1:-}"

if [[ -z "$config_file" || ! -f "$config_file" ]]; then
    echo "Usage: bash generate_fastq_stems.sh <project_config.yaml> [manifest ...]" >&2
    exit 1
fi

if ! command -v yq >/dev/null 2>&1; then
    echo "Error: mikefarah/yq version 4 is required." >&2
    exit 1
fi

shift
if (( $# > 0 )); then
    manifests=("$@")
else
    manifests=("${SCRIPT_DIR}/E1_isolates.txt")
fi

for manifest in "${manifests[@]}"; do
    if [[ ! -f "$manifest" ]]; then
        echo "Error: manifest not found: $manifest" >&2
        exit 1
    fi

    key="$(basename "$manifest" .txt)"
    fastq_dir="$(yq -r ".raw_processing.isolates.manifest_fastq_dirs.\"${key}\" // .raw_processing.isolates.fastq_dir // \"\"" "$config_file")"
    fastq_dir="${fastq_dir%/}"

    if [[ -z "$fastq_dir" || "$fastq_dir" == "null" || "$fastq_dir" == /path/to/* ]]; then
        echo "Error: configure raw_processing.isolates.manifest_fastq_dirs.${key} in $config_file" >&2
        exit 1
    fi

    temp_manifest="$(mktemp "${manifest}.tmp.XXXXXX")"
    {
        printf '%s\n' "$fastq_dir"
        tail -n +2 "$manifest"
    } > "$temp_manifest"
    mv "$temp_manifest" "$manifest"
    echo "Updated $manifest -> $fastq_dir"
done
