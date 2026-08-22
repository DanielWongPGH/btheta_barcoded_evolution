#!/bin/bash
# Export annotated isolate breseq calls in the flat layout required by Step 2.
# Usage: bash step1b_export_annotated_gds.sh <project_config.yaml> [output_dir]

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
config_file="${1:-}"
output_dir_override="${2:-}"
if [[ -z "$config_file" ]]; then
    echo "Usage: bash step1b_export_annotated_gds.sh <project_config.yaml> [output_dir]" >&2
    exit 1
fi

RAW_PROCESSING_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
source "${RAW_PROCESSING_DIR}/parse_config.sh" "$config_file" isolates

output_dir="${output_dir_override:-$ANNOTATED_GD_DIR}"
mkdir -p "$output_dir"
count=0
while IFS=$'\t' read -r barcode mouse day clone; do
    [[ "$mouse" == "mouse" ]] && continue
    clone_id="m${mouse}_day${day}_clone${clone}"
    mapfile -t matches < <(find "$BRESEQ_DIR" -mindepth 1 -maxdepth 1 -type d -name "*${clone_id}_*")
    if (( ${#matches[@]} != 1 )); then
        echo "Error: expected one breseq directory for ${clone_id}, found ${#matches[@]}" >&2
        exit 1
    fi
    source_gd="${matches[0]}/output/output.gd"
    if [[ ! -f "$source_gd" ]]; then
        echo "Error: missing annotated breseq GD: $source_gd" >&2
        exit 1
    fi
    cp "$source_gd" "${output_dir}/${clone_id}_annotated.gd"
    count=$((count + 1))
done < "$SCRIPT_DIR/clone_barcode.tsv"

echo "Exported ${count} annotated isolate GD files to ${output_dir}"
