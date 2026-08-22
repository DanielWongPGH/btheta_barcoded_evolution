#!/bin/bash
# Collect annotated breseq calls for the day-zero metagenomic libraries.
# Usage: bash step1b_collect_day0_annotated_gds.sh <project_config.yaml>

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
config_file="${1:-}"
if [[ -z "$config_file" ]]; then
    echo "Usage: bash step1b_collect_day0_annotated_gds.sh <project_config.yaml>" >&2
    exit 1
fi

RAW_PROCESSING_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
source "${RAW_PROCESSING_DIR}/parse_config.sh" "$config_file" metagenomics

mkdir -p "$DAY0_ANNOTATED_GD_DIR"
count=0
while IFS=, read -r sample_id population timepoint extraction_batch extraction_tube flagged directory filename_stem; do
    [[ "$timepoint" == "0" ]] || continue
    source_gd="${BRESEQ_DIR}/${filename_stem}/output/output.gd"
    output_gd="${DAY0_ANNOTATED_GD_DIR}/${filename_stem}_annotated.gd"
    if [[ ! -f "$source_gd" ]]; then
        echo "Error: missing day-zero breseq output: $source_gd" >&2
        exit 1
    fi
    cp "$source_gd" "$output_gd"
    count=$((count + 1))
done < <(tail -n +2 "${SCRIPT_DIR}/population_samples.csv")

echo "Collected ${count} day-zero GD files in ${DAY0_ANNOTATED_GD_DIR}"
