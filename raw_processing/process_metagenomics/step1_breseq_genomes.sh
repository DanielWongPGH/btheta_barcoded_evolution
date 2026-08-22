#!/usr/bin/bash
# step1_breseq_genomes.sh
# Orchestrator script: submits breseq jobs for metagenomics samples
#
# Usage: bash step1_breseq_genomes.sh <project_config.yaml> <population> <sequence_type> [output_dir]
#   project_config.yaml   - Pipeline configuration file
#   population    - Population name(s), comma-separated (e.g., E1_1,E1_2)
#   sequence_type - Sample type (e.g., non_stragglers, stragglers)
#   output_dir    - Optional: override output directory from config

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

config_file=$1
population=$2
sequence_type=$3
output_dir_override=$4

if [[ -z "$config_file" || -z "$population" || -z "$sequence_type" ]]; then
    echo "Usage: bash step1_breseq_genomes.sh <project_config.yaml> <population> <sequence_type> [output_dir]"
    exit 1
fi

# Load configuration
RAW_PROCESSING_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
source "${RAW_PROCESSING_DIR}/parse_config.sh" "$config_file" metagenomics

# Use override or config default
if [[ -n "$output_dir_override" ]]; then
    output_dir="${output_dir_override%/}"
else
    output_dir="$BRESEQ_DIR"
fi

# Load required modules
load_breseq_modules

# Create output directories
create_output_dirs
mkdir -p "$output_dir"

echo "Population: $population"
echo "Sequence type: $sequence_type"
echo "Output dir: $output_dir"
echo "SLURM logs: $SLURM_LOG_DIR"

for pop in $(echo $population | tr "," "\n"); do
    for i in $("${PYTHON2_BIN}" "${SCRIPT_DIR}/population_parameters.py" samples ${pop} ${sequence_type}); do
        arrI=(${i//;/ })
        sample_path=${arrI[0]}
        sample_name=${arrI[1]}

        # Append _trimmed to sample path
        if [[ $sample_path == */ ]]; then
            sample_path=${sample_path%/}_trimmed
        else
            sample_path=${sample_path}_trimmed
        fi

        echo "Sample path: $sample_path"
        echo "Sample name: $sample_name"

        # Check if already processed
        if [[ -e "${output_dir}/${sample_name}/output/evidence/evidence.gd" ]]; then
            echo "${sample_name} already processed, skipping"
            continue
        fi

        echo "Submitting breseq job for: $sample_name"

        rm -rf "${output_dir}/${sample_name}"
        mkdir -p "${output_dir}/${sample_name}"

        submit_metagenomics_job \
            --job-name="breseq-${sample_name}" \
            --partition="${SLURM_PARTITION}" \
            --time="${SLURM_BRESEQ_TIME}" \
            --output="${SLURM_LOG_DIR}/breseq-%j.log" \
            --error="${SLURM_LOG_DIR}/breseq-%j.log" \
            "${SCRIPT_DIR}/step1_breseq_genomes.sbatch" \
            "$config_file" "$sample_path" "$sample_name" "$output_dir"
    done
done

echo "All breseq jobs submitted"
