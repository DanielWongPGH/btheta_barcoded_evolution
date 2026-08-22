#!/bin/bash
# step3_sbatch_noise_estimation.sh
# Orchestrator script: submits SLURM jobs for noise estimation
#
# Usage: bash step3_sbatch_noise_estimation.sh <project_config.yaml> <lib_id_fastq_map>
#   project_config.yaml       - Pipeline configuration file
#   lib_id_fastq_map  - File with library info (first line is FASTQ directory path)

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

config_file=$1
lib_id_fastq_map=$2

if [[ -z "$config_file" || -z "$lib_id_fastq_map" ]]; then
    echo "Usage: bash step3_sbatch_noise_estimation.sh <project_config.yaml> <lib_id_fastq_map>"
    exit 1
fi

# Load configuration
RAW_PROCESSING_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
source "${RAW_PROCESSING_DIR}/parse_config.sh" "$config_file" barseq

# Get fastq_dir from first row of lib_id_fastq_map
fastq_dir=$(head -1 "$lib_id_fastq_map")
fastq_dir="${fastq_dir%/}"

num_libraries=$(wc -l < "$lib_id_fastq_map")

echo "Processing FASTQs in $fastq_dir"
echo "Output directory: $OUTPUT_DIR"

# Create output directories
create_output_dirs

for (( i=2; i<=num_libraries; i++ )); do
    line=$(sed -n "${i}p" "$lib_id_fastq_map")

    env=$(echo "$line" | awk '{print $1}')
    population=$(echo "$line" | awk '{print $2}')
    timepoint=$(echo "$line" | awk '{print $3}')

    concat_info="${env}_${population}_${timepoint}"
    if [[ "$EXECUTION_SCHEDULER" == "slurm" ]]; then
        echo "Submitting job for: $concat_info"
        sbatch --job-name="$concat_info" \
            --partition="${SLURM_PARTITION}" \
            --mem="${SLURM_MEMORY}" \
            --time="${SLURM_TIME}" \
            "${SCRIPT_DIR}/step3_srun_noise_estimate.sh" \
            "$config_file" "$fastq_dir" "$lib_id_fastq_map" "$i"
    elif [[ "$EXECUTION_SCHEDULER" == "local" ]]; then
        echo "Running locally: $concat_info"
        bash "${SCRIPT_DIR}/step3_srun_noise_estimate.sh" \
            "$config_file" "$fastq_dir" "$lib_id_fastq_map" "$i"
    else
        echo "Error: execution.scheduler must be slurm or local" >&2
        exit 1
    fi
done

echo "All jobs submitted."
