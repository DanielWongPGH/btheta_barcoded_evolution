#!/usr/bin/bash
# step0_trim_adapters.sh
# Orchestrator script: submits adapter trimming jobs for metagenomics FASTQs
#
# Usage: bash step0_trim_adapters.sh <project_config.yaml> <untrimmed_fastq_dir>
#   project_config.yaml         - Pipeline configuration file
#   untrimmed_fastq_dir - Directory containing untrimmed FASTQ files

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

config_file=$1
untrimmed_fastq_dir=$2

if [[ -z "$config_file" || -z "$untrimmed_fastq_dir" ]]; then
    echo "Usage: bash step0_trim_adapters.sh <project_config.yaml> <untrimmed_fastq_dir>"
    exit 1
fi

# Load configuration
RAW_PROCESSING_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
source "${RAW_PROCESSING_DIR}/parse_config.sh" "$config_file" metagenomics

# Load modules needed for parsing
load_common_modules

# Strip trailing slash
untrimmed_fastq_dir="${untrimmed_fastq_dir%/}"

# Create output directories
create_output_dirs

# Set trimmed output directory
trimmed_fastq_dir="${untrimmed_fastq_dir}_trimmed"
mkdir -p "$trimmed_fastq_dir"

echo "Untrimmed FASTQ dir: $untrimmed_fastq_dir"
echo "Trimmed output dir: $trimmed_fastq_dir"
echo "SLURM logs: $SLURM_LOG_DIR"

# Generate sample list from FASTQ directory
samples_file="${SCRIPT_DIR}/samples_in_dir.txt"
"${PYTHON3_BIN}" "${SCRIPT_DIR}/parse_mgx_files.py" "$untrimmed_fastq_dir" "$samples_file"

num_samples=$(wc -l < "$samples_file")
echo "Found $((num_samples - 1)) samples to process"

# Submit array job for trimming
submit_metagenomics_job --array=2-${num_samples} \
    --partition="${SLURM_PARTITION}" \
    --time="${SLURM_TRIM_TIME}" \
    --mem-per-cpu="${SLURM_TRIM_MEM}" \
    --cpus-per-task="${SLURM_TRIM_CPUS}" \
    --output="${SLURM_LOG_DIR}/trim-%A_%a.log" \
    --error="${SLURM_LOG_DIR}/trim-%A_%a.log" \
    "${SCRIPT_DIR}/step0_trim_adapters.sbatch" \
    "$config_file" "$samples_file" "$trimmed_fastq_dir"

echo "Submitted trimming jobs"
