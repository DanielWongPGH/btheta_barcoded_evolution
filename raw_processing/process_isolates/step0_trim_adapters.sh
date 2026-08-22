#!/usr/bin/bash
# step0_trim_adapters.sh
# Orchestrator script: submits adapter trimming jobs for isolate FASTQs
#
# Usage: bash step0_trim_adapters.sh <project_config.yaml> <untrimmed_fastq_dir>
#   project_config.yaml         - Pipeline configuration file
#   untrimmed_fastq_dir - Directory containing untrimmed FASTQ files

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

config_file=$1
untrimmed_fastq_dir=$2

if [[ -z "$config_file" || -z "$untrimmed_fastq_dir" ]]; then
    echo "Usage: bash step0_trim_adapters.sh <project_config.yaml> <untrimmed_fastq_dir>"
    exit 1
fi

# Load configuration
RAW_PROCESSING_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
source "${RAW_PROCESSING_DIR}/parse_config.sh" "$config_file" isolates

load_isolate_modules "$MODULE_PYTHON"

# Strip trailing slash
untrimmed_fastq_dir="${untrimmed_fastq_dir%/}"

# Create output directories
create_output_dirs

echo "Untrimmed FASTQ dir: $untrimmed_fastq_dir"
echo "Trimmed output dir: $TRIMMED_DIR"
echo "SLURM logs: $SLURM_LOG_DIR"

# Generate sample list from FASTQ directory
samples_file="${SCRIPT_DIR}/samples_in_dir.txt"
"$PYTHON_BIN" "${SCRIPT_DIR}/parse_fastq_directory.py" \
    "$untrimmed_fastq_dir" "$samples_file" \
    --clone-barcode "${SCRIPT_DIR}/clone_barcode.tsv"

num_samples=$(wc -l < "$samples_file")
echo "Found $((num_samples - 1)) samples to process"
if (( num_samples <= 1 )); then
    echo "No paired isolate samples found"
    exit 0
fi

if [[ "$EXECUTION_SCHEDULER" == "slurm" ]]; then
    sbatch --array=2-${num_samples} \
        --partition="${SLURM_PARTITION}" \
        --time="${SLURM_TRIM_TIME}" \
        --mem-per-cpu="${SLURM_TRIM_MEM}" \
        --cpus-per-task="${SLURM_TRIM_CPUS}" \
        -o "${SLURM_LOG_DIR}/trim-%j.out" \
        -e "${SLURM_LOG_DIR}/trim-%j.err" \
        "${SCRIPT_DIR}/step0_trim_adapters.sbatch" \
        "$config_file" "$samples_file"
elif [[ "$EXECUTION_SCHEDULER" == "local" ]]; then
    for (( task_id=2; task_id<=num_samples; task_id++ )); do
        SLURM_ARRAY_TASK_ID="$task_id" bash "${SCRIPT_DIR}/step0_trim_adapters.sbatch" \
            "$config_file" "$samples_file"
    done
else
    echo "Error: execution.scheduler must be slurm or local" >&2
    exit 1
fi

echo "Submitted trimming jobs"
