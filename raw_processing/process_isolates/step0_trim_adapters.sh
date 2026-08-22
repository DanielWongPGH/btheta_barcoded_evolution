#!/usr/bin/bash
# step0_trim_adapters.sh
# Orchestrator script: submits adapter trimming jobs for isolate FASTQs
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
source "${RAW_PROCESSING_DIR}/parse_config.sh" "$config_file" isolates

# Load modules needed for parsing
ml load ${MODULE_R} ${MODULE_BIOLOGY} ${MODULE_SAMTOOLS} ${MODULE_SYSTEM} ${MODULE_NCURSES}

# Strip trailing slash
untrimmed_fastq_dir="${untrimmed_fastq_dir%/}"

# Create output directories
create_output_dirs

echo "Untrimmed FASTQ dir: $untrimmed_fastq_dir"
echo "Trimmed output dir: $TRIMMED_DIR"
echo "SLURM logs: $SLURM_LOG_DIR"

# Generate sample list from FASTQ directory
samples_file="${SCRIPT_DIR}/samples_in_dir.txt"
python3 "${SCRIPT_DIR}/parse_fastq_directory.py" "$untrimmed_fastq_dir" "$samples_file"

num_samples=$(wc -l < "$samples_file")
echo "Found $((num_samples - 1)) samples to process"

# Submit array job for trimming
sbatch --array=2-${num_samples} \
    --partition="${SLURM_PARTITION}" \
    --time="${SLURM_TRIM_TIME}" \
    --mem-per-cpu="${SLURM_TRIM_MEM}" \
    --cpus-per-task="${SLURM_TRIM_CPUS}" \
    -o "${SLURM_LOG_DIR}/trim-%j.out" \
    -e "${SLURM_LOG_DIR}/trim-%j.err" \
    "${SCRIPT_DIR}/step0_trim_adapters.sbatch" \
    "$config_file" "$samples_file"

echo "Submitted trimming jobs"
