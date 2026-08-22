#!/usr/bin/bash
# step1_breseq_genomes.sh
# Orchestrator script: submits breseq jobs for isolate genome analysis
#
# Usage: bash step1_breseq_genomes.sh <project_config.yaml> [samples_file]
#   project_config.yaml  - Pipeline configuration file
#   samples_file - Optional: file with samples (default: samples_in_dir.txt)
#                  First line is FASTQ directory, subsequent lines are sample names

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

config_file=$1
samples_file=${2:-"${SCRIPT_DIR}/samples_in_dir.txt"}

if [[ -z "$config_file" ]]; then
    echo "Usage: bash step1_breseq_genomes.sh <project_config.yaml> [samples_file]"
    exit 1
fi

# Load configuration
RAW_PROCESSING_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
source "${RAW_PROCESSING_DIR}/parse_config.sh" "$config_file" isolates

load_isolate_modules "$MODULE_BIOLOGY" "$MODULE_SAMTOOLS" "$MODULE_SYSTEM" \
    "$MODULE_NCURSES" "$MODULE_BOWTIE2" "$MODULE_NUMPY"

# Create output directories
create_output_dirs

# Get trimmed FASTQ path (assumes trimming was done, appends _trimmed to original dir)
# Or use TRIMMED_DIR from config
sample_path=$(head -n 1 "$samples_file")
if [[ -d "${sample_path}_trimmed" ]]; then
    sample_path="${sample_path}_trimmed"
elif [[ -d "$TRIMMED_DIR" ]]; then
    sample_path="$TRIMMED_DIR"
fi

echo "Input FASTQ dir: $sample_path"
echo "Breseq output dir: $BRESEQ_DIR"
echo "Reference: $REFERENCE_GBK"

# Process each sample
for sample_name in $(tail -n +2 "$samples_file"); do
    # Check if already processed
    if [[ -e "${BRESEQ_DIR}/${sample_name}/output/evidence/evidence.gd" ]]; then
        echo "${sample_name} already processed, skipping"
        continue
    fi

    if [[ -e "${BRESEQ_DIR}/${sample_name}" ]]; then
        echo "Error: incomplete breseq output exists for ${sample_name}; remove it after inspection." >&2
        exit 1
    fi
    mkdir -p "${BRESEQ_DIR}/${sample_name}"

    if [[ "$EXECUTION_SCHEDULER" == "slurm" ]]; then
        echo "Submitting breseq job for: $sample_name"
        sbatch \
            --job-name="breseq-${sample_name}" \
            --partition="${SLURM_PARTITION}" \
            --time="${SLURM_BRESEQ_TIME}" \
            --mem="${SLURM_BRESEQ_MEM}" \
            -o "${SLURM_LOG_DIR}/${sample_name}_breseq.out" \
            -e "${SLURM_LOG_DIR}/${sample_name}_breseq.err" \
            "${SCRIPT_DIR}/step1_breseq_genomes.sbatch" \
            "$config_file" "$sample_path" "$sample_name"
    elif [[ "$EXECUTION_SCHEDULER" == "local" ]]; then
        bash "${SCRIPT_DIR}/step1_breseq_genomes.sbatch" \
            "$config_file" "$sample_path" "$sample_name"
    else
        echo "Error: execution.scheduler must be slurm or local" >&2
        exit 1
    fi
done

echo "All breseq jobs submitted"
