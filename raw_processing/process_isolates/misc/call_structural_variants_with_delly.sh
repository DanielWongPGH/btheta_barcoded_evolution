#!/usr/bin/bash
# call_structural_variants_with_delly.sh
# Orchestrator script: submits delly jobs for structural variant calling
#
# Usage: bash call_structural_variants_with_delly.sh <project_config.yaml> <samples_file>
#   project_config.yaml  - Pipeline configuration file
#   samples_file - File with samples (first line is header, subsequent lines are sample names)

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

config_file=$1
samples_file=$2

if [[ -z "$config_file" || -z "$samples_file" ]]; then
    echo "Usage: bash call_structural_variants_with_delly.sh <project_config.yaml> <samples_file>"
    exit 1
fi

# Load configuration
RAW_PROCESSING_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
source "${RAW_PROCESSING_DIR}/parse_config.sh" "$config_file" isolates

# Create output directories
create_output_dirs

# Path to BAM file within breseq output
path_to_bam="data/reference.bam"

echo "Breseq input dir: $BRESEQ_DIR"
echo "Delly output dir: $DELLY_DIR"

# Process each sample
for sample_name in $(tail -n +2 "$samples_file"); do
    sample_path_to_bam="${BRESEQ_DIR}/${sample_name}/${path_to_bam}"

    echo "Sample: $sample_name"
    echo "BAM path: $sample_path_to_bam"

    # Check if breseq output exists
    if [[ ! -e "$sample_path_to_bam" ]]; then
        echo "WARNING: ${sample_name} breseq BAM does not exist! Skipping..."
        echo
        continue
    fi

    # Check if already processed
    if [[ -e "${DELLY_DIR}/${sample_name}.bcf" ]]; then
        echo "${sample_name} already processed, skipping"
        echo
        continue
    fi

    echo "Submitting delly job for: $sample_name"

    sbatch \
        --job-name="delly-${sample_name}" \
        --partition="${SLURM_PARTITION}" \
        --time="${SLURM_DELLY_TIME}" \
        --mem="${SLURM_DELLY_MEM}" \
        -o "${SLURM_LOG_DIR}/delly_${sample_name}.out" \
        -e "${SLURM_LOG_DIR}/delly_${sample_name}.err" \
        "${SCRIPT_DIR}/call_structural_variants_with_delly.sbatch" \
        "$config_file" "$sample_name" "$sample_path_to_bam"
    echo
done

echo "All delly jobs submitted"
