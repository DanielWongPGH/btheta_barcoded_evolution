#!/bin/bash
# Submit Bartender jobs with or without UMI correction.
#
# Usage: bash step1_sbatch_bartender.sh <project_config.yaml> <sample_list> [-UMI [true|false]]

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

if (( $# < 2 || $# > 4 )); then
    echo "Usage: bash step1_sbatch_bartender.sh <project_config.yaml> <sample_list> [-UMI [true|false]]" >&2
    exit 1
fi

config_file="$1"
sample_list="$2"
use_umi=false

if (( $# >= 3 )); then
    if [[ "$3" != "-UMI" ]]; then
        echo "Error: unknown argument: $3" >&2
        exit 1
    fi
    use_umi="${4:-true}"
fi

case "$use_umi" in
    true|TRUE|True|1|yes|YES|Yes) use_umi=true ;;
    false|FALSE|False|0|no|NO|No) use_umi=false ;;
    *)
        echo "Error: -UMI must be true or false" >&2
        exit 1
        ;;
esac

if [[ ! -f "$sample_list" ]]; then
    echo "Error: sample list not found: $sample_list" >&2
    exit 1
fi

# Load configuration
RAW_PROCESSING_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
source "${RAW_PROCESSING_DIR}/parse_config.sh" "$config_file" barseq

# Override fastq_dir from sample_list header (first row is directory path)
fastq_dir="$(head -1 "$sample_list")"

# Strip trailing slashes
fastq_dir="${fastq_dir%/}"
OUTPUT_DIR="${OUTPUT_DIR%/}"

echo "FASTQ directory: $fastq_dir"
echo "Output directory: $OUTPUT_DIR"

# Create output directories
create_output_dirs

# Count libraries (subtract 1 for header)
num_libraries=$(( $(wc -l < "$sample_list") - 1 ))

echo "Processing $num_libraries samples from $fastq_dir (UMI=$use_umi)"

if [[ "$use_umi" == true ]]; then
    worker="${SCRIPT_DIR}/step1a_run_bartender_with_umi.sh"
    log_suffix="UMI"
else
    worker="${SCRIPT_DIR}/step1b_run_bartender_no_umi.sh"
    log_suffix="noUMI"
fi

for (( i=2; i<=num_libraries+1; i++ )); do
    sample="$(sed -n "${i}p" "$sample_list")"

    if [[ "$EXECUTION_SCHEDULER" == "slurm" ]]; then
        echo "Submitting job for: $sample"
        sbatch --partition="${SLURM_PARTITION}" \
            --mem="${SLURM_MEMORY}" \
            --time="${SLURM_TIME}" \
            -e "${SLURM_LOG_DIR}/${sample}_${log_suffix}.err" \
            -o "${SLURM_LOG_DIR}/${sample}_${log_suffix}.out" \
            "$worker" \
            "$config_file" "$fastq_dir" "$sample"
    elif [[ "$EXECUTION_SCHEDULER" == "local" ]]; then
        echo "Running locally: $sample"
        bash "$worker" "$config_file" "$fastq_dir" "$sample"
    else
        echo "Error: execution.scheduler must be slurm or local" >&2
        exit 1
    fi
done

echo "All Bartender jobs submitted (UMI=$use_umi)."
