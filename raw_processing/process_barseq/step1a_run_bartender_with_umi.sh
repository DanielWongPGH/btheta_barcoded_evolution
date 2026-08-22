#!/bin/bash
#SBATCH -p normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=00:30:00

# Worker: merge R1/R2, extract UMIs, and run Bartender.
#
# Usage: Called by step1_sbatch_bartender.sh -UMI true
#   $1 - project_config.yaml
#   $2 - fastq_dir
#   $3 - sample name

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

config_file=$1
fastq_dir=$2
sample=$3

# Load configuration
RAW_PROCESSING_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
source "${RAW_PROCESSING_DIR}/parse_config.sh" "$config_file" barseq

echo "Processing sample: $sample"
echo "FASTQ dir: $fastq_dir"
echo "Merged FASTQ dir: $MERGED_FASTQ_DIR"
echo "Bartender dir: $BARTENDER_DIR"

echo "Merging FASTQs..."
echo "$PYTHON_BIN" "${SCRIPT_DIR}/merge_barseq_fastq.py" "$fastq_dir" "$sample" "$MERGED_FASTQ_DIR"
run_barseq_command "$PYTHON_BIN" "${SCRIPT_DIR}/merge_barseq_fastq.py" "$fastq_dir" "$sample" "$MERGED_FASTQ_DIR"

load_configured_modules "$MODULE_GCC"
export PATH=${BARTENDER_PATH}:$PATH

echo "Running bartender extractor..."
run_barseq_command "$BARTENDER_EXTRACTOR" \
    -f "${MERGED_FASTQ_DIR}/merged_${sample}.fastq" \
    -o "${BARTENDER_DIR}/${sample}" \
    -u ${BARTENDER_UMI_RANGE} \
    -p "${BARTENDER_PATTERN}" \
    -q "${BARTENDER_QUALITY}"

echo "Running bartender clustering..."
run_barseq_command "$BARTENDER_CLUSTER" \
    -f "${BARTENDER_DIR}/${sample}_barcode.txt" \
    -o "${BARTENDER_DIR}/${sample}" \
    -z ${BARTENDER_SEED} \
    -d ${BARTENDER_DISTANCE}

echo "Done processing $sample"
