#!/bin/bash
#SBATCH -p normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=00:30:00

# Worker: run Bartender directly on R1 without UMI correction.
#
# Usage: Called by step1_sbatch_bartender.sh -UMI false
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

echo "Processing sample (no UMI): $sample"
echo "FASTQ dir: $fastq_dir"
echo "Bartender dir: $BARTENDER_NOUMI_DIR"

# Load required modules
load_configured_modules "$MODULE_GCC"
export PATH=${BARTENDER_PATH}:$PATH

R1_FASTQ="${fastq_dir}/${sample}_R1_001.fastq"
if [[ ! -f "$R1_FASTQ" && -f "${R1_FASTQ}.gz" ]]; then
    R1_FASTQ="${MERGED_FASTQ_DIR}/${sample}_R1_001.fastq"
    if [[ ! -f "$R1_FASTQ" ]]; then
        gzip -cd "${fastq_dir}/${sample}_R1_001.fastq.gz" > "$R1_FASTQ"
    fi
fi
if [[ ! -f "$R1_FASTQ" ]]; then
    echo "Error: R1 FASTQ not found for $sample" >&2
    exit 1
fi

echo "Running bartender extractor (no UMI)..."
run_barseq_command "$BARTENDER_EXTRACTOR" \
    -f "$R1_FASTQ" \
    -o "${BARTENDER_NOUMI_DIR}/${sample}" \
    -p "${BARTENDER_PATTERN}" \
    -q "${BARTENDER_QUALITY}"

echo "Running bartender clustering..."
run_barseq_command "$BARTENDER_CLUSTER" \
    -f "${BARTENDER_NOUMI_DIR}/${sample}_barcode.txt" \
    -o "${BARTENDER_NOUMI_DIR}/${sample}" \
    -z ${BARTENDER_SEED} \
    -d ${BARTENDER_DISTANCE}

echo "Done processing $sample"
