#!/usr/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1

# step1_breseq_genomes.sbatch
# Worker script: runs breseq on a single isolate sample
#
# Usage: Called by step1_breseq_genomes.sh
#   $1 - project_config.yaml
#   $2 - input_dir (trimmed FASTQ directory)
#   $3 - sample_name

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

config_file=$1
input_dir=$2
sample_name=$3

# Load configuration
RAW_PROCESSING_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
source "${RAW_PROCESSING_DIR}/parse_config.sh" "$config_file" isolates

echo "Processing sample: $sample_name"
echo "Input dir: $input_dir"
echo "Output dir: ${BRESEQ_DIR}/${sample_name}"
echo "Reference: $REFERENCE_GBK"

# Collect all FASTQ files for this sample
FASTQ_FILES=""
for f in $(ls ${input_dir}/${sample_name}*.fastq.gz 2>/dev/null); do
    FASTQ_FILES="${FASTQ_FILES} ${f}"
done

if [[ -z "$FASTQ_FILES" ]]; then
    echo "ERROR: No FASTQ files found for $sample_name in $input_dir"
    exit 1
fi

echo "FASTQ files: $FASTQ_FILES"

# Run breseq
echo "Running breseq..."
${BRESEQ_PATH} \
    --polymorphism-reject-indel-homopolymer-length ${BRESEQ_POLY_INDEL_HOMOPOLY} \
    --polymorphism-reject-surrounding-homopolymer-length ${BRESEQ_POLY_SURROUND_HOMOPOLY} \
    --polymorphism-score-cutoff ${BRESEQ_POLY_SCORE_CUTOFF} \
    --polymorphism-minimum-variant-coverage-each-strand ${BRESEQ_POLY_MIN_COV_STRAND} \
    -o "${BRESEQ_DIR}/${sample_name}" \
    -r "${REFERENCE_GBK}" \
    ${FASTQ_FILES} \
    > "${BRESEQ_DIR}/${sample_name}.out" \
    2> "${BRESEQ_DIR}/${sample_name}.err"

echo "Done processing $sample_name"
