#!/usr/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1

# step1_breseq_genomes.sbatch
# Worker script: runs breseq on a single metagenomics sample
#
# Usage: Called by step1_breseq_genomes.sh
#   $1 - project_config.yaml
#   $2 - input_dir (trimmed FASTQ directory)
#   $3 - sample_name
#   $4 - output_dir

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

config_file=$1
input_dir=$2
sample_name=$3
output_dir=$4

# Load configuration
RAW_PROCESSING_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
source "${RAW_PROCESSING_DIR}/parse_config.sh" "$config_file" metagenomics

echo "Processing sample: $sample_name"
echo "Input dir: $input_dir"
echo "Output dir: ${output_dir}/${sample_name}"
echo "Reference: $REFERENCE_GBK"

# Collect all FASTQ files for this sample
FASTQ_FILES=""
for f in $(ls ${input_dir}/${sample_name}*_001.fastq.gz 2>/dev/null); do
    FASTQ_FILES="${FASTQ_FILES} ${f}"
done

if [[ -z "$FASTQ_FILES" ]]; then
    echo "ERROR: No FASTQ files found for $sample_name in $input_dir"
    exit 1
fi

echo "FASTQ files: $FASTQ_FILES"

# Build breseq command with polymorphism mode
BRESEQ_CMD="${BRESEQ_PATH}"
if [[ "$BRESEQ_POLY_MODE" == "true" ]]; then
    BRESEQ_CMD="${BRESEQ_CMD} -p"
fi

echo "Running breseq..."
${BRESEQ_CMD} \
    --polymorphism-reject-indel-homopolymer-length ${BRESEQ_POLY_INDEL_HOMOPOLY} \
    --polymorphism-reject-surrounding-homopolymer-length ${BRESEQ_POLY_SURROUND_HOMOPOLY} \
    --polymorphism-score-cutoff ${BRESEQ_POLY_SCORE_CUTOFF} \
    --polymorphism-minimum-variant-coverage-each-strand ${BRESEQ_POLY_MIN_COV_STRAND} \
    -o "${output_dir}/${sample_name}" \
    -r "${REFERENCE_GBK}" \
    ${FASTQ_FILES} \
    > "${output_dir}/${sample_name}.out" \
    2> "${output_dir}/${sample_name}.err"

echo "Done processing $sample_name"
