#!/usr/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1

# call_structural_variants_with_delly.sbatch
# Worker script: runs delly for structural variant calling on a single sample
#
# Usage: Called by call_structural_variants_with_delly.sh
#   $1 - project_config.yaml
#   $2 - sample_name
#   $3 - sample_path_to_bam

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

config_file=$1
sample=$2
sample_path_to_bam=$3

# Load configuration
RAW_PROCESSING_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
source "${RAW_PROCESSING_DIR}/parse_config.sh" "$config_file" isolates

# Load required modules
ml ${MODULE_BIOLOGY} ${MODULE_SAMTOOLS}

echo "Processing sample: $sample"
echo "BAM file: $sample_path_to_bam"
echo "Output: ${DELLY_DIR}/${sample}.bcf"

# Get reference fasta from breseq output (same directory as BAM)
ref_genome="${BRESEQ_DIR}/${sample}/data/reference.fasta"

if [[ ! -e "$ref_genome" ]]; then
    echo "ERROR: Reference fasta not found at $ref_genome"
    exit 1
fi

echo "Reference: $ref_genome"

# Run delly
echo "Running delly..."
singularity exec "${DELLY_SIMG}" delly call \
    -q ${DELLY_MIN_MAPQ} \
    -g "$ref_genome" \
    -o "${DELLY_DIR}/${sample}.bcf" \
    "$sample_path_to_bam"

echo "Completed delly call for $sample"
