#!/usr/bin/bash
#SBATCH --ntasks=1
#SBATCH --no-requeue

# step0_trim_adapters.sbatch
# Worker script: trims adapters from a single sample using Trimmomatic
#
# Usage: Called by step0_trim_adapters.sh via sbatch array
#   $1 - project_config.yaml
#   $2 - samples_file (first line is directory, subsequent lines are sample names)
#   $3 - trimmed_dir

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

config_file=$1
samples_file=$2
trimmed_dir=$3

# Load configuration
RAW_PROCESSING_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
source "${RAW_PROCESSING_DIR}/parse_config.sh" "$config_file" metagenomics

# Get untrimmed directory (first line) and sample name (line matching array task ID)
untrimmed_dir=$(head -1 "$samples_file")
sample=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$samples_file")

echo "Processing sample: $sample"
echo "Input dir: $untrimmed_dir"
echo "Output dir: $trimmed_dir"

# Accept compressed or uncompressed inputs without modifying the raw files.
resolve_fastq() {
    local gz_path="$1"
    if [[ -f "$gz_path" ]]; then
        printf '%s\n' "$gz_path"
    elif [[ -f "${gz_path%.gz}" ]]; then
        printf '%s\n' "${gz_path%.gz}"
    else
        echo "ERROR: FASTQ not found: $gz_path or ${gz_path%.gz}" >&2
        return 1
    fi
}

in1="$(resolve_fastq "${untrimmed_dir}/${sample}_R1_001.fastq.gz")"
in2="$(resolve_fastq "${untrimmed_dir}/${sample}_R2_001.fastq.gz")"

out1=${trimmed_dir}/${sample}_R1_001.fastq.gz
out2=${trimmed_dir}/${sample}_R2_001.fastq.gz
out3=${trimmed_dir}/${sample}_R1_001.unpaired.fastq.gz
out4=${trimmed_dir}/${sample}_R2_001.unpaired.fastq.gz

echo "Running Trimmomatic..."
java -jar "${TRIMMOMATIC_JAR}" PE \
    "$in1" "$in2" \
    "$out1" "$out3" "$out2" "$out4" \
    ILLUMINACLIP:${TRIMMOMATIC_ADAPTERS}:${TRIM_ILLUMINA_CLIP} \
    LEADING:${TRIM_LEADING} \
    TRAILING:${TRIM_TRAILING} \
    MINLEN:${TRIM_MINLEN}

touch "${trimmed_dir}/${sample}.finished"

echo "Done trimming $sample"
