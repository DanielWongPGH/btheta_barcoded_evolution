#!/usr/bin/bash
#SBATCH --ntasks=1
#SBATCH --no-requeue

# step0_trim_adapters.sbatch
# Worker script: trims adapters from a single sample using Trimmomatic
#
# Usage: Called by step0_trim_adapters.sh via sbatch array
#   $1 - project_config.yaml
#   $2 - samples_file (first line is directory, subsequent lines are sample names)

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

config_file=$1
samples_file=$2

# Load configuration
RAW_PROCESSING_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
source "${RAW_PROCESSING_DIR}/parse_config.sh" "$config_file" isolates
load_isolate_modules "$MODULE_PYTHON"

# Get untrimmed directory (first line) and sample name (line matching array task ID)
untrimmed_dir=$(head -1 "$samples_file")
sample=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$samples_file")

echo "Processing sample: $sample"
echo "Input dir: $untrimmed_dir"
echo "Output dir: $TRIMMED_DIR"

# Define input/output files
in1=""
in2=""
for suffix in .fastq.gz .fastq; do
    [[ -f "${untrimmed_dir}/${sample}_R1_001${suffix}" ]] && in1="${untrimmed_dir}/${sample}_R1_001${suffix}"
    [[ -f "${untrimmed_dir}/${sample}_R2_001${suffix}" ]] && in2="${untrimmed_dir}/${sample}_R2_001${suffix}"
done
if [[ -z "$in1" || -z "$in2" ]]; then
    echo "Error: paired FASTQs not found for $sample" >&2
    exit 1
fi
out1=${TRIMMED_DIR}/${sample}_R1_001.fastq.gz
out2=${TRIMMED_DIR}/${sample}_R2_001.fastq.gz
out3=${TRIMMED_DIR}/${sample}_R1_001.unpaired.fastq.gz
out4=${TRIMMED_DIR}/${sample}_R2_001.unpaired.fastq.gz

echo "Running Trimmomatic..."
"$JAVA_BIN" -jar "${TRIMMOMATIC_JAR}" PE \
    "$in1" "$in2" \
    "$out1" "$out3" "$out2" "$out4" \
    ILLUMINACLIP:${TRIMMOMATIC_ADAPTERS}:${TRIM_ILLUMINA_CLIP} \
    LEADING:${TRIM_LEADING} \
    TRAILING:${TRIM_TRAILING} \
    MINLEN:${TRIM_MINLEN}

touch "${TRIMMED_DIR}/${sample}.finished"

echo "Done trimming $sample"
