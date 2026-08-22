#!/usr/bin/bash
#SBATCH -p normal
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1

# step3_srun_noise_estimate.sh
# Worker script: runs estimate_noise.py for a single library
#
# Usage: Called by step3_sbatch_noise_estimation.sh
#   $1 - project_config.yaml
#   $2 - fastq_dir
#   $3 - lib_id_fastq_map
#   $4 - row number

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

config_file=$1
fastq_dir=$2
lib_id_fastq_map=$3
row=$4

# Load configuration
RAW_PROCESSING_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
source "${RAW_PROCESSING_DIR}/parse_config.sh" "$config_file" barseq

# Get the fastq stem from the specified row
fastq_stem=$(sed -n "${row}p" "$lib_id_fastq_map")

# Set SLURM output paths (note: these are set at job submission time,
# so this is mainly for documentation)
# #SBATCH -o ${SLURM_LOG_DIR}/%x.out
# #SBATCH -e ${SLURM_LOG_DIR}/%x.err

# Load required modules
load_configured_modules "$MODULE_PYTHON" "$MODULE_VIZ" "$MODULE_MATPLOTLIB"

echo "Running noise estimation for: $fastq_stem"
echo "FASTQ dir: $fastq_dir"
echo "Bartender dir: $BARTENDER_NOUMI_DIR"
echo "Output dir: $PSEUDO_DEREP_DIR"

echo "$PYTHON_BIN" estimate_noise.py \
    --fastq "$fastq_stem" \
    --fastq_dir "$fastq_dir" \
    --bartender_dir "$BARTENDER_NOUMI_DIR" \
    -o "$PSEUDO_DEREP_DIR"

run_barseq_command "$PYTHON_BIN" "${SCRIPT_DIR}/estimate_noise.py" \
    --fastq "$fastq_stem" \
    --fastq_dir "$fastq_dir" \
    --bartender_dir "$BARTENDER_NOUMI_DIR" \
    -o "$PSEUDO_DEREP_DIR"
