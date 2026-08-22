#!/usr/bin/bash
#SBATCH --time=12:00:00
##SBATCH --mem 32G
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
config_file=$1
input_dir=$2
sample_name=$3
output_dir=$4
merged_gd_file=$5
RAW_PROCESSING_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
source "${RAW_PROCESSING_DIR}/parse_config.sh" "$config_file" metagenomics

# ${output_dir}/${sample_name}/output/JC_output.gd"

export FASTQ_FILES=""
for f in $(ls ${input_dir}/${sample_name}*.fastq.gz)
do
   export FASTQ_FILES="${FASTQ_FILES} ${f}"
done

"${BRESEQ_PATH}" -p --user-evidence-gd "$merged_gd_file" -o "${output_dir}/${sample_name}" -r "${REFERENCE_GBK}" ${FASTQ_FILES} > "${output_dir}/${sample_name}.out" 2> "${output_dir}/${sample_name}.err"
