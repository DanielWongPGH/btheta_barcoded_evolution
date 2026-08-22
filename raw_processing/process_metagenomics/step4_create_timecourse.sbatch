#!/usr/bin/bash
#SBATCH --time=2:00:00
##SBATCH --mem 10G
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
config_file=$1
start_position=$2
end_position=$3
population=$4
comma_separated_sample_names=$5
comma_separated_bam_files=$6
RAW_PROCESSING_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
source "${RAW_PROCESSING_DIR}/parse_config.sh" "$config_file" metagenomics
sample_names=$(echo $comma_separated_sample_names | tr ',' ' ')
bam_files=$(echo $comma_separated_bam_files | tr ',' ' ')




pileup_file="${TIMECOURSE_DIR}/${population}_${start_position}.pileup"
timecourse_file="${TIMECOURSE_DIR}/${population}_${start_position}_timecourse.txt"
samtools mpileup -q"${MPILEUP_QUALITY}" -f "${REFERENCE_FASTA}" -r "${REFERENCE_CHROMOSOME}:${start_position}-${end_position}" ${bam_files} > "$pileup_file"

"${PYTHON2_BIN}" "${SCRIPT_DIR}/step4_create_timecourse.py" ${start_position} ${end_position} ${population} ${sample_names} < "$pileup_file" > "$timecourse_file"

# The pileup is retained as an intermediate so a failed parser can be rerun.
