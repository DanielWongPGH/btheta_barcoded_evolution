#!/usr/bin/bash
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
config_file=$1
export population=$2
export sequence_type=$3
RAW_PROCESSING_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
source "${RAW_PROCESSING_DIR}/parse_config.sh" "$config_file" metagenomics
rebreseq_dir=${4:-$REBRESEQ_DIR}
load_breseq_modules
# output_base_dir=$3
bam_files=""
sample_names=""

export reference_file="$REFERENCE_FASTA"
## if last character in output_base_dir is /, remove it
# if [[ $output_base_dir == */ ]]; then
#     output_base_dir=${output_dir%/}
# fi
# export output_base_dir

# mkdir -p ${output_base_dir}/timecourse_files
## if last character in rebreseq_dir is /, remove it
if [[ $rebreseq_dir == */ ]]; then
    rebreseq_dir=${rebreseq_dir%/}
fi
export rebreseq_dir

SLURM_DIR=$SLURM_LOG_DIR
echo "Slurm log directory: $SLURM_DIR"
mkdir -p $SLURM_DIR

for pop in $(echo $population | tr "," "\n"); do
    echo "Creating timecourse for population: $pop"
    for i in $("${PYTHON2_BIN}" "${SCRIPT_DIR}/population_parameters.py" samples ${pop} ${sequence_type}); do
        arrI=(${i//;/ })
        export sample_path=${arrI[0]}
        export sample_name=${arrI[1]}
        export params=$sample_name

        sample_names="${sample_names},${sample_name}"
        bam_files="${bam_files},${rebreseq_dir}/${sample_name}/data/reference.bam"
    done
    #echo ${bam_files}
    # remove leading comma
    sample_names=${sample_names#*,}
    bam_files=${bam_files#*,}

    echo $sample_names


    for ((i=TIMECOURSE_CHUNK_SIZE; i<=GENOME_LENGTH; i+=TIMECOURSE_CHUNK_SIZE))
    do
        start_position=$((i - TIMECOURSE_CHUNK_SIZE))
        end_position=$i
        echo $start_position $end_position
        submit_metagenomics_job --partition="${SLURM_PARTITION}" --job-name=timecourse_${pop}_${start_position} --output=$SLURM_DIR/timecourse_${pop}_${start_position}.log \
            "${SCRIPT_DIR}/step4_create_timecourse.sbatch" "$config_file" ${start_position} ${end_position} ${pop} ${sample_names} ${bam_files}
    done
done
