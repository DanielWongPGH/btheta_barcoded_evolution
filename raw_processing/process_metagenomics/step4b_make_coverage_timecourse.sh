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
export bam_files=""
export sample_names=""

genome_bed="${SCRIPT_DIR}/genome_positions.bed"


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
    bam_files=""
    sample_names=""
    summary_jsons=""
    for i in $("${PYTHON2_BIN}" "${SCRIPT_DIR}/population_parameters.py" samples ${pop} ${sequence_type}); do
        arrI=(${i//;/ })
        export sample_path=${arrI[0]}
        export sample_name=${arrI[1]}
        export params=$sample_name

        export sample_names="${sample_name},${sample_names}"
        export bam_files="${rebreseq_dir}/${sample_name}/data/reference.bam,${bam_files}"
        export summary_jsons="${rebreseq_dir}/${sample_name}/output/summary.json,${summary_jsons}"
    done
    sample_names=${sample_names%,}
    bam_files=${bam_files%,}
    summary_jsons=${summary_jsons%,}

    submit_metagenomics_job --partition="${SLURM_PARTITION}" --job-name=coverage_${pop} --output=$SLURM_DIR/coverage_${pop}.log \
        "${SCRIPT_DIR}/step4b_make_coverage_timecourse.sbatch" "$config_file" ${genome_bed} ${pop} ${sample_names} ${bam_files} ${summary_jsons}
done
## REDO E1_8 and E1_10
