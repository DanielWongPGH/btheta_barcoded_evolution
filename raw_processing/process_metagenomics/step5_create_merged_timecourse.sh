#!/usr/bin/bash
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
config_file=$1
export population=$2
export sequence_type=$3
RAW_PROCESSING_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
source "${RAW_PROCESSING_DIR}/parse_config.sh" "$config_file" metagenomics
load_breseq_modules

export reference_file="$REFERENCE_FASTA"

mkdir -p "$TIMECOURSE_DIR"

for pop in $(echo $population | tr "," "\n"); do
    gd_files=""
    sample_names=""
    for i in $("${PYTHON2_BIN}" "${SCRIPT_DIR}/population_parameters.py" samples ${pop} ${sequence_type}); do
        arrI=(${i//;/ })
        export sample_path=${arrI[0]}
        export sample_name=${arrI[1]}
        export params=$sample_name

        export gd_files="${gd_files} ${REBRESEQ_DIR}/${sample_name}/output/evidence/evidence.gd"
    done
    # write everything to a merged file
    # Version with junctions/indels
    echo ${gd_files}
    cat "${TIMECOURSE_DIR}/${pop}_"*timecourse.txt | "${PYTHON2_BIN}" "${SCRIPT_DIR}/step5_create_breseq_timecourse.py" "${reference_file}" ${pop} ${gd_files} | bzip2 -c > "${TIMECOURSE_DIR}/${pop}_merged_timecourse.bz2"
    # Version without junctions and indels
    # cat ${base_directory}/timecourse_files/${population}_*timecourse.txt | bzip2 -c > ${population}_merged_timecourse.bz2
done
