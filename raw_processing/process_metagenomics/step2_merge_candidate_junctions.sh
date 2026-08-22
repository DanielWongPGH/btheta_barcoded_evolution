#!/usr/bin/bash
# Generates a merged GD file from junction evidence across each population.
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
config_file=$1
export population=$2 # or comma-separated list of populations
export sequence_type=$3
breseq_dir=${4:-}
RAW_PROCESSING_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
source "${RAW_PROCESSING_DIR}/parse_config.sh" "$config_file" metagenomics
breseq_dir=${breseq_dir:-$BRESEQ_DIR}

export gd_files="" 


load_breseq_modules

for pop in $(echo $population | tr "," "\n"); do
    for i in $("${PYTHON2_BIN}" "${SCRIPT_DIR}/population_parameters.py" samples ${pop} ${sequence_type}); do
        arrI=(${i//;/ })
        sample_path=${arrI[0]}
        sample_name=${arrI[1]}

        echo sample path is $sample_path
        echo sample name is $sample_name

        if [ ! -e "${breseq_dir}/${sample_name}/output/evidence/evidence.gd" ]
        then
        echo "Should not get here!"
        rm ${breseq_dir}/${sample_name}/output/JC_output.gd
        continue 
        fi

        cat ${breseq_dir}/${sample_name}/output/evidence/evidence.gd | grep 'JC\|#' > ${breseq_dir}/${sample_name}/output/JC_output.gd 

        export gd_files="${gd_files} ${breseq_dir}/${sample_name}/output/JC_output.gd"
        
    done

    echo ${gd_files}
    "${GDTOOLS_PATH}" UNION -o "${breseq_dir}/${pop}_merged_breseq_output.gd" -e ${gd_files}
done
