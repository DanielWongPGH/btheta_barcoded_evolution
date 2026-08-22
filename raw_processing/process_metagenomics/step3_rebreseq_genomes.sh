#!/usr/bin/bash
# rebreseqs samples, and putting in output_directory 
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
config_file=$1
population=$2
sequence_type=$3 #e.g. stragglers, non_stragglers, non_clone_non_stragglers
RAW_PROCESSING_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
source "${RAW_PROCESSING_DIR}/parse_config.sh" "$config_file" metagenomics
breseq_dir=${4:-$BRESEQ_DIR}
output_dir=${5:-$REBRESEQ_DIR}
if [[ $breseq_dir == */ ]]; then
    breseq_dir=${breseq_dir%/}
fi

if [[ $output_dir == */ ]]; then
    output_dir=${output_dir%/}
fi

export gd_files="" 

echo Population is $population
echo Sequence type is $sequence_type
echo Output dir is $output_dir
mkdir -p $output_dir #no err if already exists

SLURM_DIR=$SLURM_LOG_DIR
echo "Slurm log directory: $SLURM_DIR"
mkdir -p $SLURM_DIR

load_breseq_modules

for pop in $(echo $population | tr "," "\n"); do
    for i in $("${PYTHON2_BIN}" "${SCRIPT_DIR}/population_parameters.py" samples ${pop} ${sequence_type}); do
        arrI=(${i//;/ })
        sample_path=${arrI[0]}
        sample_name=${arrI[1]}

        echo sample path/name is $sample_path/$sample_name

        if [ ! -e "${breseq_dir}/${sample_name}/output/evidence/evidence.gd" ]
        then
        echo "Should not get here!"
        rm ${breseq_dir}/${sample_name}/output/JC_output.gd
        continue 
        fi

        cat ${breseq_dir}/${sample_name}/output/evidence/evidence.gd | grep 'JC\|#' > ${breseq_dir}/${sample_name}/output/JC_output.gd 

        export gd_files="${gd_files} ${breseq_dir}/${sample_name}/output/JC_output.gd"
        
    done

    merged_gd_file=${breseq_dir}/${pop}_merged_breseq_output.gd
    echo ${gd_files}
    "${GDTOOLS_PATH}" UNION -o "$merged_gd_file" -e ${gd_files}


    mkdir -p $output_dir
    if [ -e $merged_gd_file ]
    then
    for i in $("${PYTHON2_BIN}" "${SCRIPT_DIR}/population_parameters.py" samples ${pop} ${sequence_type}); do
        arrI=(${i//;/ })
        sample_path=${arrI[0]}
        sample_name=${arrI[1]}
        params=$sample_name

        if [[ $sample_path == */ ]]; then
            sample_path=${sample_path%/}_trimmed
        fi

        echo sample path is $sample_path
        echo sample name is $sample_name


        if [ ! -e "${output_dir}/${sample_name}/output/evidence/evidence.gd" ]
        then
        submit_metagenomics_job --partition="${SLURM_PARTITION}" --job-name=rebreseq-${sample_name} --output=$SLURM_DIR/rebreseq-%j.log --error=$SLURM_DIR/rebreseq-%j.log \
            "${SCRIPT_DIR}/step3_rebreseq_genomes.sbatch" "$config_file" "$sample_path" "$sample_name" "$output_dir" "$merged_gd_file"
        else
            echo ${sample_name} already rebreseqed
        fi
    done
    fi
done
