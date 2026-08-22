#!/bin/bash
# Load one raw-processing section from the unified project configuration.
# Usage: source parse_config.sh <project_config.yaml> <barseq|isolates|metagenomics>

CONFIG_FILE="${1:-}"
PIPELINE="${2:-}"

if [[ -z "$CONFIG_FILE" || -z "$PIPELINE" ]]; then
    echo "Usage: source parse_config.sh <project_config.yaml> <barseq|isolates|metagenomics>" >&2
    return 1 2>/dev/null || exit 1
fi

if [[ ! -f "$CONFIG_FILE" ]]; then
    echo "Error: config file not found: $CONFIG_FILE" >&2
    return 1 2>/dev/null || exit 1
fi

if ! command -v yq >/dev/null 2>&1; then
    echo "Error: mikefarah/yq version 4 is required." >&2
    return 1 2>/dev/null || exit 1
fi

case "$PIPELINE" in
    barseq|isolates|metagenomics) ;;
    *)
        echo "Error: unknown raw-processing pipeline: $PIPELINE" >&2
        return 1 2>/dev/null || exit 1
        ;;
esac

config_value() {
    command yq -r ".raw_processing.${PIPELINE}${1}" "$CONFIG_FILE"
}

if [[ "$(config_value '')" == "null" ]]; then
    echo "Error: missing raw_processing.${PIPELINE} in $CONFIG_FILE" >&2
    return 1 2>/dev/null || exit 1
fi

if [[ "$PIPELINE" == "barseq" ]]; then
    export FASTQ_DIR="$(config_value '.fastq_dir')"
    export OUTPUT_DIR="$(config_value '.output_dir')"
    export BARTENDER_PATH="$(config_value '.bartender_path')"
    export EXECUTION_SCHEDULER="$(config_value '.execution.scheduler')"
    export LOAD_MODULES="$(config_value '.execution.load_modules')"
    export PYTHON_BIN="$(config_value '.execution.python')"
    export SLURM_LOG_DIR="${OUTPUT_DIR}/slurm"
    export MERGED_FASTQ_DIR="${OUTPUT_DIR}/merged_fastq"
    export BARTENDER_DIR="${OUTPUT_DIR}/bartender"
    export BARTENDER_NOUMI_DIR="${OUTPUT_DIR}/bartender_noUMI"
    export PSEUDO_DEREP_DIR="${OUTPUT_DIR}/pseudo_derep"
    export SLURM_PARTITION="$(config_value '.slurm.partition')"
    export SLURM_MEMORY="$(config_value '.slurm.memory')"
    export SLURM_TIME="$(config_value '.slurm.time')"
    export MODULE_PYTHON="$(config_value '.modules.python')"
    export MODULE_VIZ="$(config_value '.modules.viz')"
    export MODULE_MATPLOTLIB="$(config_value '.modules.matplotlib')"
    export MODULE_GCC="$(config_value '.modules.gcc')"
    export BARTENDER_PATTERN="$(config_value '.bartender.pattern')"
    export BARTENDER_EXTRACTOR="$(config_value '.bartender.extractor_executable')"
    export BARTENDER_CLUSTER="$(config_value '.bartender.cluster_executable')"
    export BARTENDER_QUALITY="$(config_value '.bartender.quality_char')"
    export BARTENDER_SEED="$(config_value '.bartender.seed_length')"
    export BARTENDER_DISTANCE="$(config_value '.bartender.distance')"
    export BARTENDER_UMI_RANGE="$(config_value '.bartender.umi_range')"
    export PRIMER_R1="$(config_value '.primers.R1')"
    export PRIMER_R2="$(config_value '.primers.R2')"

    create_output_dirs() {
        mkdir -p "$SLURM_LOG_DIR" "$MERGED_FASTQ_DIR" "$BARTENDER_DIR" \
            "$BARTENDER_NOUMI_DIR" "$PSEUDO_DEREP_DIR"
    }

    load_configured_modules() {
        [[ "$LOAD_MODULES" == "true" ]] || return 0
        ml "$@"
    }

    run_barseq_command() {
        if [[ "$EXECUTION_SCHEDULER" == "slurm" ]]; then
            srun "$@"
        else
            "$@"
        fi
    }
fi

if [[ "$PIPELINE" == "isolates" ]]; then
    export FASTQ_DIR="$(config_value '.fastq_dir')"
    export OUTPUT_DIR="$(config_value '.output_dir')"
    export TRIMMED_DIR="${OUTPUT_DIR}/trimmed"
    export BRESEQ_DIR="${OUTPUT_DIR}/breseq"
    export DELLY_DIR="${OUTPUT_DIR}/delly"
    export SLURM_LOG_DIR="${OUTPUT_DIR}/slurm"
    export TRIMMOMATIC_JAR="$(config_value '.trimmomatic.jar')"
    export TRIMMOMATIC_ADAPTERS="$(config_value '.trimmomatic.adapters')"
    export BRESEQ_PATH="$(config_value '.breseq.path')"
    export DELLY_SIMG="$(config_value '.delly.singularity_image')"
    export REFERENCE_GBK="$(config_value '.reference.genbank')"
    export SLURM_PARTITION="$(config_value '.slurm.partition')"
    export SLURM_TRIM_TIME="$(config_value '.slurm.trim.time')"
    export SLURM_TRIM_MEM="$(config_value '.slurm.trim.mem_per_cpu')"
    export SLURM_TRIM_CPUS="$(config_value '.slurm.trim.cpus')"
    export SLURM_BRESEQ_TIME="$(config_value '.slurm.breseq.time')"
    export SLURM_BRESEQ_MEM="$(config_value '.slurm.breseq.mem')"
    export SLURM_DELLY_TIME="$(config_value '.slurm.delly.time')"
    export SLURM_DELLY_MEM="$(config_value '.slurm.delly.mem')"
    export MODULE_R="$(config_value '.modules.R')"
    export MODULE_BIOLOGY="$(config_value '.modules.biology')"
    export MODULE_SAMTOOLS="$(config_value '.modules.samtools')"
    export MODULE_SYSTEM="$(config_value '.modules.system')"
    export MODULE_NCURSES="$(config_value '.modules.ncurses')"
    export MODULE_BOWTIE2="$(config_value '.modules.bowtie2')"
    export MODULE_NUMPY="$(config_value '.modules.numpy')"
    export TRIM_ILLUMINA_CLIP="$(config_value '.trimmomatic_params.illumina_clip')"
    export TRIM_LEADING="$(config_value '.trimmomatic_params.leading')"
    export TRIM_TRAILING="$(config_value '.trimmomatic_params.trailing')"
    export TRIM_MINLEN="$(config_value '.trimmomatic_params.minlen')"
    export BRESEQ_POLY_INDEL_HOMOPOLY="$(config_value '.breseq_params.polymorphism_reject_indel_homopolymer_length')"
    export BRESEQ_POLY_SURROUND_HOMOPOLY="$(config_value '.breseq_params.polymorphism_reject_surrounding_homopolymer_length')"
    export BRESEQ_POLY_SCORE_CUTOFF="$(config_value '.breseq_params.polymorphism_score_cutoff')"
    export BRESEQ_POLY_MIN_COV_STRAND="$(config_value '.breseq_params.polymorphism_minimum_variant_coverage_each_strand')"
    export DELLY_MIN_MAPQ="$(config_value '.delly_params.min_mapping_quality')"

    create_output_dirs() {
        mkdir -p "$TRIMMED_DIR" "$BRESEQ_DIR" "$DELLY_DIR" "$SLURM_LOG_DIR"
    }
fi

if [[ "$PIPELINE" == "metagenomics" ]]; then
    export SCRATCH_DIR="$(config_value '.scratch_dir')"
    export EXECUTION_SCHEDULER="$(config_value '.execution.scheduler')"
    export LOAD_MODULES="$(config_value '.execution.load_modules')"
    export PYTHON3_BIN="$(config_value '.execution.python')"
    export BRESEQ_DIR="${SCRATCH_DIR}/$(config_value '.output_dirs.breseq')"
    export REBRESEQ_DIR="${SCRATCH_DIR}/$(config_value '.output_dirs.rebreseq')"
    export TIMECOURSE_DIR="${SCRATCH_DIR}/$(config_value '.output_dirs.timecourse')"
    export SLURM_LOG_DIR="${SCRATCH_DIR}/$(config_value '.output_dirs.slurm_logs')"
    export PYTHON2_BIN="$(config_value '.python2.executable')"
    export TRIMMOMATIC_JAR="$(config_value '.trimmomatic.jar')"
    export TRIMMOMATIC_ADAPTERS="$(config_value '.trimmomatic.adapters')"
    export BRESEQ_PATH="$(config_value '.breseq.path')"
    export GDTOOLS_PATH="$(config_value '.breseq.gdtools')"
    export REFERENCE_GBK="$(config_value '.reference.genbank')"
    export REFERENCE_FASTA="$(config_value '.reference.fasta')"
    export REFERENCE_CHROMOSOME="$(config_value '.reference.chromosome')"
    export SLURM_PARTITION="$(config_value '.slurm.partition')"
    export SLURM_TRIM_TIME="$(config_value '.slurm.trim.time')"
    export SLURM_TRIM_MEM="$(config_value '.slurm.trim.mem_per_cpu')"
    export SLURM_TRIM_CPUS="$(config_value '.slurm.trim.cpus')"
    export SLURM_BRESEQ_TIME="$(config_value '.slurm.breseq.time')"
    export SLURM_BRESEQ_MEM="$(config_value '.slurm.breseq.mem')"
    export SLURM_REBRESEQ_TIME="$(config_value '.slurm.rebreseq.time')"
    export SLURM_REBRESEQ_MEM="$(config_value '.slurm.rebreseq.mem')"
    export SLURM_TIMECOURSE_TIME="$(config_value '.slurm.timecourse.time')"
    export SLURM_TIMECOURSE_MEM="$(config_value '.slurm.timecourse.mem')"
    export MODULE_R="$(config_value '.modules.R')"
    export MODULE_BIOLOGY="$(config_value '.modules.biology')"
    export MODULE_SAMTOOLS="$(config_value '.modules.samtools')"
    export MODULE_SYSTEM="$(config_value '.modules.system')"
    export MODULE_NCURSES="$(config_value '.modules.ncurses')"
    export MODULE_BOWTIE2="$(config_value '.modules.bowtie2')"
    export MODULE_NUMPY="$(config_value '.modules.numpy')"
    export MODULE_SCIPY="$(config_value '.modules.scipy')"
    export MODULE_GCC="$(config_value '.modules.gcc')"
    export TRIM_ILLUMINA_CLIP="$(config_value '.trimmomatic_params.illumina_clip')"
    export TRIM_LEADING="$(config_value '.trimmomatic_params.leading')"
    export TRIM_TRAILING="$(config_value '.trimmomatic_params.trailing')"
    export TRIM_MINLEN="$(config_value '.trimmomatic_params.minlen')"
    export BRESEQ_POLY_MODE="$(config_value '.breseq_params.polymorphism_mode')"
    export BRESEQ_POLY_INDEL_HOMOPOLY="$(config_value '.breseq_params.polymorphism_reject_indel_homopolymer_length')"
    export BRESEQ_POLY_SURROUND_HOMOPOLY="$(config_value '.breseq_params.polymorphism_reject_surrounding_homopolymer_length')"
    export BRESEQ_POLY_SCORE_CUTOFF="$(config_value '.breseq_params.polymorphism_score_cutoff')"
    export BRESEQ_POLY_MIN_COV_STRAND="$(config_value '.breseq_params.polymorphism_minimum_variant_coverage_each_strand')"
    export MPILEUP_QUALITY="$(config_value '.timecourse_params.mpileup_quality')"
    export TIMECOURSE_CHUNK_SIZE="$(config_value '.timecourse_params.chunk_size')"
    export GENOME_LENGTH="$(config_value '.timecourse_params.genome_length')"

    create_output_dirs() {
        mkdir -p "$BRESEQ_DIR" "$REBRESEQ_DIR" "$TIMECOURSE_DIR" "$SLURM_LOG_DIR"
    }

    load_common_modules() {
        [[ "$LOAD_MODULES" == "true" ]] || return 0
        ml "$MODULE_R" "$MODULE_BIOLOGY" "$MODULE_SAMTOOLS" "$MODULE_SYSTEM" "$MODULE_NCURSES"
    }

    load_breseq_modules() {
        [[ "$LOAD_MODULES" == "true" ]] || return 0
        ml "$MODULE_GCC"
        ml "$MODULE_R" "$MODULE_BIOLOGY" "$MODULE_SAMTOOLS" "$MODULE_SYSTEM" \
            "$MODULE_NCURSES" "$MODULE_BOWTIE2"
        if [[ -n "$MODULE_NUMPY" ]]; then ml "$MODULE_NUMPY"; fi
    }

    submit_metagenomics_job() {
        if [[ "$EXECUTION_SCHEDULER" != "slurm" ]]; then
            echo "Error: metagenomics execution.scheduler must currently be slurm" >&2
            return 1
        fi
        sbatch "$@"
    }
fi
