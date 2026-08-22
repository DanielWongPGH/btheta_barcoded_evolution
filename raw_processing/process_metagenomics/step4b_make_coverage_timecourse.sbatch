#!/usr/bin/bash
#SBATCH --time=2:00:00
##SBATCH --mem 32G
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1


SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
config_file=$1
genome_bed=$2
population=$3
comma_separated_sample_names=$4
comma_separated_bam_files=$5
comma_separated_summary_jsons=$6
RAW_PROCESSING_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
source "${RAW_PROCESSING_DIR}/parse_config.sh" "$config_file" metagenomics
sample_names=$(echo $comma_separated_sample_names | tr ',' ' ')
bam_files=$(echo $comma_separated_bam_files | tr ',' ' ')
summary_jsons=$(echo $comma_separated_summary_jsons | tr ',' ' ')

load_common_modules

escaped_sample_names=$(echo "$comma_separated_sample_names" | sed 's/,/\\t/g')
header="#CHROM\\tPOS\\t$escaped_sample_names"
echo "$header"

echo "samtools depth -b $genome_bed -Q 10 -o ${out_file} ${bam_files}"
# samtools samtools depth -b $genome_bed -Q 20 -o ${out_file} ${bam_files}

out_file="${TIMECOURSE_DIR}/${population}_coverage_timecourse.txt"
{
  echo -e "$header"

  # Calculate total reads from breseq summary.json files
  total_reads_row="TOTAL_READS\t0"
  for json in ${summary_jsons}; do
    total=$("${PYTHON3_BIN}" -c "import json; d=json.load(open('${json}')); print(d['reads']['total_reads'])")
    total_reads_row="${total_reads_row}\t${total}"
  done
  echo -e "$total_reads_row"

  # Calculate aligned reads from breseq summary.json files
  aligned_reads_row="ALIGNED_READS\t0"
  for json in ${summary_jsons}; do
    aligned=$("${PYTHON3_BIN}" -c "import json; d=json.load(open('${json}')); print(d['reads']['total_aligned_reads'])")
    aligned_reads_row="${aligned_reads_row}\t${aligned}"
  done
  echo -e "$aligned_reads_row"

  # Calculate genome-wide mean coverage for each BAM file
  genome_wide_coverage="GENOME_WIDE_COV\t0"
  for bam in ${bam_files}; do
    mean_depth=$(samtools coverage -r $(head -n1 ${genome_bed} | cut -f1) ${bam} 2>/dev/null | awk 'NR==2 {print $7}')
    if [ -z "$mean_depth" ]; then
      # Fallback: calculate from all contigs if single contig fails
      mean_depth=$(samtools coverage ${bam} | awk 'NR>1 {sum+=$7*$3; len+=$3} END {if(len>0) print sum/len; else print 0}')
    fi
    genome_wide_coverage="${genome_wide_coverage}\t${mean_depth}"
  done
  echo -e "$genome_wide_coverage"

  samtools depth -b ${genome_bed} -Q 10 ${bam_files}  # add any flags as needed
} > ${out_file}

echo completed
