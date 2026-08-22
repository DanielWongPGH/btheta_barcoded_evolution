import os
import sys

fastq_dir = sys.argv[-2].rstrip('/')
output_name = sys.argv[-1]
all_files = os.listdir(fastq_dir)
samples = []

analyzed_samples = []
with open(f'clone_barcode.tsv', 'r') as f:
    header = next(f)
    for line in f:
        line = line.rstrip('\n').split('\t')
        mouse_day_clone = f'm{line[1]}_day{line[2]}_clone{line[3]}'
        analyzed_samples.append(mouse_day_clone)

for file in all_files:
    if 'fastq' not in file or file[0] == '.':
        continue
    sample = file.rstrip('.gz').rstrip('.fastq').rstrip('_R1_001').rstrip('_R2_001')

    if sample not in samples and sum([x in sample for x in analyzed_samples]) > 0:
        samples.append(sample)

f = open(f'{output_name.rstrip(".txt")}.txt', 'w')
f.write(f'{fastq_dir}' + '\n')
for sample in samples:
    f.write(f'{sample}' + '\n')
f.close()
