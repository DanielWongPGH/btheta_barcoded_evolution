import os
import sys

mgx_fastq_dir = sys.argv[-2].rstrip('/')
output_name = sys.argv[-1]
all_files = os.listdir(mgx_fastq_dir)
samples = []

for file in all_files:
    if 'fastq' not in file or file[0] == '.':
        continue
    sample = file.replace('.gz', '').replace('.fastq', '').replace('_R1_001', '').replace('_R2_001', '')

    if sample not in samples:
        samples.append(sample)

f = open(f'{output_name}', 'w') 
f.write(f'{mgx_fastq_dir}' + '\n')
for sample in samples:
    f.write(f'{sample}' + '\n')
f.close()
