import sys
import os
import gzip

directory = sys.argv[1]
file = sys.argv[2]
output_dir = sys.argv[3]

if directory[-1] == '/':
    directory = directory[:-1]
if output_dir[-1] == '/':
    output_dir = output_dir[:-1]

if "_R1_001.fastq" in file:
    file_stem = file.split("_R1_001.fastq")[0]
elif "_R2_001.fastq" in file:
    file_stem = file.split("_R2_001.fastq")[0]
else:
    file_stem = file

R1_file = file_stem + "_R1_001.fastq"
R2_file = file_stem + "_R2_001.fastq"

def open_fastq(filename):
    path = f'{directory}/{filename}'
    if os.path.exists(path):
        return open(path, 'r')
    if os.path.exists(path + '.gz'):
        return gzip.open(path + '.gz', 'rt')
    raise FileNotFoundError(path)

sR1 = 'TGCAGTGTCG'
sR2 = 'AAATGCTGTT'

merged_file = 'merged_' + file_stem + '.fastq'
sys.stdout.write(f'Merging {R1_file} and {R2_file} in {output_dir}/{merged_file}...\n')

if os.path.exists(f'{output_dir}/{merged_file}'):
    exit()

fout = open(f'{output_dir}/{merged_file}', 'w')
f1 = open_fastq(R1_file)
f2 = open_fastq(R2_file)

def hamming(str1, str2):
    return sum([str1[i] != str2[i] for i in range(len(str1))])

tot_lines = 0
kept_lines = 0
for z, (lineR1, lineR2) in enumerate(zip(f1,f2)):
    tot_lines += 1
    if z % 4 == 0:
        flag = False
        paired_readlines = lineR1
    elif z % 4 == 1:
        sR2_sequence_str = lineR2[8:8+len(sR2)]
        sR1_sequence_str = lineR1[8:8+len(sR1)]
        if hamming(sR2_sequence_str, sR2) > 2 or hamming(sR1_sequence_str, sR1) > 2:
            flag = True
        out_line = lineR2[:8] + lineR1
        paired_readlines += out_line
    elif z % 4 == 2:
        paired_readlines += lineR1
    elif z % 4 == 3:
        paired_readlines += lineR2[:8] + lineR1
        if flag:
            continue
        else:
            fout.write(paired_readlines)
            kept_lines += 4

sys.stdout.write(f'Tot reads: {int(tot_lines / 4)}, Kept reads: {int(kept_lines / 4)}\n')
f1.close()
f2.close()
fout.close()
