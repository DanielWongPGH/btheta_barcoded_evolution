
import argparse
import pickle
import sys
import os
import numpy as np
from pathlib import Path

from step4_process_noise_estimates import canonical_noise_id

parser = argparse.ArgumentParser(description='Combine Bartender clusters into the master barcode-count table.')
parser.add_argument('data_dir', type=Path, help='Root containing E1_barseq, E2_barseq, vitro_barseq, and rebarseq')
parser.add_argument('--noise-estimates', type=Path, required=True,
                    help='Noise-fit pickle written by step4_process_noise_estimates.py')
parser.add_argument('--output-dir', type=Path, default=Path('.'))
parser.add_argument('--relabels', type=Path, default=Path(__file__).with_name('library_mislabelings.tsv'))
args = parser.parse_args()
data_dir = args.data_dir
args.output_dir.mkdir(parents=True, exist_ok=True)

with args.noise_estimates.open('rb') as handle:
    noise_estimates = pickle.load(handle)

# some libraries were mislabeled
library_relabels = {}
day_relabels = {}

CEC_DAY = 100 #dummy variable
INT_DAY = 200 #dummy variable

with args.relabels.open() as f:
    header = next(f)
    for line in f:
        fastq_label_str, true_label_str = line.strip('\n').split('\t')
        fastq_label_items = fastq_label_str.strip('"').split(',')
        true_label_str = true_label_str.strip('"').split(',')

        fastq_label = (fastq_label_items[0], int(fastq_label_items[1]), int(fastq_label_items[2]), int(fastq_label_items[3]))
        if true_label_str[0] == 'None':
            true_label = None
        else:
            true_label = (true_label_str[0], int(true_label_str[1]), int(true_label_str[2]), int(true_label_str[3]))
        library_relabels[fastq_label] = true_label


def parse_file(directory, file):
    if directory in 'E1': #files have form E1_dayX_mY_...
        if 'tube2' in file or 'S2for9_28' in file:
            return None
        file_items = file.split('_')
        if 'Plate' in file_items[1]:
            mouse = file_items[1].replace('Plate', 'P')
            day = 0
        else:
            mouse = int(file_items[1].replace('m', ''))
            day = file_items[2].replace('day', '')
            if 'cec' in file:
                day = CEC_DAY
            else:
                day = int(day)
        rep = 1
        return ('E1', mouse, day, rep)
    
    if directory == 'E2':
        if '118-S' in file: 
            return None
        file_items = file.split('_')

        if 'day0' in file: #input library
            mouse = file_items[1]
            day = 0
            if '116-S' in file_items[1]:
                mouse = file_items[1].replace('116-S', 'S')+'delay'
                day = 0
        else:
            mouse = int(file_items[1].replace('m', ''))
            day = file_items[2].replace('day', '')
            if 'cec' in file:
                day = CEC_DAY
            else:
                day = int(day)
        rep = 1
        return ('E2', mouse, day, rep)
    
    if directory == 'vitro':
        file_items = file.split('_')
        passage = int(file_items[1].replace('passage', ''))
        if 'V' in file_items[2]: #input library
            well = file_items[2]
        else:
            well = file_items[2].replace('Plate', 'p') + '-' + file_items[3]
        rep = 1
        return ('EV', well, passage, rep)
    
    if directory == 'rebarseq':
        file_items = file.split('_')
        expt, mouse_str, day_str, source = file_items[:4]
        if 'rep2' in file:
            rep = 3
        else:
            rep = 2
        if np.any( [e in mouse_str for e in ['S', 'V', '116-S']] ): # input library or mouse/well
            mouse = mouse_str.replace('m', '')
            if '116' in mouse_str:
                mouse = mouse.lstrip('116-')+'delay'
        else:
            mouse = int(mouse_str.replace('m', ''))

        if source == 'cec': day = CEC_DAY
        elif source == 'int': day = INT_DAY
        elif source == 'vit': day = 0   # input library
        else:
            day = int(day_str.replace('d', ''))
        return (expt, mouse, day, rep)

            
## directories to bartender output (*.cluster files)
E1_dir = data_dir / 'E1_barseq' / 'bartender_noUMI'
E2_dir = data_dir / 'E2_barseq' / 'bartender_noUMI'
vitro_dir = data_dir / 'vitro_barseq' / 'bartender_noUMI'
rebarseq_dir = data_dir / 'rebarseq' / 'bartender'

### get cluster files for each library
E1_bartender_files = [file for file in os.listdir(E1_dir) if ('cluster.csv' in file)] #and ('pcr' not in file)]
E2_bartender_files = [file for file in os.listdir(E2_dir) if ('cluster.csv' in file)] #and ('pcr' not in file)]
vitro_bartender_files = [file for file in os.listdir(vitro_dir) if 'cluster.csv' in file]
rebarseq_files = [file for file in os.listdir(rebarseq_dir) if ('pcr_cluster.csv' in file)] #and ('pcr' not in file)]


### Make up look up table of barcode to reads in each library
barcode_dict = {}
library_names = []
library_tot_reads = {}


for file in rebarseq_files:
    f = open(f'{rebarseq_dir}/{file}', 'r')
    header = next(f)
    library_tuple = parse_file('rebarseq', file)
    if library_tuple == None:# or library_tuple[2] == INT_DAY:
        continue

    if library_tuple in library_relabels:
        if library_relabels[library_tuple] == None:
            print('RELABEL:', file, library_tuple, library_relabels[library_tuple])
            continue
        library_tuple = library_relabels[library_tuple]
    library_names.append(library_tuple)
    library_tot_reads[library_tuple] = [0,0, 0]
    
    tot_reads = 0
    for line in f:
        line_items = line.strip('\n').split(',')
        barcode, reads = line_items[1], int(line_items[-1])
        library_tot_reads[library_tuple][0] += reads #Nominal depth D_nominal
        library_tot_reads[library_tuple][1] += reads #D_eff = D_nominal

        if barcode not in barcode_dict:
            barcode_dict[barcode] = {}
        barcode_dict[barcode][library_tuple] = reads
    # library_tot_reads[library_tuple][1] = library_tot_reads[0]
    # for barcode in barcode_dict:
    #     if library_tuple in barcode_dict[barcode]:
    #         reads = barcode_dict[barcode][library_tuple]
    #         barcode_dict[barcode][library_tuple] = (reads, reads/tot_reads)


for experiment, directory, files in zip(['E1', 'E2', 'vitro'],
                            [E1_dir, E2_dir, vitro_dir],
                            [E1_bartender_files, E2_bartender_files, vitro_bartender_files]):
    for file in sorted(files):
        f = open(f'{directory}/{file}', 'r')
        header = next(f)
        library_tuple = parse_file(experiment, file)
        if library_tuple == None:
            continue
        else:
            (expt, mouse, day, rep) = library_tuple

        noise_id = canonical_noise_id(file.removesuffix('_cluster.csv'), library_relabels)
        if noise_id is None:
            continue
        try:
            noise_dict = noise_estimates[noise_id]
        except KeyError as error:
            raise KeyError('No step-4 noise estimate for {} ({})'.format(file, noise_id)) from error
        p, Z1, Z2, D1, D2 = [float(noise_dict[k]) for k in ['p', 'Z1', 'Z2', 'D1', 'D2']]
        if Z2 > Z1: # if somehow these were switched during inference, make Z1 the larger one
            Z1, Z2 = Z2, Z1
            D1, D2 = D2, D1
            p = 1-p


        if library_tuple in library_relabels:
            if library_relabels[library_tuple] == None:
                print('RELABEL:', file, library_tuple, library_relabels[library_tuple])
                continue
            library_tuple = library_relabels[library_tuple]
        library_names.append(library_tuple)
        library_tot_reads[library_tuple] = [0, 0, 0]

        
        excluded_barcode_count = 0
        excluded_barcodes = []
        tot_count = 0
        for line in f:
            line_items = line.strip('\n').split(',')
            barcode, reads = line_items[1], int(line_items[-1])
            library_tot_reads[library_tuple][0] += reads #nominal depth D_nominal

            if barcode not in barcode_dict:
                barcode_dict[barcode] = {}

            barcode_dict[barcode][library_tuple] = reads
            tot_count += reads

        if Z2 * D2 < 3*Z1 * D1:
            D_eff = D1 * (1-np.exp(-(Z1+Z2*D2/D1))) 
            library_tot_reads[library_tuple][1] = D_eff
        else:
            library_tot_reads[library_tuple][1] = library_tot_reads[library_tuple][0]

        ####

        # frac_excluded = excluded_barcode_count/tot_count
        # if frac_excluded > 5e-3:
        #     print(f'Excluded {excluded_barcode_count} reads ({frac_excluded*100:.2f}% of total) from {file}')
                #   , representing {len(excluded_barcodes)} barcodes, from {file}')

print(('E1', 3, 15, 1) in library_names) # False
print(('E1', 4, 15, 1) in library_names) # False
print(('E2', 3, 15, 1) in library_names) # True 
library_names.sort(key = lambda x: (x[0], x[2], x[1], x[-1])) #expt, day, mouse/well, seq_rep

print(f'All barcodes detected in resequecing run: {len(barcode_dict)}')
n_detected_counts = [len(barcode_dict[barcode].items()) for barcode in barcode_dict.keys()] 
tot_bc_counts = [np.sum([reads for key, reads in barcode_dict[barcode].items()]) for barcode in barcode_dict.keys()]
print(f'Num barcodes detected in >1 library: {np.sum(np.array(n_detected_counts) > 1)}')
for barcode in list(barcode_dict.keys()):#exclude barcodes detected in only 1 library
    if (len(barcode_dict[barcode].items()) < 2):
        for library in barcode_dict[barcode]:
            library_tot_reads[library][2] += barcode_dict[barcode][library]
        del barcode_dict[barcode]
for library in library_tot_reads:
    tot_reads = library_tot_reads[library][0]
    excluded_barcode_count = library_tot_reads[library][2]
    frac_excluded = excluded_barcode_count/tot_reads
    if frac_excluded > 5e-3 and tot_reads > 1e4:
        print(f'Excluded {excluded_barcode_count} reads ({frac_excluded*100:.2f}% of total) from {library}')


### organize output csvs
barcode_tots = [(barcode, np.sum([reads for key, reads in barcode_dict[barcode].items()])) for barcode in barcode_dict.keys()]
sorted_barcodes = sorted(barcode_tots, key=lambda x: x[1], reverse=True) # sort high to low, based on TOTAL ABUNDANCE ACROSS ALL LIBRARIES

all_cols = []
all_col_names = []
col_id_map = {}
# for library_names in [E1_names, E2_names, vitro_names, ]:
for key in library_names:
    all_cols.append(key)
    if key[-1] == CEC_DAY:
        col = '_'.join([str(e) for e in key[:3]])+' cecum'
    if key[-1] == INT_DAY:
        col = '_'.join([str(e) for e in key[:3]])+' intestine'
    else:
        col = '_'.join([str(e) for e in key])
    all_col_names.append(col)
    col_id_map[col] = key



lib_tots = np.zeros(len(all_cols))
lib_trimmed_tots = np.zeros(len(all_cols))

### write barcode reads to master csv
f = (args.output_dir / 'all_barcode_reads.csv').open('w')
f_trimmed = (args.output_dir / 'all_barcode_reads_trimmed.csv').open('w')
f_demo = (args.output_dir / 'all_barcode_reads_demo.csv').open('w')

header = ','.join(['barcode', 'tot_reads', 'column=(experiment_mouse/inoc_day/passage_sequencing-replicate)'] + all_col_names)
f.write(header + '\n')
f_trimmed.write(header + '\n')
f_demo.write(header + '\n')

for b, (barcode, tot_reads) in enumerate(sorted_barcodes):
    barcode_reads = []
    for lib in all_cols:
        if lib in barcode_dict[barcode]:
            barcode_reads.append( str(barcode_dict[barcode][lib]) )
        else:
            barcode_reads.append( '0' )

    barcode_line = ','.join([barcode, str(tot_reads),''] + barcode_reads)
    f.write(barcode_line + '\n')

    barcode_arr = np.array([int(e) for e in barcode_reads])
    lib_tots += barcode_arr

    barcode_nonzero = barcode_arr > 0
    if barcode_nonzero.sum() > 1: #include only lineages measured in >1 library
        f_trimmed.write(barcode_line + '\n')
        lib_trimmed_tots += barcode_arr
    if b < 100: # demo file
        f_demo.write(barcode_line + '\n')

misc_reads = [library_tot_reads[lib][2] for lib in all_cols]
misc_line = ','.join([str(library_tot_reads[lib][2]) for lib in all_cols])
f.write(f'Measured in one lib,{np.sum(misc_reads)},,'+misc_line+'\n')
lib_Dnoms =[library_tot_reads[lib][0] for lib in all_cols]
Dnom_line = ','.join([str(library_tot_reads[lib][0]) for lib in all_cols])
f.write(f'Nominal depth Dnom,{np.sum(lib_Dnoms)},,'+Dnom_line+'\n')
lib_Deffs =[library_tot_reads[lib][1] for lib in all_cols]
Deff_line = ','.join([str(library_tot_reads[lib][1]) for lib in all_cols])
f.write(f'Effective depth Deff,{np.sum(lib_Deffs)},,'+Deff_line+'\n')

f.close()
f_trimmed.close()
f_demo.close()
print('Complete.')
# print(E1_names[:10], E2_names[:10], vitro_names[:10])
