import numpy as np
import scipy
from scipy.special import gamma, digamma, loggamma, logsumexp
import argparse
import pickle
import matplotlib.pyplot as plt

import os
from methods.config import remote_paths

library_relabels = {}
with open('library_mislabelings.tsv', 'r') as f:
    header = next(f)
    for line in f:
        fastq_label_str, true_label_str = line.strip('\n').split('\t')
        fastq_label_items = fastq_label_str.strip('"').split(',')
        true_label_str = true_label_str.strip('"').split(',')

        fastq_label = (fastq_label_items[0], int(fastq_label_items[1]), int(fastq_label_items[2]))
        if true_label_str[0] == 'None':
            true_label = None
        else:
            true_label = (true_label_str[0], int(true_label_str[1]), int(true_label_str[2]))
        library_relabels[fastq_label] = true_label

noise_estimate_dir = remote_paths['noise_estimates']
pickled_dir = remote_paths['pickles']
# noise_estimate_dir = f'~/sherlock_mountbarseq_noise'

with open(f'{pickled_dir}/mouse_col_map.pkl', 'rb') as f:
    mouse_col_map = pickle.load(f)

with open(f'{pickled_dir}/vitro_col_map.pkl', 'rb') as f:
    well_col_map = pickle.load(f)

def parse_header(header):
    header_items = header.split('_')
    cohort = header_items[0] #should be E1, E2, vitro, vitro

    if cohort == 'vitro':
        if 'passage0' in header:
            well = header_items[1]
            passage = 0
        else:
            well = header_items[1] + '_' + header_items[2]
            passage = int(header_items[3].lstrip('passage'))
        return (cohort, well, passage)
        



    mouse = header_items[1].replace('M', '')
    if 'P' in mouse or 'S' in mouse or 'V' in mouse: #mouse, not day 0
        if 'delay' in header:
            mouse = header_items[1].lstrip('M')+'_delay' 
        day = 0
    else:
        mouse = int(mouse)
        if 'cec' in header:
            day = 'cecum'
        else:
            day = int(header_items[2].lstrip('D'))
    return (cohort, mouse, day)


if __name__ == '__main__':
    dir_files = os.listdir(noise_estimate_dir)
    noise_files = [f for f in dir_files if f.endswith('.out')]

    effective_noise_dict = {}

    print(well_col_map)

    for noise_file in noise_files:
        lib_id = parse_header(noise_file)
        if lib_id in library_relabels:
            if library_relabels[lib_id] is None and lib_id not in library_relabels.values():
                print('Skipping', lib_id)
                continue
        (cohort, mouse, day) = lib_id
        print(cohort, mouse ,day)
        assert (cohort, mouse, day) in mouse_col_map or (cohort, mouse, day) in well_col_map
        
        try:

            with open(f'{noise_estimate_dir}/{noise_file}', 'r') as f:
                header = next(f)
                (cohort, mouse, day) = parse_header(header)
                

    
                effective_noise_dict[(cohort, mouse, day)] = {}
                for line in f:
                    line_items = line.lstrip('\n').split('\t')
                    if len(line_items) < 2: continue

                    label = line_items[0]
                    data = line_items[1]

                    effective_noise_dict[(cohort, mouse, day)][label] = data

        except:
            print("TRY AGAIN:", noise_file)

with open(f'pickled_noise_inference.pkl', 'wb') as f:
    pickle.dump(effective_noise_dict, f)
