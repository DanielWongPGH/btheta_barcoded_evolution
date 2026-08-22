import pickle
import sys
from pathlib import Path

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[2]
NOTEBOOKS_DIR = PROJECT_ROOT / "notebooks"
if str(NOTEBOOKS_DIR) not in sys.path:
    sys.path.insert(0, str(NOTEBOOKS_DIR))

from methods.config import local_paths, pickled_dir

# combine in vivo barcode arrays
# f = open(f'./all_barcode_reads.csv', 'r')
load_dir = Path(local_paths['rebarseq'])
save_dir = Path(pickled_dir)
save_dir.mkdir(parents=True, exist_ok=True)
 
file_in = (load_dir / 'all_barcode_reads.csv').open('r')
header = next(file_in)
header_items = header.strip('\n').split(',')

bc_col_index = 0

barcodes = []
mega_array = []
library_col_names = header_items[3:] #sum E1 input library replicates
##
vivo_array, vivo_row_count = [], 0 # sequencing run 1
vitro_array, vitro_row_count = [], 0 # sequencing run 1
reseq_array, reseq_row_count = [], 0  #sequencing run 2

id_arr_map = {}
col_id_map = {}

## parse CSV
CEC_DAY = 100
INT_DAY = 200
for z, col_str in enumerate(header_items[3:]):
    col_items = col_str.split('_')
    expt, mouse, day, rep = col_items[0], col_items[1], int(col_items[2]), int(col_items[3]) 
    if mouse.isnumeric():
        mouse = int(mouse)
    # print(mouse)
    if day == CEC_DAY:
        day = 'cecum'
    if day == INT_DAY:
        day = 'smallint'
        
    if rep > 1:
        id_arr_map[(expt, mouse, day, rep)] = ('reseq', reseq_row_count)
        reseq_row_count += 1
    elif expt == 'EV':
        if '-' in mouse:
            mouse = mouse.replace('-', '_')
        id_arr_map[(expt, mouse, day, rep)] = ('vitro', vitro_row_count)
        vitro_row_count += 1
    else:
        id_arr_map[(expt, mouse, day, rep)] = ('vivo', vivo_row_count)
        vivo_row_count += 1
        
    col_id_map[z] = (expt, mouse, day, rep)

# generate barcode read arrays
seq_depths = []
for line in file_in:
    line_items = line.strip('\n').split(',')
    barcode = line_items[bc_col_index]
    reads = [int(float(item) )for item in line_items[3:]] 
    if barcode == 'Measured in one lib':
        print('Skipping reads assoc. with barcodes measured in one lib')
    elif barcode == "Nominal depth Dnom":
        seq_depths.append(reads)
    elif barcode == "Effective depth Deff":
        seq_depths.append(reads)
    else:
        barcodes.append(barcode)
        mega_array.append(reads)
file_in.close()


barcodes = np.array(barcodes)
mega_array  = np.array(mega_array).T
seq_depths = np.array(seq_depths).T

print(seq_depths.shape)



vivo_ids = [z for z, lib in col_id_map.items() if id_arr_map[lib][0] == 'vivo']
vitro_ids = [z for z, lib in col_id_map.items() if id_arr_map[lib][0] == 'vitro']
reseq_ids = [z for z, lib in col_id_map.items() if id_arr_map[lib][0] == 'reseq']


vivo_array = np.copy(mega_array[vivo_ids])
vivo_depths = np.copy(seq_depths[vivo_ids]) ## first index is nominal, second is effective depth
vitro_array = np.copy(mega_array[vitro_ids])
vitro_depths = np.copy(seq_depths[vitro_ids])
reseq_array = np.copy(mega_array[reseq_ids])
reseq_depths = np.copy(seq_depths[reseq_ids])
reseq_depths[:, 1] = reseq_depths[:, 1] * 0.3 ## see paper Fig. S1I for justification


with open(f'{save_dir}/barcodes.pkl', 'wb') as f:
    pickle.dump(barcodes, f)

for name, subarray, subdepths in [['vivo', vivo_array, vivo_depths], 
                                  ['vitro', vitro_array, vitro_depths], 
                                  ['reseq', reseq_array, reseq_depths]]:
    with open(f'{save_dir}/{name}_array.pkl', 'wb') as f:
        pickle.dump(subarray, f)

    with open(f'{save_dir}/{name}_depths.pkl', 'wb') as f:
        pickle.dump(subdepths, f)
    

print('Shape of read arrays:', vivo_array.shape, vitro_array.shape, reseq_array.shape)
#### IN VIVO ARRAY
# generate metadata mapping mouse/in vitro timecourses to array indices
# map from mouse to inoculum 
mouse_day0_map =  { 
    'E1': {1: ('E1','P2', 0), 2: ('E1','P2', 0), 3: ('E1','P2', 0), 
           4: ('E1','P2', 0), 5: ('E1','P2', 0),  6: ('E1','P1', 0), 
           7: ('E1','P1', 0), 8: ('E1','P1', 0), 9: ('E1','P1', 0), 
           10: ('E1','P1', 0), 11: ('E1','P1', 0), 12: ('E1','P1', 0), 
           13: ('E1','P1', 0)},
    'E2': {1: ('E2','S1', 0), 2: ('E2','S2', 0), 3: ('E2','S3', 0), 
           4: ('E2','S4', 0), 5: ('E2','S5', 0), 6: ('E2','S1', 0), 
           7: ('E2','S2', 0), 8: ('E2','S3', 0), 9: ('E2','S4', 0), 
           10: ('E2','S5', 0), 11: ('E2','S1delay', 0), 12: ('E2','S2delay', 0), 
           13: ('E2','S3delay', 0), 14: ('E2','S4delay', 0), 15: ('E2','S5delay', 0)}} 
day0_map = {}
for expt in ['E1', 'E2']:
    for mouse, day0_lib in mouse_day0_map[expt].items():
        if day0_lib not in day0_map:
            day0_map[day0_lib] = []
        day0_map[day0_lib].append([expt, mouse])

medium_to_well_map = {'glucose':[('p1_A1', 'V1'), ('p1_B1', 'V4'), ('p1_C1', 'V1'), ('p1_D1', 'V4')],
     'stachyose':[('p1_A2', 'V2'), ('p1_B2', 'V5'), ('p1_C2', 'V2'), ('p1_D2', 'V5')],
     'sucrose':[('p1_A3', 'V3'), ('p1_B3', 'V6'), ('p1_C3', 'V3'), ('p1_D3', 'V6')],
     'fructose':[('p1_A5', 'V1'), ('p1_B5', 'V4'), ('p1_C5', 'V1'), ('p1_D5', 'V4')],
     'maltose':[('p1_A6', 'V2'), ('p1_B6', 'V5'), ('p1_C6', 'V2'), ('p1_D6', 'V5')],
     'trehalose':[('p1_A7', 'V3'), ('p1_B7', 'V6'), ('p1_C7', 'V3'), ('p1_D7', 'V6')],
     'melibiose':[('p1_A9', 'V1'), ('p1_B9', 'V4'), ('p1_C9', 'V1'), ('p1_D9', 'V4')],
     'iron':[('p1_A10', 'V2'), ('p1_B10', 'V5'), ('p1_C10', 'V2'), ('p1_D10', 'V5')],
     'vitafiber':[('p1_A11', 'V3'), ('p1_B11', 'V6'), ('p1_C11', 'V3'), ('p1_D11', 'V6')],
     'raffinose':[('p1_E1', 'V1'), ('p1_F1', 'V4'), ('p1_G1', 'V1'), ('p1_H1', 'V4')],
     'lactose':[('p1_E2', 'V2'), ('p1_F2', 'V5'), ('p1_G2', 'V2'), ('p1_H2', 'V5')],
     'galactose':[('p1_E5', 'V1'), ('p1_F5', 'V4'), ('p1_G5', 'V1'), ('p1_H5', 'V4')],
     
     'UMich-01':[('p2_A1', 'V1'), ('p2_B1', 'V4'), ('p2_C1', 'V1'), ('p2_D1', 'V4')],
     'bioecolians':[('p2_A2', 'V2'), ('p2_B2', 'V5'), ('p2_C2', 'V2'), ('p2_D2', 'V5')],
     'prebiotin':[('p2_A3', 'V3'), ('p2_B3', 'V6'), ('p2_C3', 'V3'), ('p2_D3', 'V6')],
     'maltodextrin':[('p2_A5', 'V1'), ('p2_B5', 'V4'), ('p2_C5', 'V1'), ('p2_D5', 'V4')],
     'promitor':[('p2_A6', 'V2'), ('p2_B6', 'V5'), ('p2_C6', 'V2'), ('p2_D6', 'V5')],
     'vitagos':[('p2_A7', 'V3'), ('p2_B7', 'V6'), ('p2_C7', 'V3'), ('p2_D7', 'V6')],
     'precticx':[('p2_A9', 'V1'), ('p2_B9', 'V4'), ('p2_C9', 'V1'), ('p2_D9', 'V4')],
     'LC742':[('p2_A10', 'V2'), ('p2_B10', 'V5'), ('p2_C10', 'V2'), ('p2_D10', 'V5')],
     'avantafiber':[('p2_A11', 'V3'), ('p2_B11', 'V6'), ('p2_C11', 'V3'), ('p2_D11', 'V6')],

     'UMich-02b':[('p2_E1', 'V1'), ('p2_F1', 'V4'), ('p2_G1', 'V1'), ('p2_H1', 'V4')],
     'yacontrol':[('p2_E2', 'V2'), ('p2_F2', 'V5'), ('p2_G2', 'V2'), ('p2_H2', 'V5')],
     'fibersol':[('p2_E5', 'V1'), ('p2_F5', 'V4'), ('p2_G5', 'V1'), ('p2_H5', 'V4')],
     'cravingzgone':[('p2_E6', 'V2'), ('p2_F6', 'V5'), ('p2_G6', 'V2'), ('p2_H6', 'V5')],
     'bimuno':[('p2_E7', 'V3'), ('p2_F7', 'V6'), ('p2_G7', 'V3'), ('p2_H7', 'V6')],
     'sunfiber':[('p2_E9', 'V1'), ('p2_F9', 'V4'), ('p2_G9', 'V1'), ('p2_H9', 'V4')],
     'ISOT 160120':[('p2_E10', 'V2'), ('p2_F10', 'V5'), ('p2_G10', 'V2'), ('p2_H10', 'V5')],
     'wako':[('p2_E11', 'V3'), ('p2_F11', 'V6'), ('p2_G11', 'V3'), ('p2_H11', 'V6')]
     }
well_to_medium_map = {}
well_to_inoc_map = {}
for medium, wells in medium_to_well_map.items():
    for well in wells:
        well_to_medium_map[well[0]] = (medium, well[1])
        day0_lib = ('EV', well[1], 0)
        well_to_inoc_map[well[0]] = day0_lib #should be value in vitro_day0_map
        if day0_lib not in day0_map:
            day0_map[day0_lib] = []
        day0_map[day0_lib].append(['EV', well[0]])


# get cols for different arrays
# mouse_col_map, col_mouse_map = {}, {}
library_array_map = {}
timecourse_meta = {'E1':{}, 'E2':{}, 'EV':{}}
col_to_id_map = {'vivo':{}, 'vitro':{}, 'reseq':{}}


for i, col_name in enumerate(library_col_names):
    library = col_id_map[i]

    if library[1] in ['p2_E3', 'p2_F3', 'p2_G3', 'p2_H3']: #dummy wells
        continue
    lib_norep = library[:-1]
    if lib_norep in day0_map: # day 0 library, append to all relevant mouse/well library
        # print('DAY 0 lib!!', lib_norep, col_name)
        day0 = 0

        if 'delay' in lib_norep[1]:
            day0 = 14

        arr_name, arr_idx = id_arr_map[library]
        col_to_id_map[arr_name][arr_idx] = lib_norep



        day0_lib = (library[0], library[1])
        if day0_lib not in library_array_map:
            library_array_map[day0_lib] = []
        library_array_map[day0_lib].append( (arr_name, arr_idx) )

        for (expt, mouse) in day0_map[lib_norep]:
            if mouse not in timecourse_meta[expt]:
                timecourse_meta[expt][mouse] = {}
            if day0 not in timecourse_meta[expt][mouse]:
                timecourse_meta[expt][mouse][day0] = []
            timecourse_meta[expt][mouse][day0].append( (arr_name, arr_idx) )

            if (expt, mouse, day0) not in library_array_map:
                library_array_map[(expt, mouse, day0)] = []    
            library_array_map[(expt, mouse, day0)].append( (arr_name, arr_idx) )
        # print(lib_norep)
        

    else:
        expt, mouse, day = lib_norep
        if (expt, mouse, day) not in library_array_map:
            library_array_map[(expt, mouse, day)] = []

        arr_name, arr_idx = id_arr_map[library]
        col_to_id_map[arr_name][arr_idx] = lib_norep

        library_array_map[(expt, mouse, day)].append( (arr_name, arr_idx) )
        if mouse not in timecourse_meta[expt]:
            timecourse_meta[expt][mouse] = {}
        if day not in timecourse_meta[expt][mouse]:
            timecourse_meta[expt][mouse][day] = []
        timecourse_meta[expt][mouse][day].append( (arr_name, arr_idx) )


with open(f'{save_dir}/timecourse_meta.pkl', 'wb') as f:
    pickle.dump(timecourse_meta, f)

with open(f'{save_dir}/library_array_map.pkl', 'wb') as f:
    pickle.dump(library_array_map, f)

with open(f'{save_dir}/medium_to_well_map.pkl', 'wb') as f:
    pickle.dump(medium_to_well_map, f)

with open(f'{save_dir}/well_to_medium_map.pkl', 'wb') as f:
    pickle.dump(well_to_medium_map, f)

with open(f'{save_dir}/col_id_map.pkl', 'wb') as f:
    pickle.dump(col_to_id_map, f)

barcode_overlap_map = {'V1':['P1', 'S1'],
                       'V2':['P1', 'S1', 'S2'],
                       'V3':['P1', 'S2'],
                       'V4':['P2', 'S3'],
                       'V5':['P2', 'S3', 'S4'],
                       'V6':['P2', 'S4'],
                       'P1':['S1', 'S2', 'S5', 'V1', 'V2', 'V3'],
                       'P2':['S3', 'S4', 'V4', 'V5', 'V6'],
                       'S1':['P1', 'V1', 'V2'],
                       'S2':['P1', 'V2', 'V3'],
                       'S3':['P2', 'V4', 'V5'],
                       'S4':['P2', 'V5', 'V6'],
                       'S5':['P1', 'S1', 'S2', 'V1', 'V2', 'V3']
                       }

with open(f'{save_dir}/barcode_overlap_map.pkl', 'wb') as f:
    pickle.dump(barcode_overlap_map, f)

    
