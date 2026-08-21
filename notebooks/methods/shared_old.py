import pickle
import numpy as np
import sys; sys.path.insert(0, '../../')
from methods.config import *
import matplotlib.pyplot as plt
import matplotlib as mpl
import scipy
import methods.sim_barcodes as sim


### load vivo array dict
with open(f'../data/rebarseq/pickled/vivo_array.pkl', 'rb') as f:
    vivo_array = pickle.load(f)
# freq_array = np.einsum('ij, i->ij', read_array, read_array.sum(axis=1)**-1.)
with open(f'../data/rebarseq/pickled/barcodes.pkl', 'rb') as f:
    barcodes = pickle.load(f)

with open(f'../data/rebarseq/pickled/mouse_meta.pkl', 'rb') as f:
    mouse_meta = pickle.load(f)

with open(f'../data/rebarseq/pickled/mouse_col_map.pkl', 'rb') as f:
   vivo_row_ids = pickle.load(f)

with open(f'../data/rebarseq/pickled/vivo_array_pseudoderep.pkl', 'rb') as f:
    vivo_array_pseudo = pickle.load(f)

with open(f'../data/rebarseq/pickled/vivo_array_meta.pkl', 'rb') as f:
    vivo_meta_pseudo = pickle.load(f)

with open(f'../data/rebarseq/pickled/mouse_meta_pseudoderep.pkl', 'rb') as f:
    mouse_meta_pseudo = pickle.load(f)


Deff_map = vivo_meta_pseudo['Deff_map']
Deff_array_pseudo = vivo_meta_pseudo['Deff_array']


# with open(f'{data_dir}/pickled/barcode_read_array.pkl', 'rb') as f:
#     pseudo_read_array = pickle.load(f)
# freq_array = np.einsum('ij, i->ij', read_array, read_array.sum(axis=1)**-1.)

# with open(f'{data_dir}/pickled/barcodes.pkl', 'rb') as f:
#     pseudo_barcodes = pickle.load(f)

# with open(f'{data_dir}/pickled/mouse_meta.pkl', 'rb') as f:
#     pseudo_mouse_meta = pickle.load(f)

# with open(f'{data_dir}/pickled/mouse_col_map.pkl', 'rb') as f:
#     mouse_col_map = pickle.load(f)

# with open(f'{data_dir}/pickled/vitro_meta.pkl', 'rb') as f:
#     well_meta = pickle.load(f)

# with open(f'{data_dir}/pickled/vitro_row_ids.pkl', 'rb') as f:
#     well_col_map = pickle.load(f)

# with open(f'{data_dir}/pickled/medium_to_well_map.pkl', 'rb') as f:
#     medium_to_well_map = pickle.load(f)

# with open(f'{data_dir}/pickled/well_to_medium_map.pkl', 'rb') as f:
#     well_to_medium_map = pickle.load(f)

with open(f'../data/pickled/barcode_pool_assignments.pkl', 'rb') as f:
    barcode_pool_assignments = pickle.load(f)

with open(f'{data_dir}/pickled/barcode_pool_map.pkl', 'rb') as f:
    index_pool_map = pickle.load(f)

with open(f'{data_dir}/pickled/barcode_overlap_map.pkl', 'rb') as f:
    barcode_overlap_map = pickle.load(f)

### 
def get_mouse_timecourse(expt, mouse, include_pseudo_derepped=False):
    ##
    mouse_timepoints = []
    mouse_reads = []
    mouse_depths = []
    mouse_contamination_rates = []

    derepped_dict = mouse_meta[expt][mouse]
    derepped_days = sorted([key for key in derepped_dict.keys() if type(key)==int])
    pseudo_days = np.array(mouse_meta_pseudo[expt][mouse][0])
    pseudo_indices = np.array(mouse_meta_pseudo[expt][mouse][1])
    
    for day in derepped_days:
        replicates = derepped_dict[day]
        replicate_arr = vivo_array[replicates]
        replicate_depths = replicate_arr.sum(axis=1)
        if day in pseudo_days:
            if replicate_depths[0] < 10**4:
                if len(replicate_arr) == 1:
                    continue
                else:
                    summed_reads = replicate_arr[1:].sum(axis=0)
            else:
                summed_reads = replicate_arr.sum(axis=0)
        
            # pseudo_index = pseudo_indices[ np.where(day == pseudo_days)[0][0] ]
            # pseudo_depth = Deff_array_pseudo[pseudo_index]
            # pseudo_reads = vivo_array_pseudo[pseudo_index]
            # if pseudo_reads.sum() > 10*pseudo_depth: ## UMI issues
            #     if not include_pseudo_derepped:
            #         if len(replicates) > 1:
            #             summed_reads = replicate_arr[1:].sum(axis=0)
            #         else:
            #             continue  
            #     elif include_pseudo_derepped and len(replicates) == 1:
            #         summed_reads = vivo_array_pseudo[pseudo_index]
            #     elif include_pseudo_derepped and len(replicates) > 1:
            #         collated_reads = [replicate_arr[1:]]
        else:
            where_well_measured = replicate_depths > 10**4
            if where_well_measured.sum() == 0:
                continue
            summed_reads = replicate_arr[where_well_measured].sum(axis=0)
        # else: 
        #     summed_reads = replicate_arr[replicates].sum(axis=0)                
        mouse_timepoints.append(day)
        mouse_reads.append(summed_reads)
        mouse_depths.append(summed_reads.sum())

    # if include_pseudo_derepped:
    #     for day in pseudo_days:
    #         if day not in mouse_timepoints:
    #             mouse_timepoints.

    return np.array(mouse_timepoints), np.array(mouse_reads), np.array(mouse_depths)

    #


### BtVPI GENOME INFO
GENOME_LENGTH = 6260361

LE1 = 'GCCTTATAAATCTGGCTCTT'
LE2 = 'TGCTTCCGGCTTGGAAACCG'

CPS_LOCATIONS = {'cps1': (462368, 489867, +1), #start, stop, strand
                 'cps2': (573167, 595424, +1),
                 'cps3': (734869, 761399, +1),
                 'cps4': (1663708, 1685198, -1),
                 'cps5': (2032185, 2047991, -1),
                 'cps6': (2104512, 2124231, -1),
                 'cps7': (3580840, 3605962, -1),
                 'cps8': (32315, 71500, +1)}


### IN VITRO ASIDES

MONOSACCHARIDES = ['glucose', 'iron', 'fructose', 'galactose']
DISACCHARIDES = ['lactose', 'maltose', 'melibiose', 'sucrose', 'trehalose'] #['sucrose', 'maltose', 'trehalose', 'melibiose', ]
OLIGOSACCHARIDES = ['raffinose', 'stachyose']
POLYSACCHARIDES = ['avantafiber', 'bioecolians', 'bimuno', 'cravingzgone', 'fibersol',
                   'ISOT 160120', 'LC742', 'maltodextrin', 'prebiotin', 'precticx', 
                   'promitor', 'sunfiber', 'UMich-01', 'UMich-02b', 'vitafiber',
                   'vitagos', 'wako', 'yacontrol']


ordered_media = MONOSACCHARIDES + DISACCHARIDES + OLIGOSACCHARIDES + POLYSACCHARIDES

VITRO_MEDIA_COLORS = {'glucose': 'red', 'lactose': 'blue', 'vitafiber':'green'} \
                | {medium:KELLY_COLORS[m%len(KELLY_COLORS)] for m, medium in enumerate(ordered_media)  if medium not in ['glucose', 'lactose', 'vitafiber']}


medium_labels = {medium:medium for medium in ordered_media}
medium_labels['iron'] = 'gluc. + iron'


gene_coords_map = {}
gene_description = {}
gene_starts, gene_stops = [], []
gene_lst = []
with open(f'{data_dir}/reference_genome/BtVPI.ptt', 'r') as f:
    for _ in range(4):
        header = next(f)
    for line in f:
        line_items = line.strip('\n').split('\t')
        loc_str, strand = line_items[0], line_items[1]
        gene, description = line_items[5], line_items[-1]

        start, stop = [int(e) for e in loc_str.split('..')]
        gene_coords_map[gene] = (start, stop, strand)
        gene_description[gene] = description
        gene_starts.append(start)
        gene_stops.append(stop)
        gene_lst.append(gene)
gene_starts, gene_stops = np.array(gene_starts), np.array(gene_stops)
genes = np.array(gene_lst)

def find_gene_from_position(position):
    jc_start = position >= gene_starts
    jc_stop = position <= gene_stops

    if (jc_start * jc_stop).sum() > 0:
        return genes[ jc_start * jc_stop ]
    else: # intergenic
        if np.any(~jc_start) and np.any(~jc_stop):
            left_gene = genes[~jc_stop][-1]
            right_gene = genes[~jc_start][0]
        elif np.any(~jc_start):
            right_gene = genes[~jc_start][0]
            left_gene = genes[-1]
        elif np.any(~jc_stop):
            right_gene = genes[~jc_stop][-1]
            left_gene = genes[0]

        return [f'{left_gene}-{right_gene}']

##########################################
### BARCODE FREQ CALCULATIONS  ####
##########################################

def calc_freqs(read_array):
    if read_array.ndim > 1:
        return np.einsum('ij, i->ij', read_array, read_array.sum(axis=1)**-1.)
    else:
        return read_array / read_array.sum()

def maxmin_freqs(freqs0, depth0, freqs1, depth1):
    if type(freqs0) != np.ndarray:
        freqs0 = np.array(freqs0)
        freqs1 = np.array(freqs1)

    overD0 = np.full(freqs0.shape[-1], 1 / depth0)
    overD1 = np.full(freqs1.shape[-1], 1 / depth1)

    # max_f0 = np.max([freqs0, overD0], axis=0)
    # max_f0[freqs0 == 0] = 0
    # min_f0 = np.min([freqs1, overD0], axis=0)
    # maxmin_f0 = np.max([max_f0, min_f0], axis=0)

    # max_f1 = np.max([freqs1, overD1], axis=0)
    # max_f1[freqs1 == 0] = 0
    # min_f1 = np.min([freqs0, overD1], axis=0)
    # maxmin_f1 = np.max([max_f1, min_f1], axis=0)

    min_f0 = np.min([freqs1, overD0], axis=0)
    maxmin_f0 = np.max([freqs0, min_f0], axis=0)

    min_f1 = np.min([freqs0, overD1], axis=0)
    maxmin_f1 = np.max([freqs1, min_f1], axis=0)

    return maxmin_f0, maxmin_f1
def get_freqs(expt, mouse, day):
    row_id = vivo_row_ids[(expt, mouse, day)]
    freqs = freq_array[row_id]
    depth = read_array[row_id].sum()
    effective_depth = Deff_array[row_id]
    return freqs, depth

def cage_avg_freqs(expt, mice, day):
    weighted_read_array = np.zeros(read_array.shape[-1])
    weighted_depth = 0

    for mouse in mice:
        if (expt, mouse, day) in vivo_row_ids:
            row_id = vivo_row_ids[(expt, mouse, day)]
            mouse_freqs = freq_array[row_id]
            mouse_D = Deff_array[row_id]
            mouse_reads = mouse_freqs * mouse_D

            weighted_read_array += mouse_reads * mouse_D / mouse_reads.sum()
            weighted_depth += mouse_D

    weighted_freq_array  = weighted_read_array/weighted_depth
    return weighted_freq_array, weighted_depth

##########################################
### MULLER PLOTS ####
##########################################
highlight_barcode_hatches = {'GCCTTATAAATCTGGCTCTT': '\\\\\\\\', 'TGCTTCCGGCTTGGAAACCG': '////'}

def find_large_barcodes(expt, mice, max_cutoff=1e-2, final_cutoff=1,):
    all_large_indices = set()
    for mouse in mice:
        
        # mouse_rows = mouse_meta[expt][mouse][1]
        # mouse_depths = Deff_array[mouse_rows]
        # # max_freqs = np.max(freq_array[mouse_rows][mouse_depths**-1. < max_cutoff*1/10], axis=0) #exclude day 0
        # # max_freqs = np.max(freq_array[mouse_rows][mouse_depths**-1. < max_cutoff], axis=0) #exclude day 0
        # max_freqs = np.max(freq_array[mouse_rows], axis=0) #exclude day 0


        # final_timepoint_row = mouse_rows[-1]
        # final_freqs = freq_array[final_timepoint_row]

        mouse_days, mouse_reads, mouse_depths =  get_mouse_timecourse(expt, mouse)
        mouse_freqs = calc_freqs(mouse_reads)
        max_freqs = np.max(mouse_freqs, axis=0) #exclude day 0
        final_freqs = mouse_freqs[-1]

        large_indices = np.where( (max_freqs > max_cutoff) + (final_freqs > final_cutoff) > 0)[0]
        all_large_indices.update( large_indices )

    all_large_indices = np.array(list(all_large_indices))
    return all_large_indices

def order_barcodes_by_pool(large_barcodes, barcode_pool_assignments,  ordered_pool_color_tuples, num_colors_in_pool=18, cg_at_end=True, cg_all_low_freqs=False):
    ordered_indices = []
    ordered_colors = []

    cg_lst = []
    cg_colors = []

    for pool, pool_colormap in ordered_pool_color_tuples:
        pool_barcodes = barcode_pool_assignments[pool][0]
        pool_cmap = mpl.colormaps[pool_colormap]
        color_run = list(iter(pool_cmap(np.linspace(0.5, 1.0, num_colors_in_pool))))
        color_run = color_run[0::3] + color_run[1::3] + color_run[2::3] # scramble adjacent colors


        small_barcode_indices = []
        c = 0
        for barcode in pool_barcodes:
            if [barcode] in ordered_indices:
                continue
            if barcode in large_barcodes:
                ordered_indices.append( [barcode] )
                ordered_colors.append(color_run[c % num_colors_in_pool])
                c += 1
            else:
                small_barcode_indices.append(barcode)

        
        if cg_all_low_freqs: #single cg_lst
            cg_lst.extend( small_barcode_indices )
        elif cg_at_end:
            cg_lst.append( small_barcode_indices )
            cg_colors.append( pool_cmap(0.35) )
        else:
            ordered_indices.append( small_barcode_indices )
            ordered_colors.append( pool_cmap(0.35) )

    if cg_all_low_freqs:
        ordered_indices.append(cg_lst)
        ordered_colors.append('lightgrey')
    elif cg_at_end:
        ordered_indices.extend(cg_lst)
        ordered_colors.extend(cg_colors)

    return ordered_indices, ordered_colors 

# def make_muller_freqs(expt, mouse, sorted_indices, exclude_lowdepth=1000):
#     days = mouse_meta[expt][mouse][0]
#     rows = mouse_meta[expt][mouse][1]

#     included_days = []
#     freqs = []
    
#     for day, row in zip(days, rows):
#         day_freqs = []
#         mouse_freqs = freq_array[row]
#         mouse_depth = read_array[row].sum()

#         mouse_Deff = Deff_array[row]
        
#         if mouse_Deff < exclude_lowdepth: #if above threshold
#             continue
        
#         included_days.append(day)
#         for lineage_indices in sorted_indices:
#             day_freqs.append( mouse_freqs[lineage_indices].sum() ) 
#         day_freqs.append( 1 - np.sum(day_freqs) ) #leftover

#         freqs.append(day_freqs)

#     return np.array(included_days), np.array(freqs)

def make_muller_freqs(expt, mouse, sorted_indices, exclude_lowdepth=1000):
    days, mouse_reads, mouse_depths = get_mouse_timecourse(expt, mouse)

    print(days, mouse_reads[:10, :10])
    mouse_freqs = calc_freqs(mouse_reads)

    included_days = []
    muller_freqs = []
    
    for day, freqs, depth in zip(days,  mouse_freqs, mouse_depths):
        day_freqs = []
        
        if depth < exclude_lowdepth: #if above threshold
            continue
        
        included_days.append(day)
        for lineage_indices in sorted_indices:
            day_freqs.append( freqs[lineage_indices].sum() ) 
        day_freqs.append( 1 - np.sum(day_freqs) ) #leftover
        muller_freqs.append(day_freqs)


    # print(included_days)
    return np.array(included_days), np.array(muller_freqs)



def make_vitro_muller_freqs(passages, freqs, sorted_indices, exclude_lowdepth=1000):
    included_passages = []
    muller_freqs = []
    for passage, passage_freqs in zip(passages, freqs):
        stacked_freqs = []
      
        included_passages.append(passage)
        for lineage_indices in sorted_indices:
            stacked_freqs.append( passage_freqs[lineage_indices].sum() ) 
        stacked_freqs.append( 1 - np.sum(stacked_freqs) ) #leftover

        muller_freqs.append(stacked_freqs)

    return np.array(included_passages), np.array(muller_freqs)

def muller_plot(ax, times, bc_indices, freqs, colors, hatch_dict={}):
    """ freqs has dimensions of (time, bc freq).
     times and freqs.shape[0] should have same dimension. """
    cum_freqs = np.cumsum(freqs, axis=1)

    # ax.fill_between(times, cum_freqs[:, i], cum_freqs[:, i+1], facecolor=colors[i + 1])
    
    ## first lineage
    try:
        hatch = hatch_dict.get(barcodes[bc_indices[0][0]], None)
    except:
        hatch = None
    ax.fill_between(times, 0, cum_freqs[:, 0], facecolor=colors[0], lw=0, rasterized=True, hatch=hatch) #first lineage
    for i in range(1, cum_freqs.shape[1]-1):
        try:
            hatch = hatch_dict.get(barcodes[bc_indices[i][0]], None)
        except:
            hatch = None
        ax.fill_between(times, cum_freqs[:, i-1], cum_freqs[:, i], facecolor=colors[i], lw=0, rasterized=True, hatch=hatch)

    ax.fill_between(times, cum_freqs[:, -2], cum_freqs[:, -1], color='lightgrey', lw=0, rasterized=True)
    # ax.plot(times, cum_freqs[:, -6], color='black', lw=2)

##########################################
### DIVERSITY STATISTICS (ENTROPY, NUM. BARCODES DETECTED, MEAN-FITNESS) ####
##########################################

def coarse_grain(sorted_freqs, max_groups=3, min_freq=10**-2, recursive=False):
    groupings = [[] for i in range(max_groups)]
    current_freqs = [sorted_freqs[group].sum() for group in groupings]

    for i, f in enumerate(sorted_freqs):
        z = np.argmin(current_freqs)
        groupings[z].append(i)
        current_freqs = [sorted_freqs[group].sum() for group in groupings]

    if np.min(current_freqs) > min_freq:
        return groupings, np.array(current_freqs)

    else:
        if recursive:
            return coarse_grain(sorted_freqs, max_groups=len(groupings)-1, min_freq=min_freq)
        else:
            return groupings, np.array(current_freqs)

def high_frequency_expectation(sampled_days, freqs, high_frequency_indices, Deff_array, Ntau=10**5, min_day=0):
    low_frequency_indices = np.array([i not in high_frequency_indices for i in np.arange(freqs.shape[1])])
    high_frequency_sum = freqs[:, high_frequency_indices].sum(axis=1)

    ## convert high_frequency_sum to an expected mean fitness over time
    int_Xbar = -np.log((1-high_frequency_sum) / (1-high_frequency_sum[0]))
    Xbar_in_interval = np.diff(int_Xbar) / np.diff(sampled_days) # assume Xbar is roughly constant in a sampled interval (should be true for small intervals)

    Xbar_lst = []
    for day in range(min_day, np.max(sampled_days)):
        Xbar_lst.append( Xbar_in_interval[ np.where(sampled_days <= day)[0][-1] ] )

    sampled_array = [freqs[0]]
    sim_low_freqs = freqs[0, low_frequency_indices]

    integrated_Xbar = 1
    for i, day in enumerate(range(min_day+1, np.max(sampled_days)+1)):
        Xbar = Xbar_lst[i]

        sim_low_freqs =  sim.rnd.poisson(sim_low_freqs * np.exp(-Xbar) * Ntau) / Ntau 
        integrated_Xbar *= np.exp(-Xbar)
        if day in sampled_days:
            day_index = np.where(sampled_days == day)[0][0]

            expected_array = np.copy(freqs[day_index])
            expected_array[low_frequency_indices] = sim_low_freqs * (1-expected_array[high_frequency_indices].sum()) / sim_low_freqs.sum()

            sampled_reads = sim.rnd.poisson(expected_array * Deff_array[ day_index ]) 
            sampled_freqs = sampled_reads / sampled_reads.sum()

            # sampled_reads = sim.rnd.poisson(sampled_freqs * Deff_array[ day_index ]) 
            # sampled_freqs = sampled_reads / sampled_reads.sum()

            sampled_array.append(sampled_freqs)
    return np.array(sampled_array)

def calc_median_lfc(freq_array, Deff_array, barcode_groups): #day0_freqs = ref_freqs
    ref_freqs, D0 = freq_array[0], Deff_array[0]

    median_lfc_lst = []
    for freqs_t, Dt in zip(freq_array[1:], Deff_array[1:]):
        median_lfcs = []


        for group_indices in barcode_groups:
            f1_group = ref_freqs[group_indices]
            f2_group = freqs_t[group_indices]
            if np.percentile(f2_group, 50) == 0: #if most barcodes in group are unmeasured, exclude
                continue 
                    
            f1_maxmin, f2_maxmin = maxmin_freqs(f1_group, D0, f2_group, Dt)
            lfc = np.log( f2_maxmin / f1_maxmin )
            median_lfcs.append( [np.median(lfc), lfc] )


            # lfcs_in_group = lfc[group_indices]
            # noise_max = np.log( 1 / (Dt * np.min(ref_freqs[group_indices])) )  
            
            # median_lfc = np.median( lfcs_in_group )
            # if median_lfc - noise_max > 0: # if at least half of barcodes in group are measured --> THIS WILL OVERESTIMATE LFC (UNDERESTIMATE MEAN FITNESS)
            #     median_lfcs.append( median_lfc )
        median_lfc_lst.append( median_lfcs )

    return median_lfc_lst

def calc_entropy(freq_array, min_freq=1e-4):
    entropy_lst = []
    for freqs in freq_array:
        where_large = freqs > min_freq
        if where_large.sum() == 0:
            entropy_lst.append(0)
        else:
            large_freqs = freqs[where_large] / freqs[where_large].sum()
            entropy = scipy.stats.entropy(large_freqs)
            entropy_lst.append(entropy)

    return np.array(entropy_lst)

def calc_measured_frequencies(true_frequencies, depths):
    measured_frequencies = []
    for true_freqs, depth in zip(true_frequencies, depths):
        measured_counts = sim.rnd.poisson(depth * true_freqs)
        measured_frequencies.append( measured_counts / np.sum(measured_counts) )
    return np.array( measured_frequencies )

def estimate_num_barcodes_above_cutoff(expt, mouse, cutoff, min_depth=10**3):
    mouse_days, mouse_rows, _ = mouse_meta[expt][mouse]
    mouse_floors = Deff_array[mouse_rows]**-1.
    day0_freqs = freq_array[mouse_rows[0]]

    valid = (mouse_floors < np.min([cutoff, min_depth**-1.]))
    freqs = freq_array[mouse_rows][valid]
    n_barcodes_above_cutoff = (freqs > cutoff).sum(axis=1)

    return mouse_days[valid], n_barcodes_above_cutoff

def estimate_num_barcodes_above_cutoff_simulated(expt, mouse, cutoff, Ntau=10**6, min_depth=10**3):
    mouse_days, mouse_rows, _ = mouse_meta[expt][mouse]
    mouse_floors = Deff_array[mouse_rows]**-1.
    
    day0_freqs = freq_array[mouse_rows[0]]
    neutral_sim_freqs = sim.neutral_sim(day0_freqs, Ntau, mouse_days, mouse_floors**-1.)
    valid = (mouse_floors < np.min([cutoff, min_depth**-1.])) 
    freqs = neutral_sim_freqs[valid]
    n_barcodes_above_cutoff = (freqs > cutoff).sum(axis=1)

    return mouse_days[valid], n_barcodes_above_cutoff


##########################################
#### TRANSMISSION MEASUREMENTS ####
##########################################


def conditioning(donor_freqs, recip_freqs0, donor_min=0, recip_max=0, pool_assignment_bool=None):
    valid_bool = (donor_freqs >= donor_min) * (recip_freqs0 <= recip_max)
    if np.any(pool_assignment_bool) != None:
        valid_bool = valid_bool * pool_assignment_bool
    return valid_bool

# functions for prob. detect-based inference of transmission
def calc_p_detected(donor_freqs, recip_freqs, recip_floor=0, binspace=np.logspace(-5,0,51), pseudocount=0):
    density = []
    errs = []

    recip_floor = np.max([1e-8, recip_floor])
    for b, bin in enumerate(binspace[:-1]):
        where_in_bin = (donor_freqs >= bin) * (donor_freqs < binspace[b+1])

        where_nonzero = (recip_freqs[where_in_bin] > recip_floor)
        if where_in_bin.sum() >= 1:
            n = where_in_bin.sum()
            p = where_nonzero.sum() / n
            if p == 0:
                p = pseudocount/n
            density.append( p)
            if p*(1-p) > 0:
                errs.append( np.sqrt(p*(1-p)/n ) )
            else:
                errs.append( np.sqrt(1/n*(1-1/n)/n ) )
        else:
            density.append(np.nan)
            errs.append( np.nan)
    return np.array(density), np.array(errs), 1/2*(binspace[:-1] + binspace[1:])

def detection_neg_loglikelihood(Dm, f_d, f_r):
        f_present = (f_r > 0).astype(int)
        term1 = (-f_d * Dm)*(1-f_present)
        term2 = np.log(1-np.exp(-f_d * Dm))*f_present
        return -(term1+term2).sum()

def detection_neg_ll_jac(Dm, f_d, f_r):
    f_present = (f_r > 0).astype(int)
    numerator = f_d**2 * f_present
    denominator = 2*(np.cosh(Dm * f_d)-1)
    return (numerator/denominator).sum()

def migration_inference_by_pdetected(donor_data, recip_data0, recip_data1, valid_bool, recip_floor=0, binspace=np.logspace(-5,0,51), pseudocount=0):
    donor_freqs, donor_depth = donor_data
    recip_freqs0, recip_depth0 = recip_data0
    recip_freqs1, recip_depth1 = recip_data1
    
    density, errs, bins = calc_p_detected(donor_freqs[valid_bool], recip_freqs1[valid_bool],
                                          recip_floor=recip_floor, binspace=binspace, pseudocount=pseudocount)
    notnan = ~np.isnan(density) * (density > 0)

    if np.sum(recip_freqs1[valid_bool]) == 0: #no detected barcodes --> find m such that prob of zero measured barcodes == 1/2
        m = np.log(2) / np.sum(donor_freqs[valid_bool]*recip_depth1)
        print('no detected barcodes, m=',  m)
        return m, 0, density, errs, bins

    if notnan.sum() < 2:
        return np.nan, np.nan, density, errs, bins

    res = scipy.optimize.minimize(detection_neg_loglikelihood, x0 = 1, args=((donor_freqs[valid_bool], recip_freqs1[valid_bool])))
    m = res.x[0]/recip_depth1
    uncertainty = detection_neg_ll_jac(res.x[0], donor_freqs[valid_bool], recip_freqs1[valid_bool]) **-0.5 / recip_depth1

    return m, uncertainty, density, errs, bins

# functions for weighted average rate of transmission
def coarse_grain(sorted_freqs, max_groups=3, min_freq=10**-2):
    groupings = [[] for i in range(max_groups)]
    current_freqs = [sorted_freqs[group].sum() for group in groupings]

    for i, f in enumerate(sorted_freqs):
        z = np.argmin(current_freqs)
        groupings[z].append(i)
        current_freqs = [sorted_freqs[group].sum() for group in groupings]

    return groupings, np.array(current_freqs)

def migration_inference_by_ratio(donor_data, recip_data0, recip_data1, valid_bool, recip_floor=0, max_groups=2, min_freq=10**-2):
    donor_freqs, donor_depth = donor_data
    recip_freqs0, recip_depth0 = recip_data0
    recip_freqs1, recip_depth1 = recip_data1

    donor_valid = np.copy(donor_freqs[valid_bool])
    recip_valid = np.copy(recip_freqs1[valid_bool])
    recip_valid[recip_valid < recip_floor] = 0
    sorting = np.argsort(donor_valid)

    cg_donor, cg_recip = [], []
    for max_group in range(2, max_groups+1, 1):
        donor_coarse_graining = coarse_grain(donor_valid[sorting], max_groups=max_group, min_freq=min_freq)
        cg_groups, cg_donor_freqs = donor_coarse_graining

        cg_recipient_freqs1 = np.array([recip_valid[sorting][cg_group].sum() for cg_group in cg_groups])
        coarse_grained_ratios = cg_recipient_freqs1 / cg_donor_freqs

        cg_donor.extend(list(cg_donor_freqs))
        cg_recip.extend(list(cg_recipient_freqs1))

        if max_group == 2:
            m = np.median(coarse_grained_ratios)

    return m, (donor_valid, np.array(cg_donor)), (recip_valid, np.array(cg_recip)), coarse_grained_ratios

def p_transmit(f, m):
    return 1 - np.exp(-m*f)

def find_consecutive_timepoints(expt, donor_mice, recipient_mouse, min_day=11, max_dt=2, max_day=55):
    consecutive_sampled_days = set()
    for day in range(min_day, max_day+1): # first timepoint
        for mouse in donor_mice:
            if (expt, mouse, day) in vivo_row_ids:
                for dt in range(1,max_dt+1): # second timepoint = day + dt
                    if (expt, recipient_mouse, day) in vivo_row_ids and (expt, recipient_mouse, day + dt) in vivo_row_ids:
                        consecutive_sampled_days.add( (day, day+dt) )
    return sorted(consecutive_sampled_days, key=lambda x: x[0])

def measure_transmission_in_interval(expt, recipient_mouse, donor_mice, t0, t1, 
                                     pool_assignment_bool=None, coarse_grain_num_groups=2, min_barcodes=10,
                                     rescale_migration = 1 ):

    donor_ids = [vivo_row_ids[(expt,m, t0)] for m in donor_mice if (expt, m, t0) in vivo_row_ids]
    donor_freqs = np.mean( freq_array[donor_ids], axis=0)
    donor_floor = np.min([1 / Deff_array[donor_ids]] )


    recip_freqs_t0, recip_freqs_t1 = freq_array[ [vivo_row_ids[(expt,recipient_mouse, t)] for t in (t0, t1)] ]
    recip_floor_t0, recip_floor_t1 =  1 / Deff_array[ [vivo_row_ids[(expt,recipient_mouse, t)] for t in (t0, t1)] ]

    donor_data = (donor_freqs, donor_floor**-1.)
    recip_data0 = (recip_freqs_t0, recip_floor_t0**-1.)
    recip_data1 = (recip_freqs_t1, recip_floor_t1**-1.)

    # find barcodes that are valid for inference
    donor_min = 2*donor_floor
    valid_bool = conditioning(donor_freqs, recip_freqs_t0, donor_min=donor_min, recip_max=0, pool_assignment_bool=pool_assignment_bool)
    if valid_bool.sum() < min_barcodes:
        print(f'Not enough barcodes for {recipient_mouse} at {t0} to {t1}')
        return None

    pdetection_inference = migration_inference_by_pdetected(donor_data, recip_data0, recip_data1, valid_bool, recip_floor=0)
    m_detected, m_detected_err, density, errs, bins = pdetection_inference
    ratio_inference = migration_inference_by_ratio(donor_data, recip_data0, recip_data1, valid_bool, max_groups=coarse_grain_num_groups, recip_floor=0)
    m_ratio, (donor_valid, donor_cg), (recip_valid, recip_cg), ratios = ratio_inference
    m_ratio_err = np.abs(np.diff(ratios)[0])/2

    output_dict = {}
    output_dict['freq floors'] = (donor_floor, recip_floor_t0, recip_floor_t1)
    output_dict['detected inference'] = (m_detected * rescale_migration, m_detected_err * rescale_migration)
    output_dict['ratio inference'] = (m_ratio * rescale_migration, m_ratio_err * rescale_migration)
    output_dict['detected raw data'] = (bins, density)
    output_dict['ratio raw data'] = (donor_valid, recip_valid, donor_cg, recip_cg)

    return output_dict

##########################################
#### ENGRAFTMENT TRACKING ####
##########################################

def calc_interpolated_freqs(known_trajectory, remove_zeros=False, Deff_array=[]):
    """ Assumes sampled frequencies never zero! """
    days, freqs = known_trajectory
    days, freqs = np.array(days), np.array(freqs)

    interpolated_days, interpolated_freqs = [], []

    if remove_zeros:
        freqs[freqs == 0] = 1/Deff_array[freqs == 0]

    for day in range(np.min(days), np.max(days)+1):
        if day in days:
            interpolated_days.append(day)
            freq = freqs[days == day][0]
            interpolated_freqs.append(freq)
        else:
            nearest_day_before = np.max(days[days < day]) #nearest sampled day before
            nearest_day_after = np.min(days[days > day])

            freq_before = freqs[days == nearest_day_before][0]
            freq_after = freqs[days == nearest_day_after][0]

            log_slope = ( np.log(freq_after) - np.log(freq_before) ) / (nearest_day_after - nearest_day_before)
            log_interp_freq = np.log(freq_before) + log_slope*(day - nearest_day_before)

            interpolated_days.append(day)
            interpolated_freqs.append( np.exp(log_interp_freq) )
    
    return np.array(interpolated_days), np.array(interpolated_freqs)
      
def trajectory_inference(donor_trajectory, recip_init, donor_transmission, recip_transmission, no_engraftment=False):
    donor_days, donor_freqs = donor_trajectory
    donor_days, donor_freqs = np.array(donor_days), np.array(donor_freqs)
    recip_init = np.array(recip_init)
    recip_transmission  = np.array(recip_transmission)


    assert np.all(np.diff(donor_days) == 1)

    recip_freqs = np.zeros((len(recip_init), len(donor_days)))
    recip_freqs[:, 0] = recip_init

    for i, day in enumerate(donor_days[1:]):
        mean_freq = (donor_freqs[i] + np.sum(recip_freqs[:, i])) / (len(recip_init) + 1)
        effective_fitness = (donor_freqs[i+1] - donor_transmission * mean_freq) / ((1- donor_transmission)*donor_freqs[i])
        if no_engraftment:
            effective_fitness = 0

        # print(effective_fitness)
        recip_freqs[:, i+1] = (1 - recip_transmission) * effective_fitness * recip_freqs[:, i] + recip_transmission * mean_freq

    return donor_days, recip_freqs

def plot_barcode_trajectories(ax, barcode, expt, donor_mice, recip_mice, transmission_dict, color_dict, init_day,
                               plot_other_donors=False, save_seq=False, recip_theory_colors=None):
    barcode_index = np.where(barcodes == barcode)[0][0]

    donor_days_and_rows = {donor_mouse:[*mouse_meta[expt][donor_mouse]] for donor_mouse in donor_mice}
    donor_trajectories = {donor_mouse:[days, freq_array[rows][:, barcode_index], Deff_array[rows]] for donor_mouse, [days, rows, cecum] in donor_days_and_rows.items()}

    donor_interpolated_trajectories = {donor_mouse:[*calc_interpolated_freqs((days, freqs), remove_zeros=True, Deff_array=Deff)] for donor_mouse, (days, freqs, Deff) in donor_trajectories.items()}
    donor_init_freqs = sorted([(donor_mouse,freqs[days == init_day]) for donor_mouse, [days, freqs] in donor_interpolated_trajectories.items()], key = lambda x: x[1], reverse=True)
    reference_donor = donor_init_freqs[0][0] #mouse with highest freq at init day
    donor_transmission = transmission_dict[reference_donor]

    ref_donor_days, reference_freqs = donor_interpolated_trajectories[ reference_donor ]
        
    recip_init = []
    recip_transmission = []
    if len(donor_init_freqs) > 1:
        for donor, init_freq in donor_init_freqs[1:]:
            recip_init.append(init_freq)
            recip_transmission.append(transmission_dict[donor])
    recip_init.extend( [0]*len(recip_mice) ) #assume 0 freq in recipient
    recip_transmission.extend( [transmission_dict[m] for m in recip_mice] )
      
    recip_days, recip_freqs = trajectory_inference((ref_donor_days[ref_donor_days >= init_day], reference_freqs[ref_donor_days >= init_day]), recip_init, donor_transmission, recip_transmission)
    recip_days, recip_freqs_transient = trajectory_inference((ref_donor_days[ref_donor_days >= init_day], reference_freqs[ref_donor_days >= init_day]), recip_init, donor_transmission, recip_transmission, no_engraftment=True)

    #### PLOT ###
    if save_seq != False:
        save_seq.savefig(f'{plot_dir}/E2_traj_presentation_v0.pdf', dpi=300, bbox_inches='tight', transparent=True)
    
    ax.plot(donor_trajectories[reference_donor][0], donor_trajectories[reference_donor][1] + 1e-7, color=color_dict[reference_donor], alpha=1, lw=3)
    if plot_other_donors:
        for (donor, init_freq) in donor_init_freqs[1:]:
            ax.plot(donor_trajectories[donor][0], donor_trajectories[donor][1] + 1e-7, color=color_dict[donor_init_freqs[1][0]], alpha=1, lw=3, zorder=0)

    if recip_theory_colors is None:
        recip_theory_colors = [('grey', 'black')]*len(recip_mice)
    #plot only one recipient theory curves (could do more, but it's a bit crowded)
    for recip_engrafted, recip_transient, colors in zip(recip_freqs[-1:], recip_freqs_transient[-1:], recip_theory_colors[-1:]):
        ax.plot( recip_days[1:], recip_transient[1:], color=colors[0], alpha=1, linestyle='dashed', lw=1.5)
        ax.plot( recip_days[1:], recip_engrafted[1:], color=colors[1], alpha=1, linestyle='dashed', lw=1.5)


    transient_expectation = {k:{day: recip_freqs_transient[1][m] for m, day in enumerate(recip_days)} for k in range(len(recip_mice))}
    # transient_expectation[0] = 0
    # recip_mice = [m for m in [4,5,6,7] if m != donor_mouse]
    reps = []
    for m, recip_mouse in enumerate(recip_mice):
        recip_days, recip_rows = mouse_meta[expt][recip_mouse][:2]
        recip_reads = np.copy(read_array[recip_rows])
        recip_freqs = calc_freqs(recip_reads)
        # recip_freqs[recip_freqs == 0] = 1e-7

        recip_freqs = freq_array[recip_rows][:, barcode_index]
        recip_floor = Deff_array[recip_rows]**-1.

        mask = []
        for r, (day, floor) in enumerate(zip(recip_days, recip_floor)):
            if recip_freqs[r] > 0 and day >= init_day:
                mask.append(True)
            elif day in transient_expectation[m] and transient_expectation[m][day] > floor:
                mask.append(True)
            else:
                mask.append(False)

        ax.plot(recip_days[mask], recip_freqs[mask] + 1e-7, color=color_dict[recip_mouse], marker='.', markersize=5, lw=1.5, linestyle='solid')
        # ax.plot(recip_days, recip_freqs + 1e-7, color=color_dict[recip_mouse], marker='.', markersize=5, lw=1.5, linestyle='solid')

def plot_barcode_trajectories_and_save(ax, barcode, expt, donor_mice, recip_mice, transmission_dict, color_dict, init_day,
                               plot_other_donors=False, save_seq=False, recip_theory_colors=None):
    barcode_index = np.where(barcodes == barcode)[0][0]

    donor_days_and_rows = {donor_mouse:[*mouse_meta[expt][donor_mouse]] for donor_mouse in donor_mice}
    donor_trajectories = {donor_mouse:[days, freq_array[rows][:, barcode_index], Deff_array[rows]] for donor_mouse, [days, rows, cecum] in donor_days_and_rows.items()}

    donor_interpolated_trajectories = {donor_mouse:[*calc_interpolated_freqs((days, freqs), remove_zeros=True, Deff_array=Deff)] for donor_mouse, (days, freqs, Deff) in donor_trajectories.items()}
    donor_init_freqs = sorted([(donor_mouse,freqs[days == init_day]) for donor_mouse, [days, freqs] in donor_interpolated_trajectories.items()], key = lambda x: x[1], reverse=True)
    reference_donor = donor_init_freqs[0][0] #mouse with highest freq at init day
    donor_transmission = transmission_dict[reference_donor]

    ref_donor_days, reference_freqs = donor_interpolated_trajectories[ reference_donor ]
        
    recip_init = []
    recip_transmission = []
    if len(donor_init_freqs) > 1:
        for donor, init_freq in donor_init_freqs[1:]:
            recip_init.append(init_freq)
            recip_transmission.append(transmission_dict[donor])
    recip_init.extend( [0]*len(recip_mice) ) #assume 0 freq in recipient
    recip_transmission.extend( [transmission_dict[m] for m in recip_mice] )
      
    recip_days, recip_freqs = trajectory_inference((ref_donor_days[ref_donor_days >= init_day], reference_freqs[ref_donor_days >= init_day]), recip_init, donor_transmission, recip_transmission)
    recip_days, recip_freqs_transient = trajectory_inference((ref_donor_days[ref_donor_days >= init_day], reference_freqs[ref_donor_days >= init_day]), recip_init, donor_transmission, recip_transmission, no_engraftment=True)
    # print('Recip engrafted:', recip_freqs,)
    # print('Recip transient:', recip_freqs_transient)

    #### PLOT ###
    if save_seq != False:
        save_seq.savefig(f'{plot_dir}/E2_traj_presentation_v0.pdf', dpi=300, bbox_inches='tight', transparent=True)
    
    ax.plot(donor_trajectories[reference_donor][0], donor_trajectories[reference_donor][1] + 1e-7, color=color_dict[reference_donor], alpha=1, lw=3)
    if plot_other_donors:
        for (donor, init_freq) in donor_init_freqs[1:]:
            ax.plot(donor_trajectories[donor][0], donor_trajectories[donor][1] + 1e-7, color=color_dict[donor_init_freqs[1][0]], alpha=1, lw=3, zorder=0)
    if save_seq != False:
        save_seq.savefig(f'{plot_dir}/E2_traj_presentation_v1.pdf', dpi=300, bbox_inches='tight', transparent=True)


    ax.plot( recip_days[1:], recip_freqs_transient[1][1:], color='grey', alpha=1, linestyle='dashed', lw=1.5)
    if save_seq != False:
        save_seq.savefig(f'{plot_dir}/E2_traj_presentation_v2.pdf', dpi=300, bbox_inches='tight', transparent=True)

    ax.plot( recip_days[1:], recip_freqs[1][1:], color='black', alpha=1, linestyle='dashed', lw=1.5)
    if save_seq != False:
        save_seq.savefig(f'{plot_dir}/E2_traj_presentation_v3.pdf', dpi=300, bbox_inches='tight', transparent=True)




    transient_expectation = {k:{day: recip_freqs_transient[1][m] for m, day in enumerate(recip_days)} for k in range(len(recip_mice))}
    # transient_expectation[0] = 0

    # recip_mice = [m for m in [4,5,6,7] if m != donor_mouse]
    reps = []
    for m, recip_mouse in enumerate(recip_mice):
        recip_days, recip_rows = mouse_meta[expt][recip_mouse][:2]
        recip_reads = np.copy(read_array[recip_rows])
        recip_freqs = calc_freqs(recip_reads)
        # recip_freqs[recip_freqs == 0] = 1e-7

        recip_freqs = freq_array[recip_rows][:, barcode_index]
        recip_floor = Deff_array[recip_rows]**-1.

        mask = []
        for r, (day, floor) in enumerate(zip(recip_days, recip_floor)):
            if recip_freqs[r] > 0 and day >= init_day:
                mask.append(True)
            elif day in transient_expectation[m] and transient_expectation[m][day] > floor:
                mask.append(True)
            else:
                mask.append(False)

        ax.plot(recip_days[mask], recip_freqs[mask] + 1e-7, color=color_dict[recip_mouse], marker='.', markersize=5, lw=1.5, linestyle='solid')
        # ax.plot(recip_days, recip_freqs + 1e-7, color=color_dict[recip_mouse], marker='.', markersize=5, lw=1.5, linestyle='solid')
    
def calc_integrated_ratio(donor_freqs, recip_freqs, floor=10**-3):
        ratio = np.ma.masked_array(recip_freqs / donor_freqs, np.full(donor_freqs.shape[0], False))
        ratio[ ratio == 0 ] = floor
        # ratio[ np.isnan(ratio)] = floor
        ratio.mask[ ratio == 0] = True
        ratio.mask[ np.isnan(ratio) ] = True
        ratio.mask[ np.isinf(ratio)] = True
        donor_weights = np.ma.masked_array(donor_freqs, ratio.mask)
        donor_weights /= np.ma.sum(donor_weights)

        return np.ma.sum(ratio * donor_weights )
        # return np.ma.mean(ratio)
        # return np.exp( np.ma.sum( np.log(ratio)*donor_weights ) )