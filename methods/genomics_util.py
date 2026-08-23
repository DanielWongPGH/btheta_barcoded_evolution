import pickle
import numpy as np
from methods.config import *
import Bio.SeqIO as SeqIO
from Bio.Seq import Seq
import numpy as np


with open(f'{pickled_dir}/PUL_to_gene_map.pkl', 'rb') as f:
    PUL_gene_map = pickle.load(f)

with open(f'{pickled_dir}/gene_to_PUL_map.pkl', 'rb') as f:
    gene_PUL_map = pickle.load(f)


genome = Seq(SeqIO.read(f'{data_dir}/reference_genome/BtVPI.fasta', 'fasta').seq)
GENOME_LENGTH = len(genome)

gene_coords_map = {}
gene_description = {}
with open(f'{data_dir}/reference_genome/BtVPI.ptt', 'r') as f:
    for _ in range(4):
        next(f)
    for line in f:
        line_items = line.rstrip('\n').split('\t')
        start, stop = [int(value) for value in line_items[0].split('..')]
        gene_coords_map[line_items[5]] = (start, stop, line_items[1])
        gene_description[line_items[5]] = line_items[-1]

gene_array = np.array(sorted([
    [gene, start, stop] for gene, (start, stop, strand) in gene_coords_map.items()
]))
def find_gene_from_position(pos):
    gene_starts = gene_array[:,1].astype(int)
    gene_stops = gene_array[:,2].astype(int)
    gene_names = gene_array[:,0]

    start_bool = pos >= gene_starts
    stop_bool = pos <= gene_stops

    if (start_bool * stop_bool).sum() > 0:
        return gene_names[start_bool * stop_bool][0]
    
    else:
        upstream_gene = gene_names[start_bool][-1]
        downstream_gene = gene_names[stop_bool][0]
        return f'{upstream_gene}-{downstream_gene}'
    
def get_variants_from_annotated_gd_file(gd_file, freq=False, track_RA=False):
    clone_evidence = {'INS':[], 'SNP':[], 'DEL':[], 'SUB':[], 'JC':[], 'MC':[], 'RA':[]}

    prev_RA = (-1, '', 0, 0) #dummy
    with open(gd_file, 'r') as f:
        for line in f:
            line_items = line.strip('\n').split('\t')
            variant_type = line_items[0]

            if variant_type == '#=MAPPED-BASES':
                coverage = int(line_items[1]) / GENOME_LENGTH
            if variant_type not in ['INS', 'SNP', 'DEL', 'SUB', 'JC', 'MC', 'RA']:
                continue

            if variant_type in ['INS', 'SNP', 'DEL', 'SUB']:
                try:
                    position = int(line_items[4])
                    if variant_type == 'SUB':
                        mutation = line_items[6]
                        meta_tup_list = [item.split("=") for item in line_items[7:]]
                    else:
                        mutation = line_items[5]
                        meta_tup_list = [item.split("=") for item in line_items[6:]]
                        if variant_type == 'DEL':
                            mutation = '-'+mutation
                        
                    mutation_dict = {item[0]:item[1] for item in meta_tup_list}
                    ref_seq = mutation_dict['ref_seq']
                    if variant_type == 'INS':
                        mutation = f'+{mutation}'
                    elif variant_type == 'DEL':
                        mutation = f'-{ref_seq}'
                    else:
                        mutation = f'{ref_seq}->{mutation}'
                    gene_description = mutation_dict['gene_product']


                    gene = mutation_dict['locus_tag'].replace('_','').replace('[','').replace(']','').replace('|','/') #formatting
                    if 'gene_strand' in mutation_dict:
                        strand = mutation_dict['gene_strand'].replace('>', '+').replace('<', '-').replace('|','/')
                    else:
                        strand = ''
                    mutation_category = mutation_dict['mutation_category']

                    if 'snp' in mutation_category and mutation_category.replace('snp_','') not in ['intergenic', 'noncoding']:
                        codon_ref, aa_ref = mutation_dict['codon_ref_seq'], mutation_dict['aa_ref_seq']
                        codon_new, aa_new = mutation_dict['codon_new_seq'], mutation_dict['aa_new_seq']
                        aa_position = mutation_dict['aa_position']
                        mutation_description = f'{aa_ref}{aa_position}{aa_new} ({codon_ref}->{codon_new})'
                    else:
                        if 'gene_position' in mutation_dict:
                            mutation_description = f'{mutation} ({mutation_dict["gene_position"]})'
                        else:
                            mutation_description = ''

                    if mutation_category == 'snp_synonymous':
                        syn_flag = 'True'
                    elif 'intergenic' in mutation_description:
                        syn_flag = ''
                    else:
                        syn_flag = 'False'

                    if '/' in gene:
                        genes = gene.split('/')
                        if genes[0] in gene_PUL_map:
                            PUL1 = gene_PUL_map[genes[0]][0]
                        else:
                            PUL1 = ''
                        if genes[1] in gene_PUL_map:
                            PUL2 = gene_PUL_map[genes[1]][0]
                        else:
                            PUL2 = ''
                        if PUL1 == PUL2:
                            PUL = PUL1
                        else:
                            PUL = f'{PUL1}-{PUL2}'
                    else:
                        if gene in gene_PUL_map:
                            PUL = gene_PUL_map[gene][0]
                        else:
                            PUL = ''

                    pruned_mutation_dict = {'position':position, 'syn':syn_flag,
                                            'gene':gene, 'strand':strand, 'gene_description':gene_description, 'PUL':PUL, 'mutation':mutation,
                                            'mutation_category':variant_type, 'mutation_description':mutation_description}

                    if freq:
                        frequency = float(mutation_dict['frequency'])
                        pruned_mutation_dict['frequency'] = frequency

                    clone_evidence[variant_type].append( (position, variant_type, pruned_mutation_dict) )
                except:
                    pass
                    # print(gd_file, line)

            if variant_type == 'MC':
                start, end = int(line_items[4]), int(line_items[5])
                gene1 = find_gene_from_position(start)  # first gene
                gene2 = find_gene_from_position(end)  # second gene

                clone_evidence['MC'].append( ((start, end), (gene1, gene2)) )

            if variant_type == 'JC':
                start, end = int(line_items[4]), int(line_items[7])
                gene1 = find_gene_from_position(start) # first gene
                gene2 = find_gene_from_position(end) # second gene
                orientations = (int(line_items[5]), int(line_items[8]))

                evidence_meta = {e.split('=')[0]:e.split('=')[1] for e in line_items[11:]}
                JC_coverage = int(evidence_meta['coverage_minus']) + int(evidence_meta['coverage_plus'])

                try:
                    freq = float(evidence_meta['frequency'])
                except:
                    # print(line_items)
                    freq = -1

                clone_evidence['JC'].append( ((start, end), (gene1, gene2), orientations, freq, JC_coverage) )

            if variant_type == 'RA' and track_RA:
                position = int(line_items[4])
                evidence_index = int(line_items[5])
                gene = find_gene_from_position(position) # first gene
                evidence_meta = {e.split('=')[0]:e.split('=')[1] for e in line_items[8:]}
                freq = float(evidence_meta['polymorphism_frequency'])
                RA_coverage = np.sum(int(e) for e in evidence_meta['new_cov'].split('/'))
                tot_coverage = RA_coverage/1
                if tot_coverage < 5:
                    continue
                # print(position, evidence_meta)

                ref_seq, new_seq = line_items[6], line_items[7]
                if new_seq == '.': #deletion
                    mutation_category = 'DEL'
                    mutation = f'-{ref_seq}'
                elif ref_seq == '.': #insertion
                    mutation_category = 'INS'
                    mutation = f'+{new_seq}'
                else: 
                    mutation_category = 'SNP'
                    mutation = f'{ref_seq}->{new_seq}'

                # if variant_type == 'DEL' and prev_RA[0] == position-1 and prev_RA[1] == 'DEL':
                #     # TODO: merge
                # elif variant_type == 'INS' and prev_RA[0] == position-1 and prev_RA[1] == 'INS':
                #     # TODO: merge w/ previous
                # elif variant_type == 'SNP' and prev_RA[0] == position-1 and prev_RA[1] == 'SUB':
                #     # TODO: merge w/ previous
                
                syn_flag = ''
                if 'gene_product' not in evidence_meta:
                    gene = find_gene_from_position(position)
                    if '/' in gene:
                        gene1, gene2 = gene.split('/')
                        gene1_description = gene_description[gene1]
                        gene2_description = gene_description[gene2]

                        gene1_strand = gene_coords_map[gene1][2]
                        gene2_strand = gene_coords_map[gene2][2]

                        gene_description = '/'.join([gene1_description, gene2_description])
                        strand = '/'.join([gene1_strand, gene2_strand])
                    mutation_description  = mutation
                else:
                    gene_description = evidence_meta['gene_product']
                    gene = evidence_meta['locus_tag'].replace('_','').replace('[','').replace(']','').replace('|','/')
                    strand = evidence_meta['gene_strand'].replace('>', '+').replace('<', '-').replace('|','/')
                    gene_position = evidence_meta['gene_position']
                    if 'intergenic' in gene_position or 'coding' in gene_position or 'noncoding' in gene_position:
                        mutation_description = f'{mutation} ({gene_position})'
                    else:
                        ref_codon = evidence_meta['codon_ref_seq']
                        new_codon = evidence_meta['codon_new_seq']
                        ref_aa = evidence_meta['aa_ref_seq']
                        new_aa = evidence_meta['aa_new_seq']
                        
                        if new_aa == ref_aa: syn_flag = 'True'
                        else: syn_flag = 'False'
                    
                if '/' in gene:
                    gene1, gene2 = gene.split('/')
                    if gene1 in gene_PUL_map: PUL1 = gene_PUL_map[gene1][0]
                    else: PUL1 = ''
                    if gene2 in gene_PUL_map: PUL2 = gene_PUL_map[gene2][0]
                    else: PUL2 = ''
                    if PUL1 == PUL2: PUL = PUL1
                    else: PUL = f'{PUL1}-{PUL2}'
                else:
                    if gene in gene_PUL_map: PUL = gene_PUL_map[gene][0]
                    else: PUL = ''

                    
                pruned_mutation_dict = {'position':position, 'syn':syn_flag, 'frequency':frequency,
                                       'gene':gene, 'strand':strand, 'gene_description':gene_description, 'PUL':PUL, 'mutation':mutation,
                                            'mutation_category':mutation_category, 
                                            'mutation_description':mutation_description, 'coverage':RA_coverage}


                clone_evidence['RA'].append( (position, mutation_category, pruned_mutation_dict) )
    return coverage, clone_evidence






def parse_isolate_gd_filename(filename):
    file_meta = filename.strip('.gd').split('_')
    try:
        mouse, day, clone = int(file_meta[0].strip('m')), int(file_meta[1].strip('day')), int(file_meta[2].strip('clone'))
    except:
        print(f'Unable to parse file: {filename}')

    return (mouse, day, clone)


def check_preexisting(mut, preexisting_variants, exclude_genes={}):
    if mut[2] in exclude_genes: # e.g. hypermutable
        return True
    elif mut in preexisting_variants:
        return True
    return False

# def check_preexisting(mut, preexisting_variants, min_day0_reps=5, exclude_genes={}):
#     if mut[2] in exclude_genes: # e.g. hypermutable
#         return True
#     elif mut not in preexisting_variants:
#         return False
#     elif len(preexisting_variants[mut]['day0']) >= min_day0_reps or len(preexisting_variants[mut]['isolate'])>0:
#         return True
#     return False

### JC-type evidence

def process_sample_JC_evidence(sample_label, sample_evidence, JC_dict, day0_JC_coords={}, min_coverage=5, min_size=10, min_freq=0.2):
    for mut in sample_evidence:
        (start, end), genes, orientations, freq, JC_coverage = mut
        if check_JC_against_list((start, end), day0_JC_coords) \
                or JC_coverage <= min_coverage \
                or (end-start) <= min_size \
                or freq <= min_freq:
            continue
        if (start, end, orientations[0], orientations[1]) not in JC_dict:
            JC_dict[(start, end, orientations[0], orientations[1])] = []
        JC_dict[(start, end, orientations[0], orientations[1])].append( (sample_label, freq, JC_coverage) )

def check_overlap(MC1, MC2):
    mc1, mc2 = sorted([MC1, MC2], key=lambda x:x[0]) # ( ) ( )

    if mc2[0] >= mc1[1]:
        return False, [MC1, MC2]
    elif mc1[1] < mc2[1]:
        return True, [(mc1[0], mc2[1]), (mc1[0], mc2[1])] #merged
    else:
        return True, [mc1, mc1]

def check_JC_overlap(JC1, JC2, delta=2):
    if np.abs(JC1[0] - JC2[0]) < delta and np.abs(JC1[1] - JC2[1]) < delta:
        return True
    else:
        return False

def check_JC_against_list(JC, lst_of_JCs):
    for check_JC in lst_of_JCs:
        if check_JC_overlap(JC, check_JC):
            return True
    return False

###

inverted_repeats = {}
with open(f'{data_dir}/reference_genome/BtVPI_plausible_inverted_repeats.txt', 'r') as f:
    start1, stop1, start2, stop2, ir1, ir2 = None, None, None, None, None, None
    for line in f:
        line_items = line.strip().split()
        if line_items == []:
            if ir1 is not None:
                inverted_repeats[((start1, stop1), (start2, stop2))] = (ir1, ir2)
            count = 0
        if count == 2:
            start1, ir1, stop1 = int(line_items[0]), line_items[1].upper(), int(line_items[-1]) 
        if count == 4:
            start2, ir2, stop2 = int(line_items[-1]), line_items[1][::-1].upper(), int(line_items[0])

        count += 1



def find_spanning_inverted_repeats(JCcoords1, JCcoords2, delta=1):
    ''' JC = (start, stop) '''
    repeat_idx_array = np.array([[ir1[0], ir1[1]] for ir1, ir2 in inverted_repeats.keys()])
    inv_repeat_idx_array = np.array([[ir2[0], ir2[1]] for ir1, ir2 in inverted_repeats.keys()])

    start1, stop1 = JCcoords1
    start2, stop2 = JCcoords2

    start1_delta_repeat = start1 - repeat_idx_array 
    stop1_delta_repeat = stop1 - inv_repeat_idx_array

    start1_bool = (start1_delta_repeat[:, 0] > -delta) * (start1_delta_repeat[:, 1] < delta)
    stop1_bool = (stop1_delta_repeat[:, 0] > -delta) * (stop1_delta_repeat[:, 1] < delta)

    start2_delta_inv_repeat = start2 - repeat_idx_array
    stop2_delta_inv_repeat = stop2 - inv_repeat_idx_array

    start2_bool = (start2_delta_inv_repeat[:, 0] > -delta) * (start2_delta_inv_repeat[:, 1] < delta)
    stop2_bool = (stop2_delta_inv_repeat[:, 0] > -delta) * (stop2_delta_inv_repeat[:, 1] < delta)

    

    flag = (start1_bool * stop1_bool * start2_bool * stop2_bool).sum()
    return flag




### genome annotations
gene_to_PUL_map = {}
PUL_to_gene_map = {}
with open(f'{data_dir}/reference_genome/BtVPI_PULs_Martens2008.tsv', 'r') as f:
    for _ in range(4):
        header = next(f)

    PUL_num = None
    PUL_genes = []
    for line in f:
        line_items = line.strip('\n').split('\t')
        if line_items[0]: # new PUL!
            if PUL_num:
                PUL_to_gene_map[PUL_num] = PUL_genes
                for gene in PUL_genes:
                    gene_to_PUL_map[gene] = [PUL_num, [other_gene for other_gene in PUL_genes if other_gene != gene]]

            PUL_num = int(line_items[0])
            gene = line_items[2]
            PUL_genes = [gene]
        else:
            gene = line_items[2]

            if gene:
                PUL_genes.append(gene)

    # last PUL
    PUL_to_gene_map[PUL_num] = PUL_genes
    for gene in PUL_genes:
        gene_to_PUL_map[gene] = [PUL_num, [other_gene for other_gene in PUL_genes if other_gene != gene]]


CPS_positions = {'cps1': (462368, 489867, 1),
 'cps2': (573167, 595424, 1),
 'cps3': (734869, 761399, 1),
 'cps4': (1663708, 1685198, -1),
 'cps5': (2032185, 2047991, -1),
 'cps6': (2104512, 2124231, -1),
 'cps7': (3580840, 3605962, -1),
 'cps8': (32315, 71500, 1)}

gene_CPS_map = {f'BT0{x}':'cps1' for x in range(375, 401)} \
            | {f'BT0{x}':'cps2' for x in range(462, 483)} \
            | {f'BT0{x}':'cps3' for x in range(595, 616)} \
            | {f'BT{x}':'cps4' for x in range(1338, 1359)} \
            | {f'BT{x}':'cps5' for x in range(1642, 1657)} \
            | {f'BT{x}':'cps6' for x in range(1707, 1726)} \
            | {f'BT{x}':'cps7' for x in range(2862, 2887)} \
            | {f'BT00{x}':'cps8' for x in range(36, 75)}

CPS_gene_map = {'cps1': [f'BT0{x}' for x in range(375, 401)]} \
            | {'cps2': [f'BT0{x}' for x in range(462, 483)]} \
            | {'cps3': [f'BT0{x}' for x in range(595, 616)]} \
            | {'cps4': [f'BT{x}' for x in range(1338, 1359)]} \
            | {'cps5': [f'BT{x}' for x in range(1642, 1657)]} \
            | {'cps6': [f'BT{x}' for x in range(1707, 1726)]} \
            | {'cps7': [f'BT{x}' for x in range(2862, 2887)]} \
            | {'cps8': [f'BT00{x}' for x in range(36, 75)]}
