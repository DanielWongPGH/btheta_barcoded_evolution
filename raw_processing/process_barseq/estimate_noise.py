import numpy as np
import scipy
from scipy.special import gamma, digamma, loggamma, logsumexp
import argparse
import pickle
import gzip
import os
import matplotlib.pyplot as plt

sR1 = 'TGCAGTGTCG' # 10 bases of R1 primer sequence -- UMI is the 8 bases upstream
sR2 = 'AAATGCTGTT' # 10 bases of R2 primer sequence -- UMI is the 8 bases upstream

# library_relabels = {}
# with open('library_mislabelings.tsv', 'r') as f:
#     header = next(f)
#     for line in f:
#         fastq_label_str, true_label_str = line.strip('\n').split('\t')
#         fastq_label_items = fastq_label_str.strip('"').split(',')
#         true_label_str = true_label_str.strip('"').split(',')

#         fastq_label = (fastq_label_items[0], int(fastq_label_items[1]), int(fastq_label_items[2]))
#         if true_label_str[0] == 'None':
#             true_label = None
#         else:
#             true_label = (true_label_str[0], int(true_label_str[1]), int(true_label_str[2]))
#         library_relabels[fastq_label] = true_label


## parse fastq for name details process fastq

def parse_fastq(fastq, cohort):
    try:

        if cohort == 'E1':
            fastq_items = fastq.split('_')
            day_str, mouse_str = fastq_items[1], fastq_items[2]


            if 'day0' in fastq_items[0]:
                day = 0
                mouse = 'P' + fastq_items[1].strip('Set')

            else:
                mouse = int(mouse_str.replace('m', ''))
                if 'sacsecum' in day_str:
                    day = 'cecum'
                elif 'sac' in day_str:
                    if mouse > 10:
                        day = 51
                    else:
                        day = 54
                else:
                    day = int(day_str.replace('day', ''))

            return ('E1', mouse, day)

        if cohort == 'E2':
            fastq_items = fastq.split('_')

            if 'me2' == fastq_items[0]: # DAY 0 
                day = 0
                barcode_set = fastq_items[2]

                if '116' in fastq: barcode_set += '_delay'  # delay

                return ('E2', barcode_set, day)

            elif 'mouseexp2' == fastq_items[0]: # NOT DAY 0
                day_str, mouse_str = fastq_items[1], fastq_items[2]

                if 'sac' in day_str:
                    if fastq_items[3] == 'c':
                        day = 'cecum'
                    else:
                        day = 55
                else:
                    day = int(day_str.replace('day', ''))

                expt_str, mouse = mouse_str.split('m')
                expt, mouse = int(expt_str), int(mouse)

                if expt == 114:
                    mouse += 5
                elif expt == 116:
                    mouse += 10

                return ('E2', mouse, day)
            else:
                return (fastq, 0, 'ERROR')
            
        if cohort == 'vitro':
            fastq_items = fastq.split('_')
            if 'd0' in fastq:
                day = 0
                pool = 'V' + fastq_items[3].strip('S')

                return ('vitro', pool, day)
            else:
                if 'passage' in fastq:
                    passage = int(fastq_items[1].strip('passage'))
                    well = fastq_items[2]+ '_' + fastq_items[3]
                    return ('vitro', well, passage)
                elif 'day' in fastq:
                    day = int(fastq_items[1].strip('day'))
                    well = 'p'+fastq_items[0].strip('Plate') + '_' + fastq_items[2]
                    return ('vitro', well, day)

    except: return None

def hamming(str1, str2):
    return sum([str1[i] != str2[i] for i in range(len(str1))])

def process_fastq(lines_R1, lines_R2): 
    ## we use these to roughly estimate the number of templates
    ## we are strict about the UMI sequencences and positions, so miss some genuine reads
    ## but as long as seq errors are rare, this should be a small fraction of the total reads 
    bc_umi12_dict = {}
    bc_umi2_dict = {}

    LR_counts = 0
    R_counts = 0

    for i in range(len(lines_R1)//4):
        lineR1 = lines_R1[1+i*4]
        lineR2 = lines_R2[1+i*4]

        try:
            R1_index = lineR1.index(sR1)
            R2_index = lineR2.index(sR2)
            if R2_index != 8: #do not tolerate error in position of right UMI
                continue
            umiR2 = lineR2[R2_index-8:R2_index]
        except ValueError:
            continue

        bc_with_repeats = lineR1[R1_index + len(sR1) + 15: R1_index + len(sR1) + 15 + 34]
        if len(bc_with_repeats) !=  34: continue
        bc = bc_with_repeats[4:9] + bc_with_repeats[11:16] + bc_with_repeats[18:23] + bc_with_repeats[25:30]
        spacers = bc_with_repeats[0:4] + bc_with_repeats[9:11] + bc_with_repeats[16:18] + bc_with_repeats[23:25] + bc_with_repeats[30:34] 
        if spacers != 'GACTTTACAATTGC': continue #do not tolerate any errors in spacers

        if bc not in bc_umi2_dict:
            bc_umi2_dict[bc] = {}
        if umiR2 not in bc_umi2_dict[bc]:
            bc_umi2_dict[bc][umiR2] = 0
        bc_umi2_dict[bc][umiR2] += 1
        R_counts += 1

        if R1_index == 8: ## if right UMI # 
            umiR1 = lineR1[R1_index-8:R1_index]
            umi = umiR1 + umiR2
            if bc not in bc_umi12_dict:
                bc_umi12_dict[bc] = {}
            if umi not in bc_umi12_dict[bc]:
                bc_umi12_dict[bc][umi] = 0
            bc_umi12_dict[bc][umi] += 1
            LR_counts += 1

    n_paired_reads = (len(lines_R1) // 4)
    print(f'{R_counts//1000}k ({R_counts / n_paired_reads:.2f}) R reads successfully mapped')
    print(f'{LR_counts//1000}k ({LR_counts / n_paired_reads:.2f}) LR reads successfully mapped')
    return bc_umi2_dict, bc_umi12_dict, R_counts, LR_counts, n_paired_reads

def load_bartender_processed_data(bartender_dir, fastq):
    barcodes = []
    read_array = []
    with open(f'{bartender_dir.rstrip("/")}/{fastq}_cluster.csv', 'r') as f:
        header = next(f) #cluster idx, cluster barcode, cluster score, counts        
        for line in f:
            line_items = line.strip('\n').split(',')
            barcode = line_items[1]
            count = int(line_items[3])
            barcodes.append(barcode)
            read_array.append(count)

    read_array = np.array(read_array)
    barcodes = np.array(barcodes)
    freq_array = read_array / read_array.sum()

    return read_array, freq_array, barcodes

# get pruned library
def get_pruned_library_reps(fastq_parsed_library, barseq_frequencies, cutoff=1e-2, min_cutoff=1e-5):
    # check consistency of ignoring error correction
    # consider barcodes with total frequency < cutoff 
    pruned_library_counts = []
    valid_barcodes = barcodes[barseq_frequencies < cutoff] 
    for barcode in valid_barcodes:
        try: 
            pruned_library_counts.extend( list(fastq_parsed_library[barcode].values()) )  
        except:
            continue 
    pruned_library_counts = np.array(pruned_library_counts)    
    return pruned_library_counts        

## estimate noise by fitting 2-component ZTPoisson mixture model

def ZTpois_negLL_deriv(params, data, weights):
    assert data.shape == weights.shape
    mean = params

    dll_dmean = - np.sum(weights * (data/mean - 1 - 1/(np.exp(mean) - 1)))
    return np.array(dll_dmean)

def ZTpois_negLL(params, data, weights):
    assert data.shape == weights.shape
    mean = params

    return - np.sum( weights * (data*np.log(mean) - mean - np.log1p(-np.exp(-mean)) ))

def ZTpois_class_probs(params, data):
    p_class1, mean1, mean2 = params
    p_class2 = 1-p_class1


    LL1 = data*np.log(mean1) - mean1 - np.log1p(-np.exp(-mean1)) 
    LL2 = data*np.log(mean2) - mean2 - np.log1p(-np.exp(-mean2)) 
    # log_ratio = LL2 - LL1

    # print(np.exp(LL1)/np.exp(LL2))
    # logprob1_v =  np.log( p_class1 * np.exp(LL1) / (p_class1 * np.exp(LL1) + p_class2 * np.exp(LL2)) )
    logprob1 = np.log(p_class1) + LL1 - logsumexp([np.log(p_class1) + LL1, np.log(p_class2) + LL2], axis=0)

    # print(logprob1_v[:100], logprob1[:100])
    return logprob1

def ZTpois_CDF_mixture(x, params):
    class1_prob, mean1, mean2 = params

    cdf1 = (scipy.stats.poisson.cdf(x, mean1)-scipy.stats.poisson.pmf(0, mean1))/(1-scipy.stats.poisson.pmf(0, mean1)) 
    cdf2 = (scipy.stats.poisson.cdf(x, mean2)-scipy.stats.poisson.pmf(0, mean2))/(1-scipy.stats.poisson.pmf(0, mean2))

    return class1_prob * cdf1 + (1-class1_prob) * cdf2
    
def EM_step(params, data):
    p_class1, mean1, mean2  = params

    # get probabilities of each data point belonging to each component
    logweights_class1 = ZTpois_class_probs(params, data)
    logweights_class2 = np.log1p(-np.exp(logweights_class1))

    p_class1_new = np.exp(logsumexp(logweights_class1)) / data.shape[0]
    if p_class1_new >= 1:
        p_class1_new = 1-1e-6
    elif p_class1_new <= 0:
        p_class1_new = 1e-6
    
    res_class1 = scipy.optimize.minimize(ZTpois_negLL, x0=[mean1], args=(data, np.exp(logweights_class1)), bounds=[(1e-5, 10**5)], tol=1e-6, jac=ZTpois_negLL_deriv)
    res_class2 = scipy.optimize.minimize(ZTpois_negLL, x0=[mean2], args=(data, np.exp(logweights_class2)), bounds=[(1e-5, 10**5)], tol=1e-6, jac=ZTpois_negLL_deriv)

    mean1_new = res_class1.x[0]
    mean2_new = res_class2.x[0]

    print(p_class1_new, mean1_new, mean2_new)

    return np.array([p_class1_new, mean1_new, mean2_new])

def run_EM(data, init=(0.5, 1, 0.9, 50, 0.9), niter=100, err=1e-6):
    params = init
    for i in range(niter):
        old_params = np.copy(params)
        params = EM_step(params, data)

        if np.all(np.abs(old_params - params)/params < err):
            print('Converged after', i, 'iterations')
            break
    return params

def open_fastq(path):
    if os.path.exists(path):
        return open(path)
    if os.path.exists(path + '.gz'):
        return gzip.open(path + '.gz', 'rt')
    raise FileNotFoundError(path)


if __name__ == '__main__':
    # fastq file, cohort
    parser = argparse.ArgumentParser()
    parser.add_argument('--fastq')
    parser.add_argument('--bartender_dir')
    parser.add_argument('-f', '--fastq_dir')
    parser.add_argument('-o', '--output_dir')
    # parser.add_argument('--lib_fastq_map')

    args = parser.parse_args()
    fastq, fastq_dir, output_dir, bartender_dir = args.fastq, args.fastq_dir.rstrip('/'), args.output_dir.rstrip('/'), args.bartender_dir.rstrip('/')

    print(f'Processing {fastq}.')

    with open_fastq(f'{fastq_dir}/{fastq}_R1_001.fastq') as f:
        lines_R1 = f.readlines()
    with open_fastq(f'{fastq_dir}/{fastq}_R2_001.fastq') as f:
        lines_R2 = f.readlines()

    num_reads = len(lines_R1)//4
    if num_reads != len(lines_R2)//4:
        print('Number of forward and reverse reads not equal.')
    print(f'Loaded {fastq} fastq.')

    bc_rUMI_dict, bc_frUMI_dict, R_mapped, LR_mapped, n_paired_reads = process_fastq(lines_R1, lines_R2)
    print(f'Got barcode-UMI pairs from {fastq}.')
    barseq_reads, barseq_freqs, barcodes = load_bartender_processed_data(bartender_dir, fastq)
    nominal_depth = barseq_reads.sum()


    data = get_pruned_library_reps(bc_rUMI_dict, barseq_freqs, cutoff=1e-2) 
    mean_above_10, fano_above_10 = np.mean(data[data > 10]), np.var(data[data > 10])/np.mean(data[data > 10])
    mean_below_10, fano_below_10 = np.mean(data[data <= 10]), np.var(data[data <= 10])/np.mean(data[data <= 10])
    p = np.sum(data > 10) / data.shape[0]
    print("data mean and variance:", mean_above_10, fano_above_10, mean_below_10, fano_below_10)

    p = np.max([p, 1e-6])
    init = (p, 10, 1)
    params = run_EM(data, init=init, niter=1000)
    p, Z1, Z2 = params

    D1 = nominal_depth / (Z1 + Z2 * (1 - p)/p * (1-np.exp(-Z1))/(1-np.exp(-Z2)) )
    D2 = (nominal_depth - Z1 * D1) / Z2

    D_var = (D1 * Z1 + D2 * Z2)**2 / (Z1*(1+Z1)*D1 + Z2*(1+Z2)*D2)
    D_det = D1 * (1 - np.exp(-Z1)) + D2 * (1 - np.exp(-Z2))

    print("Fit params:", params)

    fout = open(f'{output_dir}/{fastq}_noise.out', 'w')

    fout.write(fastq+'\n')
    fout.write(f'Nominal Depth\t{nominal_depth}\n')
    fout.write(f'Dvar\t{D_var}\n')
    fout.write(f'Ddet\t{D_det}\n')
    fout.write(f'p\t{p}\n')
    fout.write(f'Z1\t{Z1}\n')
    fout.write(f'Z2\t{Z2}\n')
    fout.write(f'D1\t{D1}\n')
    fout.write(f'D2\t{D2}\n')
    fout.write(f'\n')
    fout.write(f'R-UMI mapped\t{R_mapped}\n')
    fout.write(f'LR-UMI mapped\t{LR_mapped}\n')
    fout.write(f'paired reads\t{n_paired_reads}\n')


    data_all = get_pruned_library_reps(bc_rUMI_dict, barseq_freqs, cutoff=1) 
    bincount_bc_rUMI = np.bincount(data_all)
    string_bc_rUMI = ','.join([str(x) for x in bincount_bc_rUMI])
    fout.write(f'n. barcode-rUMI reads\t{string_bc_rUMI}\n')

    # data_Rcutoff = get_pruned_library_reps(bc_rUMI_dict, barseq_freqs, cutoff=1e-2)
    bincount_bc_rUMI_cutoff = np.bincount(data)
    string_bc_rUMI_cutoff = ','.join([str(x) for x in bincount_bc_rUMI_cutoff])
    fout.write(f'n. barcode-rUMI reads (<1%)\t{string_bc_rUMI_cutoff}\n')

    # data_LR = get_pruned_library_reps(bc_frUMI_dict, barseq_freqs, cutoff=1)
    # bincount_bc_frUMI = np.bincount(data_LR)
    # string_bc_frUMI = ','.join([str(x) for x in bincount_bc_frUMI])
    # fout.write(f'n. barcode-frUMI reads\t{string_bc_frUMI}\n')
    fout.close()


        # # params[2] /= 2
    fig, ax = plt.subplots()
    ax.hist(data, bins=np.logspace(0, 4, 1000), density=True, alpha=0.5, cumulative=True, histtype='step', label='Observed')
    ax.hist(data_all, bins=np.logspace(0, 4, 1000), density=True, alpha=0.5, cumulative=True, histtype='step', label='Observed (no cutoff)')
    x = np.logspace(0, 4, 1000)
    ax.plot(x, ZTpois_CDF_mixture(x, params), color='red', lw=1, label='Fitted')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('Reads per barcode-UMI')
    ax.set_ylabel('Num. barcode-UMI')


    ax.set_title(fastq)
    ax.legend()
    fig.savefig(f'{output_dir}/{fastq}_noise_plot.pdf')
    plt.close(fig)
    
    
