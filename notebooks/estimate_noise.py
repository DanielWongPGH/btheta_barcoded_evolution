import numpy as np
import scipy
from scipy.special import gamma, digamma, loggamma, logsumexp
import argparse
import pickle
import matplotlib.pyplot as plt

sR = 'aaatgctgttccatcactgg'.upper()
sL = 'tgcagtgtcgaaagaaacaaa'.upper()
# sRc = reverse_complement(sR)
# sLc = reverse_complement(sL)

## process fastq
def parse_fastq(fastq, cohort):
    if cohort == 'E1':
        fastq_items = fastq.split('_')
        day_str, mouse_str = fastq_items[1], fastq_items[2]

        mouse = int(mouse_str.replace('m', ''))

        if 'day0' in fastq_items[0]:
            day = 0
        elif 'sacsecum' in day_str:
            day = 'cec'
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

            if 'd0' == fastq_items[1]:
                expt = 1134 # 113 & 114
            else:
                expt = int(fastq_items[1])
            barcode_set = int(fastq_items[2].replace('S', ''))

            return (expt, barcode_set, day)

        elif 'mouseexp2' == fastq_items[0]: # NOT DAY 0
            day_str, mouse_str = fastq_items[1], fastq_items[2]

            if 'sac' in day_str:
                if fastq_items[3] == 'c':
                    day = 'cec'
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

def process_fastq(lines_R1, lines_R2):
    bc_umi12_dict = {}
    bc_umi2_dict = {}

    LR_counts = 0
    R_counts = 0

    for i in range(len(lines_R1)//4):
        lineR1 = lines_R1[1+i*4]
        lineR2 = lines_R2[1+i*4]

        try:
            R1_index = lineR1.index(sL)
            R2_index = lineR2.index(sR)

            if R2_index != 8:
                continue

        except:
            continue


        if R1_index == 8:
            umiR1 = lineR1[R1_index-8:R1_index]
        else:
            umiR1 = ''
        umiR2 = lineR2[R2_index-8:R2_index]

        umi = umiR1 + umiR2

        bc_with_repeats = lineR1[R1_index + len(sL) + 4: R1_index + len(sL) + 4 + 34]
        if len(bc_with_repeats) !=  34:
            continue

        # ex. barcodes with repeats
        # GACT AGTCG TT TAACT AC CAGGC AA AGAAG TTGC
        # GACT TCCGA TT CCAGG AC CGGTC AA CTCCG TTGC
        # GACT ACCAG TT GGTTA AC TGGGG AA GTATT TTGC
        # GACT CAGAT TT TTAGC AC CACAC AA ATTCA TTGC

        bc = bc_with_repeats[4:9] + bc_with_repeats[11:16] + bc_with_repeats[18:23] + bc_with_repeats[25:30]
        spacers = bc_with_repeats[0:4] + bc_with_repeats[9:11] + bc_with_repeats[16:18] + bc_with_repeats[23:25] + bc_with_repeats[30:34] 
        if spacers != 'GACTTTACAATTGC':
            continue

        if bc not in bc_umi2_dict:
            bc_umi2_dict[bc] = {}
        if umi not in bc_umi2_dict[bc]:
            bc_umi2_dict[bc][umiR2] = 0
        bc_umi2_dict[bc][umiR2] += 1
        R_counts += 1

        # if R1_index == 8: ## if right UMI # 
        #     if bc not in bc_umi12_dict:
        #         bc_umi12_dict[bc] = {}
        #     if umi not in bc_umi12_dict[bc]:
        #         bc_umi12_dict[bc][umi] = 0
        #     bc_umi12_dict[bc][umi] += 1
        #     LR_counts += 1

    n_paired_reads = (len(lines_R1) // 4)
    print(f'{R_counts//1000}k ({R_counts / n_paired_reads:.2f}) R reads successfully mapped')
    print(f'{LR_counts//1000}k ({LR_counts / n_paired_reads:.2f}) LR reads successfully mapped')
    return bc_umi2_dict, bc_umi12_dict, R_counts, LR_counts, n_paired_reads

def load_bartender_processed_data(pickled_dir):
    with open(f'{pickled_dir}/barcode_read_array.pkl', 'rb') as f:
        read_array = pickle.load(f)

    freq_array = np.einsum('ij, i->ij', read_array, read_array.sum(axis=1)**-1.)

    with open(f'{pickled_dir}/barcodes.pkl', 'rb') as f:
        barcodes = pickle.load(f)

    # with open(f'{pickled_dir}/pickled/mouse_meta.pkl', 'rb') as f:
    #     mouse_meta = pickle.load(f)

    with open(f'{pickled_dir}/vivo_row_ids.pkl', 'rb') as f:
        vivo_row_ids = pickle.load(f)

    # with open(f'{pickled_dir}/pickled/barcode_pool_assignments.pkl', 'rb') as f:
    #     barcode_pool_assignments = pickle.load(f)

    # with open(f'{pickled_dir}/pickled/barcode_pool_map.pkl', 'rb') as f:
    #     index_pool_map = pickle.load(f)

    return read_array, freq_array, barcodes, vivo_row_ids

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

if __name__ == '__main__':
    # fastq file, cohort
    parser = argparse.ArgumentParser()
    parser.add_argument('--fastq')
    parser.add_argument('--cohort')
    parser.add_argument('--pickled_dir')
    parser.add_argument('-d', '--data_dir')
    parser.add_argument('-o', '--output_dir')

    args = parser.parse_args()
    fastq, cohort, fastq_dir, output_dir, pickled_dir = args.fastq, args.cohort, args.data_dir.rstrip('/'), args.output_dir.rstrip('/'), args.pickled_dir.rstrip('/')

    expt, mouse, day = parse_fastq(fastq, cohort)
    if day == "ERROR":
        print("Unparsed:", fastq)

    fout = open(f'{output_dir}/{expt}_M{mouse}_D{day}_noise.out', 'w')


    with open(f'{fastq_dir}/{fastq}_R1_001.fastq') as f:
        lines_R1 = f.readlines()

    with open(f'{fastq_dir}/{fastq}_R2_001.fastq') as f:
        lines_R2 = f.readlines()

    num_reads = len(lines_R1)//4
    if num_reads != len(lines_R2)//4:
        print('Number of forward and reverse reads not equal.')
    print(f'Loaded {fastq} fastq.')


    bc_rUMI_dict, bc_frUMI_dict, R_mapped, LR_mapped, n_paired_reads = process_fastq(lines_R1, lines_R2)
    print(f'Got barcode-UMI pairs from {fastq}.')

    read_array, freq_array, barcodes, vivo_row_ids = load_bartender_processed_data(pickled_dir)
    nominal_depth = read_array[vivo_row_ids[(expt, mouse, day)]].sum()
    barseq_freqs = freq_array[vivo_row_ids[(expt, mouse, day)]]

    data = get_pruned_library_reps(bc_rUMI_dict, barseq_freqs, cutoff=1e-2) 
    mean_above_10, fano_above_10 = np.mean(data[data > 10]), np.var(data[data > 10])/np.mean(data[data > 10])
    mean_below_10, fano_below_10 = np.mean(data[data <= 10]), np.var(data[data <= 10])/np.mean(data[data <= 10])
    p = np.sum(data > 10) / data.shape[0]
    print("data mean and variance:", mean_above_10, fano_above_10, mean_below_10, fano_below_10)

    p = np.max([p, 1e-6])
    init = (p, 10, 1)
    params = run_EM(data, init=init, niter=1000)
    p, Z1, Z2 = params

    # D1 = nominal_depth / (Z1 + Z2 * (1 - p)/p )
    # D2 = (nominal_depth - Z1 * D1) / Z2
    D1 = nominal_depth / (Z1 + Z2 * (1 - p)/p * (1-np.exp(-Z1))/(1-np.exp(-Z2)) )
    D2 = (nominal_depth - Z1 * D1) / Z2
    # D1 = nominal_depth / (Z1/(1-np.exp(-Z1)) + Z2 / (1-np.exp(-Z2)) * (1 - p)/p )
    # D2 = (nominal_depth - Z1/(1-np.exp(-Z1)) * D1) / (Z2 / (1-np.exp(-Z2)))

    D_var = (D1 * Z1 + D2 * Z2)**2 / (Z1*(1+Z1)*D1 + Z2*(1+Z2)*D2)
    D_det = D1 * (1 - np.exp(-Z1)) + D2 * (1 - np.exp(-Z2))

    print("Fit params:", params)

    fout.write(f'{expt}_M{mouse}_D{day}\n')
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


    data_R = get_pruned_library_reps(bc_rUMI_dict, barseq_freqs, cutoff=1) 
    bincount_bc_rUMI = np.bincount(data_R)
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
    ax.hist(data_R, bins=np.logspace(0, 4, 1000), density=True, alpha=0.5, cumulative=True, histtype='step', label='Observed (no cutoff)')
    x = np.logspace(0, 4, 1000)
    ax.plot(x, ZTpois_CDF_mixture(x, params), color='red', lw=1, label='Fitted')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('Reads per barcode-UMI')
    ax.set_ylabel('Num. barcode-UMI')
    ax.set_title(f'{expt}_M{mouse}_D{day}')
    ax.legend()
    fig.savefig(f'{output_dir}/{expt}_M{mouse}_D{day}_noise_plot.pdf')
    plt.close(fig)
    
    

