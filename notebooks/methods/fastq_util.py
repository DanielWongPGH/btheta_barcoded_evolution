
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

def pair_R1R2_files(fastq_lst):
    fastq_root_set = set()

    for fastq in fastq_lst:
        fastq_root = fastq.replace('_R1_001.fastq', '').replace('_R2_001.fastq', '').replace('.gz', '')
        if fastq_root not in fastq_root_set:
            R1_file = fastq_root + '_R1_001.fastq'
            R2_file = fastq_root + '_R2_001.fastq'

            if R1_file in fastq_lst and R2_file in fastq_lst:
                fastq_root_set.add(fastq_root)
    return fastq_root_set

rev_comp_map = {'G':'C', 'C':'G', 'A':'T', 'T':'A'}
def reverse_complement(string):
    rc = ''
    for N in string[::-1]:
        rc += rev_comp_map[N]
    return rc

sR = 'aaatgctgttccatcactgg'.upper()
sL = 'tgcagtgtcgaaagaaacaaa'.upper()
sRc = reverse_complement(sR)
sLc = reverse_complement(sL)