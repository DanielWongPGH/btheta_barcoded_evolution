"""
This script processes barcodes for a single mouse.
"""
import sys
import subprocess

fastq_directory = sys.argv[-2] # directory of fastq files, belonging to single mouse (possibly including day 0)
if fastq_directory[-1] == '/':
    fastq_directory = fastq_directory.rstrip('/')
process = subprocess.Popen(['ls', fastq_directory], stdout=subprocess.PIPE)
fastq_string = str(process.communicate()[0])
fastq_lst = fastq_string.strip("b'").strip("\\n'").split('\\n')
print(fastq_lst)

directory = sys.argv[-1]
if directory[-1] == '/':
    directory = directory.rstrip('/')

days = list()
for fastq in fastq_lst:
    fastq_items = fastq.split('_')
    for item in fastq_items:
        if 'day' in item:
            try:
                day = int(item.strip('day'))
            except:
                day = 100 ## sac set to day 100
            break
    days.append(day)


    extractor_command_lst = ['bartender_extractor_com',
                   '-f', f'{fastq_directory}/{fastq}', # input file
                   '-o', f'{directory}/{day}_extracted', # output file
                   '-p', 'GACT[5]TT[5]AC[5]AA[5]TTGC', # barcode pattern
                   #'-u', '13,8',
                   '-q', '?'] # quality score

    print(extractor_command_lst)

    subprocess.run(extractor_command_lst)

    cluster_command_lst = ['bartender_single_com',
                   '-f', f'{directory}/{day}_extracted_barcode.txt', # input file
                   '-o', f'{directory}/{day}_putative', # output file
                   '-z', '5', # z-score threshold (default = 5)
                   '-d', '2'] # max distance between cluster-able barcodes (default = 2)


    subprocess.run(cluster_command_lst)

days.sort()

time_point_command_string = ''
for day in days:
    time_point_command_string += f'{directory}/{day}_putative_cluster.csv,'
    time_point_command_string += f'{directory}/{day}_putative_quality.csv,'

time_point_command_string.rstrip(',')

time_point_command_lst = ['bartender_combiner_com',
                          '-f', time_point_command_string, #input file
                          '-o', f'{directory}/{day}_longitudinal'] #output file

subprocess.run(time_point_command_lst)

print("Returned {}/longitudinal".format(directory))
