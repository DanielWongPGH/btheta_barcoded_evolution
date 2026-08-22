"""Generate population_samples.csv from the metagenomics FASTQ manifests.

The historical filename parser lives in the archive-deposition utilities and
is supplied explicitly with --parser-root rather than assumed to be at a
machine-specific location.
"""
import argparse
from pathlib import Path
import sys


SCRIPT_DIR = Path(__file__).resolve().parent

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--e1-manifest", type=Path,
                    default=SCRIPT_DIR / "E1_MGX_fastq_stems.txt")
parser.add_argument("--e2-manifest", type=Path,
                    default=SCRIPT_DIR / "E2_MGX_fastq_stems.txt")
parser.add_argument("--output", type=Path,
                    default=SCRIPT_DIR / "population_samples.csv")
parser.add_argument(
    "--parser-root", type=Path, required=True,
    help="directory containing prep_fastq_for_SRA/parse_fastqs_for_SRA.py",
)
args = parser.parse_args()

sys.path.insert(0, str(args.parser_root.resolve()))
import prep_fastq_for_SRA.parse_fastqs_for_SRA as pf


def read_manifest(path):
    with path.open() as handle:
        fastq_dir = next(handle).strip().rstrip("/")
        stems = [line.strip() for line in handle if line.strip()]
    if not fastq_dir:
        raise ValueError("Manifest has no FASTQ directory header: {}".format(path))
    return fastq_dir, stems


E1_dir, E1_stems = read_manifest(args.e1_manifest)
E2_dir, E2_stems = read_manifest(args.e2_manifest)

fastq_set = set()
for expt, stems in [('E1', E1_stems), ('E2', E2_stems)]:
    for stem in stems:
        fastq_set.add((expt, stem))


all_samples = {}
for expt, fastq in fastq_set:
    dtype = 'metagenomics'
    (expt, mouse, day, site, index, rep) = pf.parse_mgx_file(expt, fastq)

    key = (expt, mouse, day, site)

    if key not in all_samples:
        all_samples[key] = []
    all_samples[key].append( (fastq, index, dtype, f'rep{rep}') )
    


# ##  get well to invitro set 
# well_to_bc_set = {}
# with open('../invitro_well_identities.tsv') as f:
#     header = next(f) #plate \t well \t carbon source \t barcode inoc \t notes
#     for line in f:
#         line_items = line.split('\t')
#         plate = line_items[0]
#         well_num = line_items[1]

#         well = f'p{plate}_{well_num}'
#         bc = line_items[3]

#         well_to_bc_set[well] = bc


## string longitudinal samples together, including day 0
timecourse_map = {}
for key, files in all_samples.items():
    (expt, mouse, day, site) = key
    if day == 0:
        continue
    if site == 'cecum':
        day = 100
    if site == 'smallint':
        day = 200

    host = (expt, mouse)
    if host not in timecourse_map:
        timecourse_map[host] = {}

    # get file that is rep1
    rep1_files = [f for f in files if 'rep1' in f]
    timecourse_map[host][day] = rep1_files[0][0]

print(timecourse_map.keys())

## Add in day 0 for each replicate
for sample, fastqs in all_samples.items():
    if sample[2] == 0:
        print(sample, fastqs)

for host in list(timecourse_map.keys()):
    # case by case
    (expt, mouse) = host
    if expt == 'E1':
        if mouse < 6: day0_sample = (expt, 'P2', 0, 'vitro')
        else: day0_sample = (expt, 'P1', 0, 'vitro')

    elif expt == 'E2':
        if mouse > 15: #skip isolator 118
            del timecourse_map[host]
            continue
        elif mouse < 11:
            set_num = f'S{(mouse-1)%5 + 1}'
            day0_sample = (expt, set_num, 0, 'vitro')
        else: # delay colonized, set day0 timepoint as community
            # set_num = f'S{(mouse-11)%5 + 1}' # set co
            day0_sample = (expt, 'comm', 0, 'vitro')

        print((expt, mouse), day0_sample)
    else: #in vitro
        barcode_inoc = well_to_bc_set[mouse]
        if barcode_inoc == '':
            print(f'Deleting timecourse for {host} because empty well')
            del timecourse_map[host]
            continue
        day0_sample = (expt, barcode_inoc, 0, 'vitro')

    day0_files = all_samples[day0_sample]
    rep1_files = [f for f in day0_files if 'rep1' in f]
    timecourse_map[host][0] = rep1_files[0][0]


## out
header = f'sampleID,Population,Timepoint,ExtractionBatch,ExtractionTube,Flagged?,Directory,FilenameStem\n'
out_file = args.output.open("w")
out_file.write(header)

sample_host_map = {}
for host, timecourse in timecourse_map.items():
    (expt, mouse) = host
    if expt == 'E1':
        host_id = f'E1_{mouse}'
    elif expt == 'E2':
        host_id = f'E2_{mouse}'
    else:
        host_id = f'vitro_{mouse}'

    for day, fastq in timecourse.items():
        if (fastq, day) not in sample_host_map:
            sample_host_map[(fastq, day)] = []
        sample_host_map[(fastq, day)].append(host_id)
    print(host_id)

sorted_samples = sorted(sample_host_map.keys())
for (sample, day) in sorted_samples:
    host_ids = sample_host_map[(sample, day)]
    host_str = ';'.join(host_ids)

    if 'E1' in host_str:
        directory = E1_dir 
    elif 'E2' in host_str:
        directory = E2_dir
    else:
        raise ValueError("No FASTQ directory configured for experiment: {}".format(host_str))

    out_line = f'{sample},{host_str},{day},NA,NA,0,{directory},{sample}\n'
    out_file.write(out_line)
# sorted_hosts = sorted(timecourse_map.keys())
# for host in sorted_hosts:
#     timecourse = timecourse_map[host]
#     (expt, isolator, mouse) = host
#     sampled_days = sorted(timecourse.keys())

#     out_line = f'{expt}\t{isolator}\t{mouse}\t'
#     for day in sampled_days:
#         fastq = timecourse[day]
#         out_line += f'\t({day},{fastq})'
#     out_line += '\n'
#     out_file.write(out_line)
# out_file.close()


# out_filename = "all_mgx_timecourse_lists.tsv"
# header = f'expt\tisolator\tmouse_or_well\t(day, fastq)\n'
# out_file = open(out_filename, "w")
# out_file.write(header)

# sorted_hosts = sorted(timecourse_map.keys())
# for host in sorted_hosts:
#     timecourse = timecourse_map[host]
#     (expt, isolator, mouse) = host
#     sampled_days = sorted(timecourse.keys())

#     out_line = f'{expt}\t{isolator}\t{mouse}\t'
#     for day in sampled_days:
#         fastq = timecourse[day]
#         out_line += f'\t({day},{fastq})'
#     out_line += '\n'
#     out_file.write(out_line)
# out_file.close()


# ### 
# out_filename = "all_population_samples.csv"
# out_file = open(out_filename, "w")
# out_file.close()
