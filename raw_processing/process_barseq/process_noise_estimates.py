import argparse
import pickle
from pathlib import Path


def parse_header(header):
    items = header.strip().split('_')
    cohort = items[0]
    if cohort == 'vitro':
        if 'passage0' in header:
            return cohort, items[1], 0
        return cohort, items[1] + '_' + items[2], int(items[3].lstrip('passage'))

    mouse = items[1].replace('M', '')
    if any(label in mouse for label in ('P', 'S', 'V')):
        if 'delay' in header:
            mouse = items[1].lstrip('M') + '_delay'
        day = 0
    else:
        mouse = int(mouse)
        day = 'cecum' if 'cec' in header else int(items[2].lstrip('D'))
    return cohort, mouse, day


def load_relabels(path):
    relabels = {}
    with path.open() as handle:
        next(handle)
        for line in handle:
            source, target = line.rstrip('\n').split('\t')
            source_items = source.strip('"').split(',')
            source_id = source_items[0], int(source_items[1]), int(source_items[2])
            target_items = target.strip('"').split(',')
            relabels[source_id] = None if target_items[0] == 'None' else (
                target_items[0], int(target_items[1]), int(target_items[2])
            )
    return relabels


def main():
    parser = argparse.ArgumentParser(description='Aggregate per-library noise estimates.')
    parser.add_argument('noise_estimate_dir', type=Path)
    parser.add_argument('--relabels', type=Path, default=Path(__file__).with_name('library_mislabelings.tsv'))
    parser.add_argument('--output', type=Path, default=Path('pickled_noise_inference.pkl'))
    args = parser.parse_args()

    relabels = load_relabels(args.relabels)

    effective_noise = {}
    for noise_file in sorted(args.noise_estimate_dir.glob('*_noise.out')):
        library_id = parse_header(noise_file.name)
        library_id = relabels.get(library_id, library_id)

        with noise_file.open() as handle:
            header_id = parse_header(next(handle))
            header_id = relabels.get(header_id, header_id)
            if library_id != header_id:
                raise ValueError('Filename and header disagree in {}: {} != {}'.format(
                    noise_file, library_id, header_id
                ))
            if library_id is None:
                continue
            if library_id in effective_noise:
                raise ValueError('Duplicate noise estimate for {}'.format(library_id))

            values = {}
            for line in handle:
                fields = line.lstrip('\n').split('\t')
                if len(fields) >= 2:
                    values[fields[0]] = fields[1]
        effective_noise[library_id] = values

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open('wb') as handle:
        pickle.dump(effective_noise, handle)
    print('Wrote {} libraries to {}'.format(len(effective_noise), args.output))


if __name__ == '__main__':
    main()
