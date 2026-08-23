#!/usr/bin/env python3
"""Validate and collect per-library Bar-seq noise-model fits.

The output is keyed by the canonical three-part library identifiers used by
the analysis notebooks: ``(experiment, sample, timepoint)``. Step 5 imports
the parsing helpers here so it applies the same relabeling rules.
"""

import argparse
import pickle
from pathlib import Path


CECUM_DAY = 100
REQUIRED_PARAMETERS = ('p', 'Z1', 'Z2', 'D1', 'D2')


def load_relabels(path):
    relabels = {}
    with path.open() as handle:
        next(handle)
        for line in handle:
            source, target = line.rstrip('\n').split('\t')
            source = [field.strip() for field in source.strip('"').split(',')]
            target = [field.strip() for field in target.strip('"').split(',')]
            source_id = (source[0], int(source[1]), int(source[2]), int(source[3]))
            relabels[source_id] = None if target[0] == 'None' else (
                target[0], int(target[1]), int(target[2]), int(target[3])
            )
    return relabels


def parse_fastq_library(stem):
    """Parse a noise-file stem into the four-part count-table identifier."""
    fields = stem.split('_')
    experiment = fields[0]
    if experiment == 'E1':
        if 'tube2' in stem or 'S2for9_28' in stem:
            return None
        if 'Plate' in fields[1]:
            return ('E1', fields[1].replace('Plate', 'P'), 0, 1)
        day = CECUM_DAY if 'cec' in stem else int(fields[2].removeprefix('day'))
        return ('E1', int(fields[1].removeprefix('m')), day, 1)
    if experiment == 'E2':
        if '118-S' in stem:
            return None
        if 'day0' in stem:
            sample = fields[1]
            if '116-S' in sample:
                sample = sample.replace('116-S', 'S') + 'delay'
            return ('E2', sample, 0, 1)
        day = CECUM_DAY if 'cec' in stem else int(fields[2].removeprefix('day'))
        return ('E2', int(fields[1].removeprefix('m')), day, 1)
    if experiment == 'EV':
        passage = int(fields[1].removeprefix('passage'))
        sample = fields[2] if 'V' in fields[2] else fields[2].replace('Plate', 'p') + '-' + fields[3]
        return ('EV', sample, passage, 1)
    raise ValueError('Unrecognized Bar-seq noise-file stem: {}'.format(stem))


def canonical_noise_id(stem, relabels):
    """Return the notebook-facing three-part ID for a noise-file stem."""
    library_id = parse_fastq_library(stem)
    if library_id is None:
        return None
    library_id = relabels.get(library_id, library_id)
    if library_id is None:
        return None
    experiment, sample, timepoint, _replicate = library_id
    if experiment == 'EV':
        experiment = 'vitro'
        sample = sample.replace('-', '_')
    elif isinstance(sample, str):
        sample = sample.replace('delay', '_delay')
    if timepoint == CECUM_DAY:
        timepoint = 'cecum'
    return experiment, sample, timepoint


def read_noise_file(path):
    stem = path.name.removesuffix('_noise.out')
    with path.open() as handle:
        header = next(handle).strip()
        if header != stem:
            raise ValueError('Filename and header disagree in {}: {} != {}'.format(path, stem, header))
        values = {}
        for line in handle:
            if '\t' in line:
                label, value = line.split('\t', 1)
                values[label] = value
    missing = [parameter for parameter in REQUIRED_PARAMETERS if parameter not in values]
    if missing:
        raise ValueError('{} is missing required parameters: {}'.format(path, ', '.join(missing)))
    for parameter in REQUIRED_PARAMETERS:
        float(values[parameter])
    return values


def main():
    parser = argparse.ArgumentParser(description='Validate and aggregate per-library noise estimates.')
    parser.add_argument('noise_estimate_dir', type=Path)
    parser.add_argument('--relabels', type=Path, default=Path(__file__).with_name('library_mislabelings.tsv'))
    parser.add_argument('--output', type=Path, default=Path('pickled_noise_inference.pkl'))
    args = parser.parse_args()
    relabels = load_relabels(args.relabels)
    effective_noise = {}
    for noise_file in sorted(args.noise_estimate_dir.glob('*_noise.out')):
        library_id = canonical_noise_id(noise_file.name.removesuffix('_noise.out'), relabels)
        if library_id is None:
            continue
        if library_id in effective_noise:
            raise ValueError('Duplicate noise estimate for {}'.format(library_id))
        effective_noise[library_id] = read_noise_file(noise_file)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open('wb') as handle:
        pickle.dump(effective_noise, handle)
    print('Wrote {} libraries to {}'.format(len(effective_noise), args.output))


if __name__ == '__main__':
    main()
