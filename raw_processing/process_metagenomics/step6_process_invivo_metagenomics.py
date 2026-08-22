#!/usr/bin/env python3
"""Build notebook-ready in-vivo metagenomic mutation and coverage indexes.

Run after the merged mutation and coverage timecourses have been copied to the
configured local metagenomics directory. The output pickles are consumed by
``notebooks/4b_invivo_metagenomics.ipynb``.
"""

import argparse
import os
import pickle
import re
import sys
import tempfile
from pathlib import Path

import numpy as np


PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from project_config import load_config


EXCLUDE_GENES = {
    "BT1040/BT1041", "BT1042", "BT3239/BT3240", "BT4054/BT4055", "BT2672/BT2673",
} ## seem to have a lot of low-frequency mutations that are likely spurious,

MOB_GENES = ["BT{}".format(gene) for gene in range(3134, 3155)]
FIG4_SNP_GENES = [
    "BT0317", "BT0867", "BT1046", "BT1573", "BT3082", "BT3465", "BT4245",
    "BT4247", "BT4295", "BT0900",
]
PLOT_GENES = [
    "BT0867", "BT3465", "BT2660", "BT2661", "BT4014", "BT4015", "BT4540",
    "BT4543", "BT4520", "BT4523",
] + MOB_GENES + FIG4_SNP_GENES

COVERAGE_SPANS = {
    "BT0867": (1066079, 1066079),
    "BT3465": (1065108, 1068370),
    "BT3134-BT3154": (4001887 - 25, 4017683 + 25),
    "BT4014/BT4015": (5233175 - 25, 5233371 + 25),
    "BT2660/BT2661": (3316290 - 25, 3316499 + 25),
    "BT4540/BT4541": (5961716 - 25, 5962851 + 25),
    "BT4540/BT4543": (5961162 - 25, 5964550 + 25),
    "BT4520/BT4523": (5939745 - 25, 5943532 + 25),
}

MUTATION_FILE_RE = re.compile(r"^(E[12])_m?(\d+)_merged_timecourse$")
COVERAGE_FILE_RE = re.compile(r"^(E[12])_m?(\d+)_coverage_timecourse\.txt$")


def parse_gene_coordinates(ptt_path):
    """Return the historical PTT-derived gene coordinate map."""
    genes = {}
    with ptt_path.open() as handle:
        for _ in range(4):
            next(handle)
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            start, stop = (int(value) for value in fields[0].split(".."))
            genes[fields[5]] = (start, stop)
    return genes


def find_gene_from_position(position, sorted_genes):
    containing = [gene for gene, start, stop in sorted_genes if start <= position <= stop]
    if containing:
        return containing[0]

    upstream = [gene for gene, start, _ in sorted_genes if start <= position]
    downstream = [gene for gene, _, stop in sorted_genes if stop >= position]
    if upstream and downstream:
        return "{}-{}".format(upstream[-1], downstream[0])
    return None


def load_mutation_timecourse(path, plot_gene_locations, sorted_genes):
    mutations = {}
    with path.open() as handle:
        for line_number, line in enumerate(handle, start=1):
            fields = line.rstrip("\n").split(",")
            if len(fields) != 6:
                raise ValueError("{}:{}: expected six comma-separated fields".format(path, line_number))
            _, position, mutation, days, reads, coverage = fields
            position = int(position)
            if not np.any(np.abs(plot_gene_locations - position) <= 10000):
                continue

            gene = find_gene_from_position(position, sorted_genes)
            if gene in EXCLUDE_GENES:
                continue

            mutation = mutation.strip()
            days_array = np.array([float(value) for value in days.split() if value], dtype=int)
            reads_array = np.array([float(value) for value in reads.split()])
            coverage_array = np.array([float(value) for value in coverage.split()])
            if not (len(days_array) == len(reads_array) == len(coverage_array)):
                raise ValueError("{}:{}: mismatched timecourse lengths".format(path, line_number))

            mutations[position] = (mutation, reads_array, coverage_array)
            if gene is not None:
                mutations.setdefault(gene, []).append((position, mutation, reads_array, coverage_array))
            if "junction" in mutation:
                try:
                    junction = tuple(int(value) for value in mutation.split("_")[1:])
                except ValueError:
                    pass
                else:
                    mutations[junction] = (position, reads_array, coverage_array)
            mutations["days"] = days_array
    return mutations


def sample_time(sample_name):
    parts = sample_name.split("_")
    if len(parts) < 4:
        raise ValueError("Cannot determine timepoint from coverage header: {}".format(sample_name))
    value = parts[3]
    if "cecum" in value:
        return 100
    if "smallint" in value:
        return 200
    return int(value.replace("day", ""))


def load_coverage_timecourse(path):
    coverage = {
        key: {"single_site_coverages": [], "genome_span": span}
        for key, span in COVERAGE_SPANS.items()
    }
    with path.open() as handle:
        header = next(handle).rstrip("\n").split("\t")
        time_to_column = {}
        for index, sample in enumerate(header[2:], start=2):
            time_to_column[sample_time(sample)] = index
        sampled_times = np.array(sorted(time_to_column))
        sorted_columns = [time_to_column[time] for time in sampled_times]

        for line in handle:
            fields = line.rstrip("\n").split("\t")
            position = int(fields[1])
            for key, (start, stop) in COVERAGE_SPANS.items():
                if start <= position <= stop:
                    coverage[key]["single_site_coverages"].append(
                        np.array([float(fields[column]) for column in sorted_columns])
                    )

    for key, values in coverage.items():
        values["times"] = sampled_times
        values["average_coverage"] = np.mean(values.pop("single_site_coverages"), axis=0)
    return coverage


def write_pickle(value, output_path):
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile("wb", dir=output_path.parent, delete=False) as handle:
        pickle.dump(value, handle, protocol=pickle.HIGHEST_PROTOCOL)
        temporary_path = Path(handle.name)
    os.replace(temporary_path, output_path)


def main():
    config = load_config()
    default_input_dir = Path(config["local"]["mgx_timecourses"])
    if (default_input_dir / "mgx_timecourse").is_dir():
        default_input_dir /= "mgx_timecourse"
    default_output_dir = Path(config["local"]["pickles"])

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, default=default_input_dir)
    parser.add_argument("--output-dir", type=Path, default=default_output_dir)
    parser.add_argument(
        "--reference-ptt", type=Path,
        default=Path(config["local"]["data"]) / "reference_genome" / "BtVPI.ptt",
    )
    args = parser.parse_args()

    genes = parse_gene_coordinates(args.reference_ptt)
    missing_genes = sorted(set(PLOT_GENES) - set(genes))
    if missing_genes:
        raise ValueError("Plot genes missing from {}: {}".format(args.reference_ptt, ", ".join(missing_genes)))
    plot_gene_locations = np.array([genes[gene] for gene in PLOT_GENES])
    sorted_genes = sorted((gene, start, stop) for gene, (start, stop) in genes.items())

    mutation_timecourses = {}
    coverage_timecourses = {}
    for path in sorted(args.input_dir.iterdir()):
        mutation_match = MUTATION_FILE_RE.match(path.name)
        if mutation_match:
            key = (mutation_match.group(1), int(mutation_match.group(2)))
            mutation_timecourses[key] = load_mutation_timecourse(path, plot_gene_locations, sorted_genes)
            continue
        coverage_match = COVERAGE_FILE_RE.match(path.name)
        if coverage_match:
            key = (coverage_match.group(1), int(coverage_match.group(2)))
            coverage_timecourses[key] = load_coverage_timecourse(path)

    if not mutation_timecourses:
        raise ValueError("No merged mutation timecourses found in {}".format(args.input_dir))
    if not coverage_timecourses:
        raise ValueError("No coverage timecourses found in {}".format(args.input_dir))

    write_pickle(mutation_timecourses, args.output_dir / "metagenomic_mutation_timecourses.pkl")
    write_pickle(coverage_timecourses, args.output_dir / "metagenomic_coverage_timecourses.pkl")
    print("Wrote {} mutation and {} coverage timecourses to {}".format(
        len(mutation_timecourses), len(coverage_timecourses), args.output_dir
    ))


if __name__ == "__main__":
    main()
