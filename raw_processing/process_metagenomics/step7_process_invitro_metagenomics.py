#!/usr/bin/env python3
"""Process annotated in-vitro breseq GD files into Table S6 and notebook data.

The compact pickle written by this step is loaded by the retained
``step7_process_invitro_metagenomics.ipynb`` diagnostics notebook.
"""

import argparse
import os
import pickle
import sys
import tempfile
from pathlib import Path

import numpy as np


PROJECT_ROOT = Path(__file__).resolve().parents[2]
NOTEBOOKS_DIR = PROJECT_ROOT / "notebooks"
if str(NOTEBOOKS_DIR) not in sys.path:
    sys.path.insert(0, str(NOTEBOOKS_DIR))
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from project_config import load_config
from methods import genomics_util as breseq
from methods import shared


EXCLUDED_WELLS = {"p2_E3", "p2_F3", "p2_G3", "p2_H3"}
ANNOTATED_SETS = ["S{}_annotated".format(index) for index in range(1, 7)]


def find_annotated_gds(gd_directory):
    """Map (well, passage) to a unique annotated GD file."""
    gd_files = {}
    for vitro_set in ANNOTATED_SETS:
        set_directory = gd_directory / vitro_set
        if not set_directory.is_dir():
            raise FileNotFoundError("Annotated GD directory not found: {}".format(set_directory))
        for path in sorted(set_directory.iterdir()):
            if path.suffix != ".gd" or "d0" in path.name:
                continue
            fields = path.name.split("_")
            if len(fields) < 5:
                raise ValueError("Unexpected annotated GD filename: {}".format(path))
            passage = int(fields[2].replace("passage", ""))
            well = "{}_{}".format(fields[3], fields[4].replace("_annotated.gd", ""))
            if well in EXCLUDED_WELLS:
                continue
            key = (well, passage)
            if key in gd_files:
                raise ValueError("Duplicate annotated GD for {}: {} and {}".format(key, gd_files[key], path))
            gd_files[key] = path
    if not gd_files:
        raise ValueError("No annotated GD files found beneath {}".format(gd_directory))
    return gd_files


def filter_mutations_in_sample(
    all_mutations, sample, evidence, preexisting_variants, exclude_genes=None,
    sv_min_frequency=0.0, sv_min_size=10, sv_min_coverage=5, polymorphism_minimum=0.0,
):
    exclude_genes = exclude_genes or set()
    for evidence_type in ("SNP", "SUB", "DEL", "INS", "RA"):
        for position, _, mutation_data in evidence[evidence_type]:
            mutation_information = [
                mutation_data[key] for key in (
                    "gene", "strand", "gene_description", "PUL", "mutation",
                    "mutation_description", "syn",
                )
            ]
            variant_type = mutation_data["mutation_category"]
            frequency = mutation_data["frequency"]
            if frequency < polymorphism_minimum:
                continue

            key = tuple([position, variant_type] + mutation_information)
            if breseq.check_preexisting(key, preexisting_variants, exclude_genes=exclude_genes):
                continue
            if evidence_type == "RA" and key in all_mutations and sample in all_mutations[key]:
                continue
            if evidence_type == "RA":
                key = tuple([position, "RA"] + mutation_information)
            all_mutations.setdefault(key, {})[sample] = frequency

    for (start, end), genes, orientations, frequency, junction_coverage in evidence["JC"]:
        if frequency < sv_min_frequency or end - start < sv_min_size or junction_coverage < sv_min_coverage:
            continue
        key = ((start, end), "JC", genes, orientations)
        if not breseq.check_preexisting(key, preexisting_variants, exclude_genes=exclude_genes):
            all_mutations.setdefault(key, {})[sample] = (frequency, junction_coverage)

    for (start, end), genes in evidence["MC"]:
        if end - start < sv_min_size:
            continue
        key = ((start, end), "MC", genes)
        if not breseq.check_preexisting(key, preexisting_variants, exclude_genes=exclude_genes):
            all_mutations.setdefault(key, {})[sample] = 0


def merge_mutations_in_wells(all_mutations):
    filtered = {
        medium: {well: {} for well, _ in shared.medium_to_well_map[medium]}
        for medium in shared.ordered_media
    }
    mutation_flags = {}
    for mutation, sample_frequencies in all_mutations.items():
        sampled_wells = {well for well, _ in sample_frequencies}
        for well, passage in sample_frequencies:
            medium, _ = shared.well_to_medium_map[well]
            filtered[medium][well].setdefault(mutation, {})[passage] = sample_frequencies[(well, passage)]

        for well in sampled_wells:
            medium, _ = shared.well_to_medium_map[well]
            passage_frequencies = filtered[medium][well][mutation]
            if 15 not in passage_frequencies or 31 not in passage_frequencies:
                del filtered[medium][well][mutation]
                continue
            mutation_flags[mutation] = True

    all_filtered = {}
    for medium, wells in filtered.items():
        for well, mutations in wells.items():
            for mutation, passages in mutations.items():
                if mutation_flags.get(mutation):
                    all_filtered.setdefault(mutation, {}).update(
                        {(well, passage): frequency for passage, frequency in passages.items()}
                    )
    return filtered, all_filtered


def mutation_table_row(mutation):
    mutation_type = mutation[1]
    if mutation_type == "JC":
        location_1, location_2 = mutation[0]
        orientation_1, orientation_2 = mutation[-1]
        gene_1, gene_2 = mutation[2]
        genes, descriptions, puls = [], [], []
        for gene in (gene_1, gene_2):
            split_genes = gene.split("-")
            genes.append("/".join(split_genes))
            descriptions.append("/".join(shared.gene_description[item] for item in split_genes))
            puls.extend(str(breseq.gene_PUL_map[item][0]) for item in split_genes if item in breseq.gene_PUL_map)
        return (
            str(location_1), " - ".join(genes), " - ".join(descriptions),
            ",".join(sorted(set(puls))), mutation_type,
            "{}_{}_{}_{}".format(location_1, orientation_1, location_2, orientation_2),
        )
    if mutation_type == "MC":
        return None
    return (
        str(mutation[0]), str(mutation[2]), str(mutation[4]), str(mutation[5]),
        mutation_type, str(mutation[-2]),
    )


def write_table(filtered, output_path):
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile("w", dir=output_path.parent, delete=False) as handle:
        header = ["position", "gene", "gene_description", "PUL", "mut_type", "mutation", "freq. in (Medium, Well)"]
        well_order = [(medium, well) for medium in shared.ordered_media for well, _ in shared.medium_to_well_map[medium]]
        handle.write("\t".join(header + ["({}, {})".format(medium, well) for medium, well in well_order]) + "\n")
        for mutation in sorted(filtered, key=lambda value: value[0] if isinstance(value[0], int) else value[0][0]):
            row = mutation_table_row(mutation)
            if row is None:
                continue
            values = list(row)
            for medium, well in well_order:
                trajectory = filtered[mutation]
                passage_values = {passage: frequency for (sample_well, passage), frequency in trajectory.items() if sample_well == well}
                if mutation[1] == "JC":
                    formatted = ", ".join("{}:{:.2f}".format(passage, passage_values[passage][0]) for passage in sorted(passage_values))
                else:
                    formatted = ", ".join("{}:{:.2f}".format(passage, passage_values[passage]) for passage in sorted(passage_values))
                values.append(formatted)
            handle.write("\t".join(values) + "\n")
        temporary_path = Path(handle.name)
    os.replace(temporary_path, output_path)


def write_pickle(value, output_path):
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile("wb", dir=output_path.parent, delete=False) as handle:
        pickle.dump(value, handle, protocol=pickle.HIGHEST_PROTOCOL)
        temporary_path = Path(handle.name)
    os.replace(temporary_path, output_path)


def main():
    config = load_config()
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gd-directory", type=Path, default=Path(config["local"]["invitro_gd"]))
    parser.add_argument("--output-dir", type=Path, default=Path(config["local"]["pickles"]))
    parser.add_argument("--table-output", type=Path, default=Path(config["local"]["tables"]) / "tableS4_invitro_mutations.tsv")
    args = parser.parse_args()

    with (Path(config["local"]["pickles"]) / "preexisting_variants.pkl").open("rb") as handle:
        preexisting_variants = pickle.load(handle)

    annotated_gds = find_annotated_gds(args.gd_directory)
    well_mutations = {
        medium: {well: {} for well, _ in shared.medium_to_well_map[medium]}
        for medium in shared.medium_to_well_map
    }
    all_mutations = {}
    for (well, passage), gd_path in annotated_gds.items():
        medium, _ = shared.well_to_medium_map[well]
        coverage, evidence = breseq.get_variants_from_annotated_gd_file(gd_path, freq=True, track_RA=False)
        well_mutations[medium][well][passage] = {"coverage": coverage, "evidence": evidence}
        filter_mutations_in_sample(all_mutations, (well, passage), evidence, preexisting_variants)

    filtered_well_mutations, filtered_mutations = merge_mutations_in_wells(all_mutations)
    output = {
        "well_mutations_map": well_mutations,
        "all_mutations_in_vitro": all_mutations,
        "filtered_well_mutations": filtered_well_mutations,
        "all_mutations": filtered_mutations,
    }
    write_pickle(output, args.output_dir / "invitro_metagenomic_mutations.pkl")
    write_table(filtered_mutations, args.table_output)
    print("Processed {} annotated GD files into {} retained mutations".format(len(annotated_gds), len(filtered_mutations)))


if __name__ == "__main__":
    main()
