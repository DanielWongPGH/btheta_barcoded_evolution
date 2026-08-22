#!/usr/bin/env python3
"""Identify isolate mutations and candidate drivers from annotated breseq GD files."""

import argparse
import copy
import pickle
import sys
from pathlib import Path

import numpy as np
import scipy.stats


PROJECT_ROOT = Path(__file__).resolve().parents[2]
NOTEBOOKS_DIR = PROJECT_ROOT / "notebooks"
for path in (PROJECT_ROOT, NOTEBOOKS_DIR):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

from project_config import CONFIG_FILE, load_config
from methods import genomics_util as genomics
from methods import shared


EXCLUSION_LIST = {
    "BT1040/BT1041", "BT1042", "BT3239/BT3240", "BT4054/BT4055",
    "BT2672/BT2673",
}
SIMPLE_TYPES = {"SNP", "SUB", "DEL", "INS"}
STRUCTURAL_TYPES = {"insertion", "IR-inversion", "SV", "unpaired_JC"}
SAMPLE_DAYS = (9, 17, 36, 51)


def load_pickle(path):
    with path.open("rb") as handle:
        return pickle.load(handle)


def dump_pickle(value, path):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("wb") as handle:
        pickle.dump(value, handle)


def filter_mutations_in_sample(
    all_mutations, sample, evidence, preexisting_variants, exclude_genes=None,
    sv_min_freq=0.1, sv_min_size=10, sv_min_cov=5, polymorphism_min=None,
):
    exclude_genes = set() if exclude_genes is None else exclude_genes
    for mutation_type in SIMPLE_TYPES:
        for position, _, mutation in evidence[mutation_type]:
            if polymorphism_min is not None and mutation["frequency"] < polymorphism_min:
                continue
            frequency = mutation["frequency"] if polymorphism_min is not None else 1
            information = [mutation[key] for key in (
                "gene", "strand", "gene_description", "PUL", "mutation",
                "mutation_description", "syn",
            )]
            key = tuple([position, mutation_type] + information)
            if genomics.check_preexisting(key, preexisting_variants, exclude_genes=exclude_genes):
                continue
            all_mutations.setdefault(key, {})[sample] = frequency

    for mutation in evidence["JC"]:
        (start, end), genes, orientations, frequency, coverage = mutation
        if frequency < sv_min_freq or end - start < sv_min_size or coverage < sv_min_cov:
            continue
        key = ((start, end), "JC", genes, orientations)
        if not genomics.check_preexisting(key, preexisting_variants, exclude_genes=exclude_genes):
            all_mutations.setdefault(key, {})[sample] = (frequency, coverage)

    for mutation in evidence["MC"]:
        (start, end), genes = mutation
        if end - start < sv_min_size:
            continue
        key = ((start, end), "MC", genes)
        if not genomics.check_preexisting(key, preexisting_variants, exclude_genes=exclude_genes):
            all_mutations.setdefault(key, {})[sample] = 0


def load_isolate_evidence(gd_dir, clone_barcode_map):
    evidence_map = {}
    for gd_file in sorted(gd_dir.glob("*.gd")):
        clone = genomics.parse_isolate_gd_filename(gd_file.name)
        if clone not in clone_barcode_map:
            continue
        coverage, evidence = genomics.get_variants_from_annotated_gd_file(gd_file)
        evidence_map[clone] = {
            "barcode": clone_barcode_map[clone], "coverage": coverage,
            "evidence": evidence,
        }
    return evidence_map


def identify_preexisting_variants(
    day0_gd_dir, isolate_evidence, min_day0_frequency=0.1,
    min_clone_prevalence=0.9,
):
    day0_mutations = {}
    for gd_file in sorted(day0_gd_dir.glob("*.gd")):
        replicate = gd_file.name.replace("_day0_annotated.gd", "")
        _, evidence = genomics.get_variants_from_annotated_gd_file(gd_file, freq=True)
        filter_mutations_in_sample(
            day0_mutations, replicate, evidence, [],
            polymorphism_min=min_day0_frequency,
        )

    day0_mutations = {
        mutation: replicates for mutation, replicates in day0_mutations.items()
        if len(replicates) >= 2
    }
    all_isolate_mutations = {}
    for clone, clone_data in isolate_evidence.items():
        filter_mutations_in_sample(
            all_isolate_mutations, clone, clone_data["evidence"], [],
        )

    minimum_count = min_clone_prevalence * len(isolate_evidence)
    preexisting = {
        mutation for mutation, clones in all_isolate_mutations.items()
        if len(clones) >= minimum_count
    }
    preexisting.update(day0_mutations)
    return list(preexisting)


def collect_filtered_mutations(isolate_evidence, preexisting_variants):
    mutations = {}
    for clone, clone_data in isolate_evidence.items():
        filter_mutations_in_sample(
            mutations, clone, clone_data["evidence"], preexisting_variants,
        )
    return mutations


def pair_junctions(all_mutations):
    junctions = sorted(
        ((mutation, clones) for mutation, clones in all_mutations.items()
         if mutation[1] == "JC"),
        key=lambda item: item[0][0],
    )
    paired = {}
    for index, (mutation, mutation_clones) in enumerate(junctions[:-1]):
        start, stop = mutation[0]
        orientations = np.asarray(mutation[-1])
        for next_mutation, next_clones in junctions[index + 1:]:
            next_start, next_stop = next_mutation[0]
            next_orientations = np.asarray(next_mutation[-1])
            shared_clones = set(mutation_clones).intersection(next_clones)
            if not shared_clones:
                continue
            if len(shared_clones) < 0.5 * min(len(mutation_clones), len(next_clones)):
                continue
            if not (abs(start - next_start) < 50 or abs(stop - next_stop) < 50):
                continue

            inversion = genomics.find_spanning_inverted_repeats(
                mutation[0], next_mutation[0], 3
            )
            mutation_type = "SV"
            if (inversion and abs(start - next_start) < 50
                    and abs(stop - next_stop) < 50
                    and (orientations * next_orientations).sum() == -2):
                mutation_type = "IR-inversion"
            elif abs(start - stop) > 10000 or abs(next_start - next_stop) > 10000:
                mutation_type = "insertion"

            key = tuple([tuple(sorted((mutation, next_mutation), key=lambda value: value[0]))])
            paired[key] = [list(shared_clones), mutation_type]

    flattened = {junction for pair in paired for junction in pair[0]}
    merged = [(pair, clones, mutation_type) for pair, (clones, mutation_type) in paired.items()]
    for mutation, clones in junctions:
        if mutation not in flattened:
            merged.append(((mutation,), list(clones), "unpaired_JC"))
    return sorted(merged, key=lambda item: item[0][0])


def adjust_bh(pvalues):
    adjusted = {}
    previous = 1.0
    ranked = sorted(pvalues.items(), key=lambda item: item[1], reverse=True)
    total = len(ranked)
    for reverse_rank, (key, value) in enumerate(ranked, start=1):
        rank = total - reverse_rank + 1
        previous = min(previous, value * total / rank, 1.0)
        adjusted[key] = previous
    return adjusted


def estimate_sv_pvalues(
    merged_junctions, clone_barcode_map, barcode_clone_map, rng,
    max_samples=1_000_000, batch_size=100_000,
):
    results = {}
    raw_pvalues = {}
    for junction, clones, mutation_type in merged_junctions:
        barcodes = set(clone_barcode_map[clone] for clone in clones)
        structural_variant = (*junction, mutation_type)
        states = []
        for barcode, barcode_data in barcode_clone_map.items():
            for day in SAMPLE_DAYS:
                if day not in barcode_data:
                    continue
                day_clones = barcode_data[day]
                present = sum(clone in clones for clone in day_clones)
                states.append((len(day_clones) - present, present))
        states = np.asarray(states)
        if not len(states) or (states[:, 1] > 0).sum() <= 1:
            results[structural_variant] = (clones, barcodes, np.nan, False)
            continue

        clones_per_barcode = states.sum(axis=1)
        clones_with_junction = states[:, 1].sum()
        distribution = scipy.stats.multivariate_hypergeom(
            clones_per_barcode, clones_with_junction, seed=rng
        )
        observed = distribution.logpmf(states[:, 1])
        extreme = 0
        sampled = 0
        while sampled < max_samples and extreme / max(sampled, 1) < 0.05:
            size = min(batch_size, max_samples - sampled)
            simulated = distribution.rvs(size=size, random_state=rng)
            extreme += np.count_nonzero(distribution.logpmf(simulated) <= observed)
            sampled += size
        raw_pvalue = max(1 / max_samples, extreme / sampled)
        raw_pvalues[structural_variant] = raw_pvalue
        results[structural_variant] = (clones, barcodes, raw_pvalue, True)

    for structural_variant, adjusted in adjust_bh(raw_pvalues).items():
        clones, barcodes, _, tested = results[structural_variant]
        results[structural_variant] = (clones, barcodes, adjusted, tested)
    return results


def add_sv_frequencies(sv_pvalues, all_mutations):
    processed = {}
    for structural_variant, (clones, barcodes, pvalue, tested) in sv_pvalues.items():
        if structural_variant[-1] == "unpaired_JC":
            junction = structural_variant[0]
            frequencies = [all_mutations[junction][clone][0] for clone in clones]
        else:
            first, second = structural_variant[0]
            frequencies = []
            for clone in clones:
                first_frequency, first_coverage = all_mutations[first][clone]
                second_frequency, second_coverage = all_mutations[second][clone]
                frequencies.append(
                    (first_coverage + second_coverage)
                    / (first_coverage / first_frequency + second_coverage / second_frequency)
                )
        processed[structural_variant] = (
            clones, barcodes, pvalue, tested, np.median(frequencies)
        )
    return processed


def select_drivers(
    ordered_barcodes, barcode_clone_map, all_mutations, processed_svs,
    fdr_alpha=0.05,
):
    drivers = {}
    for barcode, barcode_clones in ordered_barcodes:
        clones = barcode_clone_map[barcode]["all"]
        barcode_mutations = {}
        for mutation, mutation_clones in all_mutations.items():
            if mutation[1] in {"JC", "MC"}:
                continue
            detected = [clone for clone in mutation_clones if clone in clones]
            if detected:
                barcode_mutations[mutation] = {
                    day: [clone for clone in detected if clone[1] == day]
                    for day in SAMPLE_DAYS
                }
                barcode_mutations[mutation]["all"] = detected

        for structural_variant, data in processed_svs.items():
            sv_clones, sv_barcodes, pvalue, _, median_frequency = data
            detected = [clone for clone in sv_clones if clone in clones]
            if not detected or median_frequency < 0.5:
                continue
            if not np.isnan(pvalue) and pvalue > fdr_alpha and len(sv_barcodes) > 1:
                continue
            barcode_mutations[structural_variant] = {
                day: [clone for clone in detected if clone[1] == day]
                for day in SAMPLE_DAYS
            }
            barcode_mutations[structural_variant]["all"] = detected

        for mutation, clone_data in barcode_mutations.items():
            if len(clone_data["all"]) <= 0.5 * len(clones):
                continue
            drivers.setdefault(mutation, {})[barcode] = clone_data
    return drivers


def describe_junction(structural_variant):
    if structural_variant[-1] == "unpaired_JC":
        junction = structural_variant[0]
        paired_evidence = ""
        mutation_type = "unpaired JC"
    else:
        junction, paired = structural_variant[0]
        paired_evidence = "{}_{}_{}_{}".format(
            paired[0][0], paired[-1][0], paired[0][1], paired[-1][1]
        )
        mutation_type = structural_variant[-1]
    locus1, locus2 = junction[2]
    descriptions = []
    pul_labels = []
    for locus in (locus1, locus2):
        genes = locus.split("-")
        descriptions.append("/".join(shared.gene_description[gene] for gene in genes))
        pul_labels.append("/".join(
            str(genomics.gene_PUL_map[gene][0])
            for gene in genes if gene in genomics.gene_PUL_map
        ))
    mutation_description = "{}_{}_{}_{}".format(
        junction[0][0], junction[-1][0], junction[0][1], junction[-1][1]
    )
    return (
        f"{locus1} - {locus2}", " - ".join(descriptions),
        " - ".join(pul_labels).strip(" -"), mutation_type,
        mutation_description, paired_evidence,
    )


def write_driver_table(path, drivers, ordered_barcodes, barcode_clone_map):
    rows = []
    for mutation in drivers:
        if mutation[-1] in STRUCTURAL_TYPES:
            if mutation[-1] == "unpaired_JC":
                position = mutation[0][0][0]
            else:
                position = mutation[0][0][0][0]
            fields = describe_junction(mutation)
        else:
            position = mutation[0]
            fields = (
                mutation[2], mutation[4], str(mutation[5]), mutation[1],
                mutation[7], "",
            )
        rows.append((position, mutation, fields))

    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        handle.write(
            "position\tlocus\tlocus description\tPUL\tmutation or evidence\t"
            "mutation description\tpaired evidence (if junction)\tn_barcodes\t"
            "barcode (detected/sampled)\n"
        )
        for position, mutation, fields in sorted(rows, key=lambda row: row[0]):
            barcode_labels = []
            for barcode, _ in ordered_barcodes:
                if barcode in drivers[mutation]:
                    detected = len(drivers[mutation][barcode]["all"])
                    sampled = len(barcode_clone_map[barcode]["all"])
                    barcode_labels.append(f"{barcode} ({detected}/{sampled})")
            handle.write("\t".join(map(str, (
                position, *fields, len(drivers[mutation]), ", ".join(barcode_labels)
            ))) + "\n")


def write_sv_table(path, processed_svs, barcode_clone_map):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        handle.write(
            "locus1\tlocus2\ttype\tn_barcodes represented\tn_clones\t"
            "frac. of clones in barcodes represented\tFDR-corrected pval\t"
            "median_freq\tbarcodes\n"
        )
        for structural_variant, data in processed_svs.items():
            clones, barcodes, pvalue, _, median_frequency = data
            total = sum(len(barcode_clone_map[barcode]["all"]) for barcode in barcodes)
            if structural_variant[-1] == "unpaired_JC":
                locus1, locus2 = structural_variant[0], ""
            else:
                locus1, locus2 = structural_variant[0]
            pvalue_label = "nan" if np.isnan(pvalue) else f"{pvalue:.2e}"
            handle.write("\t".join(map(str, (
                locus1, locus2, structural_variant[-1], len(barcodes), len(clones),
                f"{len(clones) / total:.2f}", pvalue_label,
                f"{median_frequency:.2f}", ", ".join(barcodes),
            ))) + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=CONFIG_FILE)
    parser.add_argument("--isolate-gd-dir", type=Path)
    parser.add_argument("--day0-gd-dir", type=Path)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--max-samples", type=int, default=1_000_000)
    parser.add_argument("--batch-size", type=int, default=100_000)
    args = parser.parse_args()

    config = load_config(args.config)["local"]
    data_dir = Path(config["data"])
    pickled_dir = Path(config["pickles"])
    tables_dir = Path(config["tables"])
    isolate_gd_dir = args.isolate_gd_dir or data_dir / "isolate_wgs/E1_clones_gd_annotated"
    day0_gd_dir = args.day0_gd_dir or data_dir / "reference_genome/day0_gd_annotated"
    for directory in (isolate_gd_dir, day0_gd_dir):
        if not directory.is_dir():
            raise FileNotFoundError(directory)

    barcode_clone_map = load_pickle(pickled_dir / "barcode_clone_map.pkl")
    clone_barcode_map = load_pickle(pickled_dir / "clone_barcode_map.pkl")
    ordered_barcodes = sorted(
        ((barcode, clones["all"]) for barcode, clones in barcode_clone_map.items()
         if len(barcode) == 20),
        key=lambda item: len(item[1]), reverse=True,
    )

    isolate_evidence = load_isolate_evidence(isolate_gd_dir, clone_barcode_map)
    preexisting = identify_preexisting_variants(day0_gd_dir, isolate_evidence)
    mutations = collect_filtered_mutations(isolate_evidence, preexisting)
    merged_junctions = pair_junctions(mutations)
    sv_pvalues = estimate_sv_pvalues(
        merged_junctions, clone_barcode_map, barcode_clone_map,
        np.random.default_rng(args.seed), args.max_samples, args.batch_size,
    )
    processed_svs = add_sv_frequencies(sv_pvalues, mutations)
    drivers = select_drivers(
        ordered_barcodes, barcode_clone_map, mutations, processed_svs
    )

    dump_pickle(preexisting, pickled_dir / "preexisting_variants.pkl")
    dump_pickle(sv_pvalues, pickled_dir / "SV_pvals.pkl")
    dump_pickle(mutations, pickled_dir / "all_mutations_in_isolates.pkl")
    dump_pickle(drivers, pickled_dir / "all_vivo_drivers.pkl")
    write_driver_table(
        tables_dir / "tableS1_putative_driver_mutations.tsv",
        drivers, ordered_barcodes, barcode_clone_map,
    )
    write_sv_table(
        tables_dir / "tableS2_isolates_structural_variants.tsv",
        processed_svs, barcode_clone_map,
    )
    print(f"Processed {len(isolate_evidence)} isolates")
    print(f"Retained {len(mutations)} mutations and {len(drivers)} candidate drivers")


if __name__ == "__main__":
    main()
