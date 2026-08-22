#!/usr/bin/env python3
"""Create an isolate sample manifest from paired FASTQs and clone metadata."""

import argparse
import csv
import re
from pathlib import Path


FASTQ_RE = re.compile(r"(?P<stem>.+)_(?P<read>R[12])_001\.fastq(?:\.gz)?$")


def clone_tokens(clone_barcode_path):
    tokens = set()
    with clone_barcode_path.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            tokens.add("m{}_day{}_clone{}".format(
                row["mouse"], row["day"], row["clone"]
            ))
    return tokens


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("fastq_dir", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument(
        "--clone-barcode", type=Path,
        default=Path(__file__).with_name("clone_barcode.tsv"),
    )
    args = parser.parse_args()

    tokens = clone_tokens(args.clone_barcode)
    reads = {}
    for path in args.fastq_dir.iterdir():
        match = FASTQ_RE.fullmatch(path.name)
        if not match:
            continue
        stem, read = match.group("stem", "read")
        if not any(re.search(r"(?:^|_){}(?:_|$)".format(re.escape(token)), stem)
                       for token in tokens):
            continue
        reads.setdefault(stem, set()).add(read)

    incomplete = sorted(stem for stem, observed in reads.items()
                        if observed != {"R1", "R2"})
    if incomplete:
        raise SystemExit("Missing mate FASTQ for: {}".format(", ".join(incomplete)))

    samples = sorted(reads)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w") as handle:
        handle.write(str(args.fastq_dir.resolve()) + "\n")
        handle.write("\n".join(samples))
        if samples:
            handle.write("\n")
    print("Wrote {} paired isolate samples to {}".format(len(samples), args.output))


if __name__ == "__main__":
    main()
