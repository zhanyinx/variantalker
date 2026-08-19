#!/usr/bin/env python3

from collections import defaultdict
import argparse

def _parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i", 
        "--input", 
        required=True,
        help="Input VCF file"
    )
    parser.add_argument(
        "-p", 
        "--patient", 
        required=True,
        help="Patient identifier to use in output file names"
    )
    parser.add_argument(
        "-c", 
        "--chunk-size", 
        type=int, 
        required=True,
        help="Maximum number of variants per chunk"
    )

    args = parser.parse_args()
    return args

def main():
    """
    Read the input VCF, keep only PASS variants, group variants by
    chromosome and position, and create balanced chunks ensuring that
    variants from the same position are placed in different chunks.
    """
    args = _parse_args()

    header = []
    groups = defaultdict(list)

    # Group variants by chromosome and position
    with open(args.input) as f:
        for line in f:
            if line.startswith("#"):
                header.append(line)
                continue

            fields = line.rstrip().split("\t")

            # Skip non-PASS variants
            if fields[6] != "PASS":
                continue

            groups[(fields[0], fields[1])].append(line)

    if not groups:
        print("No PASS variants found.")
        return

    # Determine the number of chunks needed based on the largest group
    n_chunks = max(len(v) for v in groups.values())

    chunks = [[] for _ in range(n_chunks)]
    singletons = []

    # Distribute variants into chunks, putting variants at the same position in different chunks
    for vars_ in groups.values():

        if len(vars_) == 1:
            singletons.extend(vars_)

        else:
            for i, v in enumerate(vars_):
                chunks[i].append(v)

    sizes = [len(c) for c in chunks]

    # Distribute singletons (variants at unique positions) to balance the chunk sizes
    for v in singletons:
        idx = sizes.index(min(sizes))
        chunks[idx].append(v)
        sizes[idx] += 1

    # Split oversized chunks if needed
    final_chunks = []

    for chunk in chunks:
        if len(chunk) <= args.chunk_size:
            final_chunks.append(chunk)

        else:
            for i in range(0, len(chunk), args.chunk_size):
                final_chunks.append(chunk[i:i + args.chunk_size])

    for i, chunk in enumerate(final_chunks, 1):
        out = f"{args.patient}_chunk_{i}.vcf"
        with open(out, "w") as o:

            o.writelines(header)
            o.writelines(chunk)

        print(f"{out}\t{len(chunk)} variants")

if __name__ == "__main__":
    main()