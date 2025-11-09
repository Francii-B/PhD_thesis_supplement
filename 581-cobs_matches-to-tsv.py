#!/usr/bin/env python3
import os
import gzip
import argparse
import csv

#Example usage: 581-cobs_matches-to-tsv.py <input_dir> <output_tsv>

def process_gz_file(filepath):
    """Extract pairs of (*, _) lines from a gzipped file, removing prefixes."""
    results = []
    current_star = None

    with gzip.open(filepath, 'rt', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('*'):
                # Remove '*' and get first tab-separated field
                current_star = line.lstrip('*').split('\t')[0]
            elif line.startswith('_') and current_star is not None:
                # Remove '_' and get first tab-separated field
                underscore_val = line.lstrip('_').split('\t')[0]
                results.append((current_star, underscore_val))
    return results


def main():
    parser = argparse.ArgumentParser(description="Convert .gz files into a combined TSV file.")
    parser.add_argument("input_dir", help="Directory containing .gz files")
    parser.add_argument("output_file", help="Output TSV file path")
    args = parser.parse_args()

    all_results = []

    # Iterate over .gz files in directory
    for filename in sorted(os.listdir(args.input_dir)):
        if filename.endswith(".gz"):
            filepath = os.path.join(args.input_dir, filename)
            file_results = process_gz_file(filepath)
            all_results.extend(file_results)

    # Write output TSV
    with open(args.output_file, 'w', newline='', encoding='utf-8') as out_f:
        writer = csv.writer(out_f, delimiter='\t')
       # writer.writerow(['start_line', 'underscore_line'])
        writer.writerows(all_results)


if __name__ == "__main__":
    main()
