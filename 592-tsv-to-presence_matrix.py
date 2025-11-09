#!/usr/bin/env python3
import os
import sys
import pandas as pd
import glob

def main(input_dir, output_file):
    # Collect all TSV files
    tsv_files = glob.glob(os.path.join(input_dir, "*.tsv"))
    if not tsv_files:
        print("No TSV files found in the directory.")
        sys.exit(1)

    subject_sets = {}
    all_subjects = set()

    # Process each file
    for file in tsv_files:
        # Extract query name from filename (without extension)
        query_name = os.path.splitext(os.path.basename(file))[0]

        # Read first two columns only
        df = pd.read_csv(file, sep="\t", usecols=[0, 1], header=None, names=["query", "subject"])

        # Collect subjects for this query
        subjects = set(df["subject"].astype(str).tolist())
        subject_sets[query_name] = subjects

        # Keep track of all subjects seen
        all_subjects.update(subjects)

    # Build the matrix
    all_subjects = sorted(all_subjects)
    matrix = pd.DataFrame(0, index=all_subjects, columns=subject_sets.keys())

    for query, subjects in subject_sets.items():
        matrix.loc[list(subjects), query] = 1

    # Save output
    matrix.to_csv(output_file, sep="\t")

    print(f"Matrix saved to {output_file}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python 592-tsv-to-presence_matrix.py <input_dir> <output_file>")
        sys.exit(1)
    input_dir = sys.argv[1]
    output_file = sys.argv[2]
    main(input_dir, output_file)

