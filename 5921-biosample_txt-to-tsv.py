import pandas as pd
import re
import argparse

#Example usage: python 5921-biosample_txt-to-tsv.py input.txt output.tsv

# --- Command line arguments ---
parser = argparse.ArgumentParser(description="Convert bacterial records txt to TSV.")
parser.add_argument("input_file", help="Path to the input txt file")
parser.add_argument("output_file", help="Path to the output tsv file")
args = parser.parse_args()

# Read file
with open(args.input_file, "r", encoding="utf-8") as f:
    content = f.read().strip()

# Split records by empty lines
records = re.split(r"\n\s*\n", content)

data = []
all_columns = set(["isolate_name"])  # start with isolate_name column

for rec in records:
    lines = rec.strip().split("\n")
    entry = {}

    # --- 1. Isolate name from first line (exclude leading number)
    first_line = lines[0]
    isolate_name = re.sub(r"^\d+:\s*", "", first_line).strip()
    entry["isolate_name"] = isolate_name

    # --- 2. Parse other fields
    for line in lines[1:]:
        line = line.strip()

        # Attributes section
        if line.startswith("/"):
            key, val = line.split("=", 1)
            key = key.strip("/ ").lower()
            val = val.strip('"')
            entry[key] = val
            all_columns.add(key)

        # Identifiers, Organism, etc.
        elif ":" in line and not line.startswith("Accession"):
            field, val = line.split(":", 1)
            field = field.strip()
            val = val.strip()

            # Split on semicolon if multiple
            parts = [p.strip() for p in val.split(";")]
            for p in parts:
                if ":" in p:
                    subfield, subval = p.split(":", 1)
                    subfield = subfield.strip()
                    subval = subval.strip()
                    entry[subfield] = subval
                    all_columns.add(subfield)
                else:
                    entry[field] = val
                    all_columns.add(field)

        # Accession line → Accession + ID
        elif line.startswith("Accession"):
            parts = line.split("\t")
            for p in parts:
                if ":" in p:
                    key, val = p.split(":", 1)
                    entry[key.strip()] = val.strip()
                    all_columns.add(key.strip())

    data.append(entry)

# Build DataFrame
columns_sorted = ["isolate_name"] + sorted(all_columns - {"isolate_name"})
df = pd.DataFrame(data, columns=columns_sorted)

# Save to TSV
df.to_csv(args.output_file, sep="\t", index=False)
