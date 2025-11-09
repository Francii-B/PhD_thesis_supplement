import os
import sys
import csv
import gzip
import shutil

def COBS_out(directory):
    output_rows = [["Filename", "matched_batch", "matched_genomes", "total_matches"]]
    for filename in os.listdir(directory):
        file_path = os.path.join(directory, filename)
        if os.path.isfile(file_path) and filename.endswith(".gz"):
            try:
                #read file, line by line and count the matches
                with gzip.open(file_path, "rb") as compressed_file:
                    with open("temp_file.txt", "wb") as uncompressed_file:
                        shutil.copyfileobj(compressed_file, uncompressed_file)
                with open("temp_file.txt", "r", newline="", encoding="utf-8") as file:
                    matched_genomes = set()
                    n_matches = 0
                    for line in file:
                        row = line.strip().split("\t")
                        if row[0][0] != "*":
                             matched_genomes.add(row[0].split("_")[1]) #count the number of matched genomes
                             n_matches += 1
                    output_rows.append([filename, int(n_matches > 0), len(matched_genomes), n_matches])
                os.remove("temp_file.txt")
            except Exception as e:
                        print(f"Reading error for '{filename}': {e}")
    return output_rows


def process_tsv_files(directory):
    """
    Read the input based on the file type
    """
    
    #Check the existance of the directory
    if not os.path.isdir(directory):
        print(f"Error: '{directory}' not valid.")
        sys.exit(1)

    #Reading based on file type
    return COBS_out(directory)



def write_output(output_data, output_file):
    """
    Write the output to tsv file.
    """
    try:
        with open(output_file, "w", newline="", encoding="utf-8") as file:
            writer = csv.writer(file, delimiter="\t")
            writer.writerows(output_data)
        print(f"Output file: '{output_file}'")
    except Exception as e:
        print(f"Writing error: {e}")

def main():
    if len(sys.argv) != 3:
        print("Usage: python 581-count_Phylign_matches.py <directory> <output_path_and_name>")
        sys.exit(1)

    directory = sys.argv[1]
    output_name = sys.argv[2]
    output_data = process_tsv_files(directory)
    write_output(output_data, output_name)

if __name__ == "__main__":
    main()
