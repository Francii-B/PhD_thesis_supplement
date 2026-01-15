# PhD thesis supplement material
This repository contains all queries and scripts described in my PhD thesis. Prefix numbers correspond to the thesis sections in which these scripts are mentioned. The following list provides a detailed overview of the repository contents:

 - ```query```: Contains the query files used in the experiments.
 - ```COBS/```,```Fulgor/```,```Themisto/```,```Metagraph/```
    - ```53-index```: Snakemake workflow for building the indexes.
    - ```54-query```: Snakemake workflow for performing *k*-mer matching.
 - ```54-collect_search_time.R```: Collects search times from each index type, computes the total CPU time, and generates the corresponding plot (Figure 24).
 - ```56-Pareto_optimization.R```: Performs and plots the Pareto optimization analysis.
 - ```581-cobs_matches-to-tsv.py```: Converts *k*-mer matching output (COBS-like format) into a .TSV file, reporting only query–genome pairs.
 - ```581-count_Phylign_matches-per_batch.py```: Count the number of *k*-mer matching matches and matched genomes for each batch. Takes COBS-like files as input.
 - ```581-Phylign-Fulgor_saturation_curve.R```: Performs the saturation curve analysis and generates the corresponding plot to identify the optimal threshold for Phylign–Fulgor.
 - ```581-common_hits.R```: Computes the set of common hits among Phylign, Phylign–Fulgor, and LexicMap, and generates the corresponding UpSet plots.
 - ```592-Phylign-Fulgor_sam-to-tsv.py```: converts the output .SAM output files from Phylign and Phylign–Fulgor into .TSV format. The .SAM file must not contain blank lines or lines starting with "==>".
 - ```592-tsv-to-presence_matrix.py```: Converts the previously generated .TSV files into a presence/absence matrix.
 - ```592-extract_rows_from_presence_matrix.R```: Extracts rows from the generated matrix corresponding to *E. coli* and UPEC genomes.
 - ```592-plasmid_gene_diffusion-per_genome.R```: Counts the number of genes found in each UPEC genome, generates a boxplot, and performs the Wilcoxon rank-sum test.
 - ```592-plasmid_prevalence_across_UPEC-per_gene.R```: Counts the number of UPEC genomes in which each gene was found.
 - ```5921-biosample_txt-to-tsv.py```: Converts *E. coli* BioSample metadata from .TXT to .TSV.