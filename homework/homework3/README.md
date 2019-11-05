### CSCI 5481 Homework 3
###### Brian Cooper

## Background Information
The file `homework3.py` implements the <a href="">Nei-Saitou neighbor-joining algorithm</a> for phylogeny construction.

## Running the Program
The program can be run as follows:<br>
`python homework3.py -i path_to_sequence_file.fna`

The sequence file must be in FASTA format. A sample is provided in the repository (`hw3.fna`).

## Output
`homework3.py` outputs a file, `distances.txt`, which contains a matrix that describes the pairwise genetic distances between each sequence in the input sequence file. The first row and column contain the sequence identifiers, while the rest of the matrix contains the distances.

Beyond the distances file, the program outputs two more files: `edges.txt`, which describes the edges in the generated phylogeny tree, and `tree.txt`, which contains tree edges in NEWICK format.

## Visualizing
Professor Dan Knights at University of Minnesota prepared two R scripts to visualize the generated tree from the output text files. They are `hw3-plot-edges.r` and `hw3-plot-newick.r` and can be run as follows:
- `Rscript hw3-plot-newick.r tree.txt hw3-tip-labels.txt`
- `Rscript hw3-plot-edges.r edges.txt hw3-tip-labels.txt`

A sample output is in the `tree.pdf` file. Note that these scripts require R and the "ape" and "RColorBrewer" packages.
