# CSCI 5481 Homework 3
# University of Minnesota
# Code by Brian Cooper
#
# This file implements the Nei-Saitou neighbor-joining algorithm for phylogeny
#   construction. Bootstrap estimation is supported.

import argparse
# import matplotlib.pyplot as plt
import numpy as np
# import string_utils as string
import sys

# ============================================================================ #
# Command-line argument parser                                                 #
# ============================================================================ #
def make_arg_parser():
    parser = argparse.ArgumentParser(prog='homework3.py',
                      formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-i","--sequence",
                      default=argparse.SUPPRESS,
                      required=True,
                      help="Path to sequence (FASTA) file [required]")

    return parser

# ============================================================================ #
# Read in FASTA sequence file                                                  #
# ============================================================================ #
def read(file):
        # Initialize sequence identifer and sequence lists
        ids = []
        seqs = []

        with open(file) as f:
            # Parse each line individually
            for line in f:
                if line.startswith(">"):
                    # Sequence identifier
                    ids.append(line.replace(">","").rstrip())
                else:
                    # Sequence
                    seqs.append(line.rstrip())

        return ids, seqs

# ============================================================================ #
# Calculate genetic distance between two sequences (s1 and s2)                 #
# ============================================================================ #
def distance(s1, s2):
    diff = 0
    length = len(s1)

    for i in range(length):
        if s1[i] != s2[i]:
            diff += 1

    if diff == 0:
        return 0
    else:
        return diff / length

# ============================================================================ #
# Build distance matrix containing pairwise sequence distances                 #
#   Write resultant matrix to file                                             #
# ============================================================================ #
def matrix(ids, seqs):
    length = len(ids) + 1

    # Generate empty matrix
    matrix = [[(0) for x in range(length)] for y in range(length)]
    matrix[0][0] = ""

    # Fill matrix with sequence identifiers
    for i in range(1, length):
        matrix[0][i] = ids[i-1]
    for j in range(1, length):
        matrix[j][0] = ids[j-1]

    # Fill matrix with pairwise sequence distances
    for m in range(1, length):
        for n in range(1, length):
            matrix[m][n] = distance(seqs[m-1], seqs[n-1])

    # Write matrix to file
    with open("distances.txt", "w") as f:
        for row in matrix:
            f.write("\t".join(map(str, row)) + "\t" + "\n")

# ============================================================================ #
# Main function                                                                #
# ============================================================================ #
if __name__ == "__main__":
    # Handle command line parameters
    parser = make_arg_parser()
    args = parser.parse_args()

    # Read FASTA sequence file
    ids, seqs = read(args.sequence)

    # Calculate genetic distance between each pair of sequences in FASTA file
    #   Write resultant distance matrix to file `distances.txt`
    matrix = matrix(ids, seqs)
