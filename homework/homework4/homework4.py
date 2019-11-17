# CSCI 5481 Homework 4
# University of Minnesota
# Code by Brian Cooper

# Imports
#   argparse: handle command-line arguments
import argparse

# ============================================================================ #
# Command-line argument parser                                                 #
# ============================================================================ #
def make_arg_parser():
    parser = argparse.ArgumentParser(prog='homework4.py',
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
                ids.append(line.replace(">", "").rstrip())
            else:
                # Sequence
                seqs.append(line.rstrip())

    return ids, seqs

# ============================================================================ #
# Main function                                                                #
# ============================================================================ #
if __name__ == "__main__":
    # Handle command line input
    parser = make_arg_parser()
    args = parser.parse_args()

    # Read FASTA sequence file
    ids, seqs = read(args.sequence)
