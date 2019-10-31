# CSCI 5481 Homework 3
# University of Minnesota
# Code by Brian Cooper
#
# This file implements the Nei-Saitou neighbor-joining algorithm for phylogeny
#   construction. Bootstrap estimation is supported.

import argparse
import matplotlib.pyplot as plt
import numpy as np
import string_utils as string
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


# ============================================================================ #
# Main function                                                                #
# ============================================================================ #
if __name__ == "__main__":
    # Handle command line parameters
    parser = make_arg_parser()
    args = parser.parse_args()

    # Read FASTA sequence file
    q = read(args.sequence)
