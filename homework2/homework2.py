# CSCI 5481 Homework 2
# University of Minnesota
# Code by Brian Cooper
#
# This file implements a modified version (anchored) of the Needleman-Wunsch
#   algorithm for sequence alignment.

import numpy as np
import argparse

# "constant" score values (not really constant, but functionally)
MATCH = 1
MISMATCH = -3
GAP = -2

# ============================================================================ #
# Command-line argument parser                                                 #
# ============================================================================ #
def make_arg_parser():
    parser = argparse.ArgumentParser(prog='homework2.py',
                          formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-q","--query",
                      default=argparse.SUPPRESS,
                      required=True,
                      help="Path to query fasta [required]")
    parser.add_argument("-r","--reference",
                      default=argparse.SUPPRESS,
                      required=True,
                      help="Path to reference fasta [required]")
    parser.add_argument("-m","--match",
                      default=None,
                      required=False,
                      help="Path to matches file [optional]")

    return parser

# ============================================================================ #
# Read in FASTA sequence file                                                  #
# ============================================================================ #
def read(file):
	parsed = ""
	with open(file) as f:
		next(f).rstrip()
		for lines in f:
			parsed += lines.rstrip()

	return parsed

# ============================================================================ #
# Create similarity matrix for use by Needleman-Wunsch function                #
# ============================================================================ #
def create_matrix(q, r):
    # Initialize similarity matrix (zero matrix)
    V = np.zeros((len(q)+1, len(r)+1))

    # Fill row 1 and column 1 with gap penalty values
    for i in range(len(q)+1):
        V[i][0] = GAP*i
    for j in range(len(r)+1):
        V[0][j] = GAP*j

    # Fill matrix with scores
    for i in range(1, len(q)+1):
        for j in range(1, len(r)+1):
            # Determine match or mismatch
            if (q[i-1] == r[j-1]):
                diag = V[i-1][j-1] + MATCH
            else:
                diag = V[i-1][j-1] + MISMATCH

            left = V[i][j-1] + GAP
            above = V[i-1][j] + GAP

            V[i][j] = max(diag, left, above)

    return V

# ============================================================================ #
# Needleman-Wunsch algorithm                                                   #
#   (param) V: similarity matrix                                               #
# ============================================================================ #
def needleman_wunsch(V):
    # Initialize alignments as empty strings
    alignQ = ""
    alignR = ""

    # Length of sequences
    i = len(q)
    j = len(r)

    # Process entire sequences
    while i > 0 or j > 0:
        if q[i-1] == r[j-1]:
            match = MATCH
        else:
            match = MISMATCH

        if (i > 0 and j > 0):
            # Match found
            if V[i-1][j-1] == (V[i][j] - match):
                alignQ = q[i-1] + alignQ
                alignR = r[j-1] + alignR
                i -= 1
                j -= 1
            # From left
            elif V[i-1][j] == (V[i][j] - GAP):
                alignQ = q[i-1] + alignQ
                alignR = "-" + alignR
                i -= 1
            # From above
            elif V[i][j-1] == (V[i][j] - GAP):
                alignQ = "-" + alignQ
                alignR = r[j-1] + alignR
                j -= 1
        elif i > 0:
            alignQ = q[i-1] + alignQ
            alignR = "-" + alignR
            i -= 1
        elif j > 0:
            alignQ = "-" + alignQ
            alignR = r[j-1] + alignR
            j -= 1

    # Final score, at bottom right of matrix
    finalScore = V[-1][-1]

    return alignQ, alignR, finalScore

# ============================================================================ #
# Main function                                                                #
# ============================================================================ #
if __name__ == "__main__":
    # Handle command line parameters
    parser = make_arg_parser()
    args = parser.parse_args()

    # Read FASTA sequence files (query and reference)
    q = read(args.query)
    r = read(args.reference)

    # Create and fill similarity matrix for Needleman-Wunsch
    V = create_matrix(q, r)

    # Perform Needleman-Wunsch algorithm
    alignQ, alignR, score = needleman_wunsch(V)

    # Output alignments and score
    print("Query alignment:     " + alignQ)
    print("Reference alignment: " + alignR)
    print("\nFinal score: " + str(score))
