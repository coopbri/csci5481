# CSCI 5481 Homework 2
# University of Minnesota
# Code by Brian Cooper
#
# This file implements a modified version (anchored) of the Needleman-Wunsch
#   algorithm for sequence alignment.

import numpy as np

# "constant" score values (not really constant, but functionally)
MATCH = 1
MISMATCH = -3
GAP = -2

# ============================================================================ #
# Create similarity matrix for use by Needleman-Wunsch function                #
# ============================================================================ #
def createMatrix():
    # Initialize similarity matrix (zero matrix)
    V = numpy.zeros(len(query)+1, len(ref)+1)

    # Fill row 1 and column 1 with gap penalty values
    for i in range(len(q)+1):
        V[i][0] = GAP*i
    for j in range(len(r)+1):
        V[0][j] = GAP*j

    # Fill matrix with scores
    for i in range(1, len(q)+1):
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
# def needlemanWunsch(V):


# ============================================================================ #
# Main function                                                                #
# ============================================================================ #
if __name__ == "__main__":
    # Read FASTA sequence files (query and reference)
    # q =
    # r =
