# CSCI 5481 Homework 2
# University of Minnesota
# Code by Brian Cooper
#
# This file implements the Needleman-Wunsch algorithm for sequence alignment. If
#   a match file is specified, a modified version (anchored) is executed.

import argparse
import matplotlib.pyplot as plt
import numpy as np
import string_utils as string
import sys

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
                      help="Path to query sequence (FASTA) [required]")
    parser.add_argument("-r","--reference",
                      default=argparse.SUPPRESS,
                      required=True,
                      help="Path to reference sequence (FASTA) [required]")
    parser.add_argument("-m","--match",
                      default=None,
                      required=False,
                      help="Path to matches file [optional]")
    parser.add_argument("-p","--permute",
                      default=None,
                      required=False,
                      action='store_true',
                      help="Perform a random sequence permutation experiment"
                      + " 10,000 times [optional]")

    return parser

# ============================================================================ #
# Read in FASTA sequence file or match file                                    #
# ============================================================================ #
def read(file, match=False):
    # Match file (for anchored NW)
    if match:
        # Load file (with NumPy for convenience)
        m = np.loadtxt(file)
        # Human sequence indices (first two columns)
        human = m[:, 0:2] - 1
        # Fly sequence indicies (second two columns)
        fly = m[:, 2:4] - 1
        return human, fly

    # FASTA sequence file
    else:
        parsed = ""
        with open(file) as f:
            next(f).rstrip()
            # Parse each line individually
            for line in f:
                parsed += line.rstrip()
        return parsed

# ============================================================================ #
# Create similarity matrix for use by Needleman-Wunsch function                #
#   (param) q: query sequence                                                  #
#   (param) r: reference sequence                                              #
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

            # Score from cell to left
            left = V[i][j-1] + GAP

            # Score from row above
            above = V[i-1][j] + GAP

            # Cell is maximum of diagonal score, left score, and above score
            V[i][j] = max(diag, left, above)

    return V

# ============================================================================ #
# Needleman-Wunsch algorithm                                                   #
#   (param) V: similarity matrix                                               #
#   (param) q: query sequence                                                  #
#   (param) r: reference sequence                                              #
# ============================================================================ #
def needleman_wunsch(V, q, r):
    # Initialize alignments as empty strings
    alignQ = ""
    alignR = ""

    # Length of sequences
    i = len(q)
    j = len(r)

    # Process both sequences entirely
    while i > 0 or j > 0:
        # Determine match or mismatch score
        diff = MATCH if q[i-1] == r[j-1] else MISMATCH

        # Both sequences not fully processed
        if (i > 0 and j > 0):
            # Match or mismatch found
            if V[i-1][j-1] == (V[i][j] - diff):
                alignQ = q[i-1] + alignQ
                alignR = r[j-1] + alignR
                i -= 1
                j -= 1
            # From left
            elif V[i-1][j] == (V[i][j] - GAP):
                alignQ = q[i-1] + alignQ
                # Fill gap on reference
                alignR = "-" + alignR
                i -= 1
            # From above
            elif V[i][j-1] == (V[i][j] - GAP):
                # Fill gap on query
                alignQ = "-" + alignQ
                alignR = r[j-1] + alignR
                j -= 1
        # Query sequence not fully processed
        elif i > 0:
            alignQ = q[i-1] + alignQ
            # Fill gap on reference
            alignR = "-" + alignR
            i -= 1
        # Reference sequence not fully processed
        elif j > 0:
            # Fill gap on query
            alignQ = "-" + alignQ
            alignR = r[j-1] + alignR
            j -= 1

    # Final score, at bottom right of matrix
    finalScore = V[-1][-1]

    return alignQ, alignR, finalScore

# ============================================================================ #
# Anchored Needleman-Wunsch algorithm                                          #
#   (param) q: query sequence                                                  #
#   (param) r: reference sequence                                              #
#   (param) qMatch: query match indices                                        #
#   (param) rMatch: reference match indicies                                   #
# ============================================================================ #
def anchored_nw(q, r, qMatch, rMatch):
    finalScore = 0.0

    alignQ = q[:int(qMatch[0,0])]
    alignR = r[:int(rMatch[0,0])]
    numMatches = np.shape(qMatch)[0]

    # Iterate over number of matches
    for n in range(numMatches):
        # Store indices for each sequence
        qStart = int(qMatch[n][0])
        rStart = int(rMatch[n][0])
        qEnd = int(qMatch[n][1])
        rEnd = int(rMatch[n][1])

        # Generate matrix based on this segment
        V = create_matrix(q[qStart:qEnd], r[rStart:rEnd])

        # Perform alignment on this segment
        _, _, score = needleman_wunsch(V, q[qStart:qEnd], r[rStart:rEnd])

        # Update final score with evaluated segment score
        finalScore += score

        # Update final sequence alignments with processed segments
        alignQ += q[qStart:qEnd]
        alignR += r[rStart:rEnd]

        # Continue to next segment
        if n < numMatches - 1:
            qNext = int(qMatch[n+1][0])
            rNext = int(rMatch[n+1][0])

            # Generate matrix based on this segment
            V = create_matrix(q[qEnd:qNext], r[rEnd:rNext])

            # Perform alignment on this segment
            qSeqSegment, rSeqSegment, score = needleman_wunsch(V,
                q[qEnd:qNext], r[rEnd:rNext])

            # Update final score with evaluated segment score
            finalScore += score

            # Update final sequence alignments with processed segments
            alignQ += qSeqSegment
            alignR += rSeqSegment

        # Reached end of matches, append rest of sequence
        elif n == numMatches - 1:
            alignQ += q[qEnd:]
            alignR += r[rEnd:]

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

    # Match file provided
    if args.match:
        # Process provided match file
        qMatch, rMatch = read(args.match, match=True)
        # Perform Needleman-Wunsch algorithm (anchored)
        alignQ, alignR, score = anchored_nw(q, r, qMatch, rMatch)
    # No match file provided
    else:
        # Perform random permutation experiment
        if args.permute:
            # List of scores for all experiment runs
            scoreList = []

            # Execute experiment 10,000 times
            for i in range(10001):
                # Generate random permutation of amino acids in provided sequences
                qRand = string.shuffle(q)
                rRand = string.shuffle(r)

                # Create similarity matrix based on random permutation
                V = create_matrix(qRand, rRand)

                # Perform Needleman-Wunsch alignment on random permutation
                alignQ, alignR, score = needleman_wunsch(V, qRand, rRand)

                # Add score to list for histogram
                scoreList.append(score)

            # Generate histogram of scores
            plt.title('Score Histogram (Random Permutations)')
            plt.hist(scoreList, 10, histtype="bar", align="mid",
                color="#66ffff", label="Final alignment score",
                edgecolor="black")
            plt.xlabel("Score")
            plt.ylabel("Frequency")
            plt.legend()
            plt.show()
        # Perform regular Needleman-Wunsch algorithm
        else:
            # Create and fill similarity matrix for Needleman-Wunsch
            V = create_matrix(q, r)
            # Perform Needleman-Wunsch algorithm (standard)
            alignQ, alignR, score = needleman_wunsch(V, q, r)

    # Output alignments and score if not doing random permutation experiment
    if not args.permute:
        print("Query alignment:     " + alignQ)
        print("\nReference alignment: " + alignR)
        print("\nFinal score: " + str(score))
